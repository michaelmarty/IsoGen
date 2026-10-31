"""Tests for fractional heavy-isotope enrichment of ATOM distributions.

The reference calculations are independent of the native FFT code: a NumPy
FFT reference, brute-force convolution, and (optionally) IsoSpecPy.
"""

import itertools
import re

import numpy as np
import pytest

import isogen
from isogen import isogenwrapper

# (atomic number, nominal shift of the heavy isotope from the lightest isotope)
HEAVY = {"2H": ("H", 1), "13C": ("C", 1), "15N": ("N", 1), "18O": ("O", 2), "34S": ("S", 2)}
NATURAL_LEN = {"H": 2, "C": 2, "N": 2, "O": 3, "S": 5}


def natural_vector(element):
    """Natural per-atom probabilities recovered from the native table."""
    v = isogen.isodist(element, type="ATOM", isolen=NATURAL_LEN[element], dist_only=True)
    return v.astype(float) / v.sum()


def parse_formula(formula):
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def mixed_vector(element, fractions):
    p = natural_vector(element) * (1.0 - sum(f for _, f in fractions))
    p = np.pad(p, (0, max(0, 3 - len(p))))
    for isotope, f in fractions:
        p[HEAVY[isotope][1]] += f
    return p


def reference(formula, enrichment, isolen):
    """Direct polynomial convolution in float64."""
    counts = parse_formula(formula)
    total = np.array([1.0])
    for element, count in counts.items():
        fr = [(iso, f) for iso, f in enrichment.items() if HEAVY[iso][0] == element]
        p = mixed_vector(element, fr)
        for _ in range(count):
            total = np.convolve(total, p)
    out = np.zeros(isolen)
    n = min(isolen, len(total))
    out[:n] = total[:n]
    return out / out.max()


CASES = [
    ("C6H12O6", {"2H": 0.05}),
    ("C6H12O6", {"13C": 0.3}),
    ("C20H30N5O8", {"15N": 0.98}),
    ("C20H30N5O8", {"18O": 0.1}),
    ("C20H30N5O8", {"2H": 0.05, "13C": 0.02, "15N": 0.5, "18O": 0.2}),
    ("H2O", {"2H": 1.0}),
    ("C100H160N30O30S", {"2H": 0.1, "13C": 0.01}),
    ("C50H80N10O10S3", {"34S": 0.25, "13C": 0.5}),
]


def test_none_empty_and_zero_match_natural():
    base = isogen.isodist("C6H12O6", type="ATOM", isolen=64, dist_only=True)
    for options in (None, {}, {"2H": 0.0}, {"2H": 0.0, "13C": 0.0}):
        got = isogen.isodist("C6H12O6", type="ATOM", isolen=64, dist_only=True,
                             heavy_isotopes=options)
        np.testing.assert_array_equal(got, base)


@pytest.mark.parametrize("formula,enrichment", CASES)
def test_matches_numpy_reference(formula, enrichment):
    got = isogen.isodist(formula, type="ATOM", isolen=64, dist_only=True,
                         heavy_isotopes=enrichment)
    np.testing.assert_allclose(got, reference(formula, enrichment, 64), atol=2e-5, rtol=1e-3)


def test_deuterium_alias_and_mass_axis():
    a = isogen.isodist("C6H12O6", type="ATOM", isolen=8, heavy_isotopes={"D": 0.05})
    b = isogen.isodist("C6H12O6", type="ATOM", isolen=8, heavy_isotopes={"2H": 0.05})
    np.testing.assert_array_equal(a, b)
    assert a.shape == (8, 2)
    natural = isogen.isodist("C6H12O6", type="ATOM", isolen=8)
    np.testing.assert_allclose(a[:, 0], natural[:, 0])


def test_full_labeling_is_shifted_distribution():
    # Every H is 2H: envelope of H2 equals the natural envelope of a heavy-only
    # molecule, i.e. one peak at +2 (no M+0/M+1 from hydrogen).
    got = isogen.isodist("H2", type="ATOM", isolen=6, dist_only=True,
                         heavy_isotopes={"2H": 1.0})
    np.testing.assert_allclose(got, [0, 0, 1, 0, 0, 0], atol=1e-6)


def test_single_atom_closed_form():
    f = 0.05
    got = isogen.isodist("C", type="ATOM", isolen=2, dist_only=True, heavy_isotopes={"13C": f})
    p = natural_vector("C")
    m0, m1 = (1 - f) * p[0], (1 - f) * p[1] + f
    np.testing.assert_allclose(got, [1.0, m1 / m0], rtol=1e-5)


def test_binomial_two_hydrogens():
    f = 0.2
    got = isogen.isodist("H2", type="ATOM", isolen=3, dist_only=True, heavy_isotopes={"2H": f})
    p1 = (1 - f) * natural_vector("H")[1] + f
    p0 = (1 - f) * natural_vector("H")[0]
    expected = np.array([p0 ** 2, 2 * p0 * p1, p1 ** 2])
    np.testing.assert_allclose(got, expected / expected.max(), rtol=1e-5)


def test_mean_shift_is_linear_in_fraction():
    isolen = 40

    def mean(f):
        d = isogen.isodist("C10H20N2O4", type="ATOM", isolen=isolen, dist_only=True,
                           heavy_isotopes={"2H": f}).astype(float)
        return (d * np.arange(isolen)).sum() / d.sum()

    m0, m1, m2 = mean(0.0), mean(0.1), mean(0.2)
    assert m1 - m0 == pytest.approx(m2 - m1, rel=1e-3)
    assert m1 - m0 == pytest.approx(20 * 0.1 * (1 - natural_vector("H")[1]), rel=1e-3)


def test_short_isolen_matches_head_of_long_run():
    enrichment = {"13C": 0.4, "2H": 0.1}
    long = isogen.isodist("C200H300N50O60", type="ATOM", isolen=200, dist_only=True,
                          heavy_isotopes=enrichment)
    short = isogen.isodist("C200H300N50O60", type="ATOM", isolen=4, dist_only=True,
                           heavy_isotopes=enrichment)
    # Both are normalized by the full-envelope peak, so truncation must not wrap.
    np.testing.assert_allclose(short, long[:4], rtol=1e-3)
    assert short.max() < 1e-6


def test_offset_shifts_output():
    a = isogenwrapper.fft_gen_atom_isodist("C6H12O6", 8, 0, heavy_isotopes={"2H": 0.1})
    b = isogenwrapper.fft_gen_atom_isodist("C6H12O6", 8, 2, heavy_isotopes={"2H": 0.1})
    np.testing.assert_allclose(b[2:], a[:6], rtol=1e-6)
    assert not b[:2].any()


def test_compat_api_accepts_heavy_isotopes():
    a = isogenwrapper.isogen_atom("C6H12O6", 8, heavy_isotopes={"2H": 0.1})
    b = isogenwrapper.fft_atom_formula_to_dist("C6H12O6", 8, heavy_isotopes={"2H": 0.1})
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("options,error", [
    ({"3H": 0.1}, ValueError),
    ({"XX": 0.1}, ValueError),
    ({"H": 0.1}, ValueError),
    ({"13C": -0.1}, ValueError),
    ({"13C": 1.5}, ValueError),
    ({"13C": float("nan")}, ValueError),
    ({"18O": 0.7, "17O": 0.5}, ValueError),
    ({"15N": 0.1}, ValueError),
    ({"D": 0.1, "2H": 0.1}, ValueError),
    ({"13C": "0.1"}, TypeError),
    ({"13C": True}, TypeError),
    ({13: 0.1}, TypeError),
    ([("13C", 0.1)], TypeError),
])
def test_invalid_heavy_isotopes(options, error):
    with pytest.raises(error):
        isogen.isodist("C6H12O6", type="ATOM", isolen=8, heavy_isotopes=options)


def test_heavy_isotopes_rejected_for_other_inputs():
    with pytest.raises(ValueError):
        isogen.isodist("PEPTIDE", type="PEPTIDE", isolen=8, heavy_isotopes={"2H": 0.1})
    with pytest.raises(ValueError):
        isogen.isodist(1000.0, type="PEPTIDE", isolen=8, heavy_isotopes={"2H": 0.1})
    with pytest.raises(ValueError):
        isogen.isodist("C6H12O6", type="ATOM", method="NN", isolen=8,
                       heavy_isotopes={"2H": 0.1})


def test_brute_force_multinomial_water():
    f_h, f_o = 0.3, 0.2
    ph = mixed_vector("H", [("2H", f_h)])[:2]
    po = mixed_vector("O", [("18O", f_o)])
    expected = np.zeros(8)
    for h1, h2, o in itertools.product(range(2), range(2), range(3)):
        expected[h1 + h2 + o] += ph[h1] * ph[h2] * po[o]
    got = isogen.isodist("H2O", type="ATOM", isolen=8, dist_only=True,
                         heavy_isotopes={"2H": f_h, "18O": f_o})
    np.testing.assert_allclose(got, expected / expected.max(), atol=1e-6, rtol=1e-4)


def isospec_reference(formula, enrichment, isolen):
    IsoSpecPy = pytest.importorskip("IsoSpecPy")
    counts = parse_formula(formula)
    elements = list(counts)
    masses, probs = [], []
    isotope_masses = {"H": [1.00782503207, 2.0141017778], "C": [12.0, 13.0033548378],
                      "N": [14.0030740048, 15.0001088982],
                      "O": [15.99491461956, 16.99913170, 17.9991610],
                      "S": [31.97207117, 32.97145891, 33.96786700, 34.96903, 35.96708071]}
    for element in elements:
        fr = [(iso, f) for iso, f in enrichment.items() if HEAVY[iso][0] == element]
        p = mixed_vector(element, fr)
        n = NATURAL_LEN[element]
        masses.append(isotope_masses[element][:n])
        probs.append(list(p[:n]))
    spec = IsoSpecPy.IsoTotalProb(
        prob_to_cover=0.999999, atomCounts=[counts[e] for e in elements],
        isotopeMasses=masses, isotopeProbabilities=probs)
    mono = sum(counts[e] * masses[i][0] for i, e in enumerate(elements))
    out = np.zeros(isolen)
    for m, p in zip(spec.masses, spec.probs):
        k = int(round(m - mono))
        if k < isolen:
            out[k] += p
    return out / out.max()


@pytest.mark.parametrize("formula,enrichment", [
    ("C6H12O6", {"2H": 0.05}),
    ("C20H30N5O8", {"2H": 0.05, "13C": 0.02, "15N": 0.5, "18O": 0.2}),
    ("C40H60N10O12", {"13C": 0.5}),
])
def test_matches_isospec(formula, enrichment):
    got = isogen.isodist(formula, type="ATOM", isolen=64, dist_only=True,
                         heavy_isotopes=enrichment)
    np.testing.assert_allclose(got, isospec_reference(formula, enrichment, 64),
                               atol=2e-4, rtol=2e-3)

def test_cli_heavy_isotopes(capsys):
    from isogen.__main__ import main
    assert main(["dist", "C6H12O6", "--type", "ATOM", "--isolen", "4",
                 "--heavy-isotopes", "2H=0.05,13C=0.02"]) == 0
    rows = [line.split(",") for line in capsys.readouterr().out.strip().splitlines()[1:]]
    expected = isogen.isodist("C6H12O6", type="ATOM", isolen=4, dist_only=True,
                              heavy_isotopes={"2H": 0.05, "13C": 0.02})
    np.testing.assert_allclose([float(r[1]) for r in rows], expected, rtol=1e-6)


def test_large_counts_match_reference():
    enrichment = {"13C": 0.3, "2H": 0.05}
    got = isogen.isodist("C3000H4500N800O900", type="ATOM", isolen=1500, dist_only=True,
                         heavy_isotopes=enrichment)
    assert np.isfinite(got).all() and got.max() == pytest.approx(1.0)
    mean = (got * np.arange(1500)).sum() / got.sum()
    nat = {el: natural_vector(el) for el in "CHNO"}
    shift = {"C": 0.7 * nat["C"][1] + 0.3, "H": 0.95 * nat["H"][1] + 0.05,
             "N": nat["N"][1], "O": nat["O"][1] + 2 * nat["O"][2]}
    expected = 3000 * shift["C"] + 4500 * shift["H"] + 800 * shift["N"] + 900 * shift["O"]
    assert mean == pytest.approx(expected, rel=2e-3)
