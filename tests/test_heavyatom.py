"""Tests for calc_heavy_percent: external references and isodist round trips."""

import numpy as np
import pytest

import isogen
from isogen import calc_heavy_percent
from isogen.heavyatom import HeavyPercentResult

# IUPAC-style natural abundances, written independently of the C library.
NATURAL = {
    "H": [0.999885, 0.000115],
    "C": [0.9893, 0.0107],
    "N": [0.99636, 0.00364],
    "O": [0.99757, 0.00038, 0.00205],
}
HEAVY_SHIFT = {"2H": ("H", 1), "13C": ("C", 1), "15N": ("N", 1), "18O": ("O", 2)}
ABS_TOL = 5e-5  # percent; limited by the float32 forward model


def reference_dist(counts, heavy, fraction, length):
    """Independent forward model built from numpy convolutions."""
    element, shift = HEAVY_SHIFT[heavy]
    out = np.array([1.0])
    for symbol, count in counts.items():
        pmf = np.array(NATURAL[symbol], dtype=float)
        if symbol == element:
            pmf = (1 - fraction) * pmf
            pmf = np.pad(pmf, (0, max(0, shift + 1 - pmf.size)))
            pmf[shift] += fraction
        for _ in range(count):
            out = np.convolve(out, pmf)
    out = np.pad(out, (0, max(0, length - out.size)))[:length]
    return out / out.max()


def forward(formula, heavy, fraction, length=40):
    return isogen.isodist(formula, type="ATOM", isolen=length, dist_only=True,
                          heavy_isotopes={heavy: fraction}).astype(float)


def test_independent_reference_matches_isodist_forward_model():
    counts = {"C": 6, "H": 12, "O": 6}
    for heavy, fraction in [("2H", 0.25), ("13C", 0.6), ("18O", 0.1)]:
        np.testing.assert_allclose(
            forward("C6H12O6", heavy, fraction), reference_dist(counts, heavy, fraction, 40), atol=2e-6)


@pytest.mark.parametrize("fraction", [0.0, 0.1, 0.35, 0.7, 1.0])
@pytest.mark.parametrize("heavy", ["2H", "13C", "15N", "18O"])
def test_recovers_independent_reference(heavy, fraction):
    counts = {"C": 10, "H": 14, "N": 5, "O": 7}
    measured = reference_dist(counts, heavy, fraction, 40)
    percent, fit = calc_heavy_percent(measured, "C10H14N5O7", heavy)
    assert percent == pytest.approx(100 * fraction, abs=ABS_TOL)
    assert fit.shape == (40,)


@pytest.mark.parametrize("formula", ["C6H12O6", "C10H14N5O7P", "C40H80NO8P", "C100H160N30O40S2"])
@pytest.mark.parametrize("fraction", [0.0, 0.003, 0.05, 0.5, 0.95, 1.0])
def test_round_trip_with_isodist(formula, fraction):
    for heavy in ("2H", "13C"):
        measured = forward(formula, heavy, fraction, length=300)
        result = calc_heavy_percent(measured, formula, heavy)
        assert result.percent == pytest.approx(100 * fraction, abs=ABS_TOL)
        np.testing.assert_allclose(result.fit_dist, measured, atol=1e-5)
        assert result.rmse < 1e-5


def test_round_trip_with_sulfur_and_phosphorus_formula():
    measured = forward("C20H30N5O8P2S2", "15N", 0.42, length=40)
    assert calc_heavy_percent(measured, "C20H30N5O8P2S2", "15N").percent == pytest.approx(42.0, abs=ABS_TOL)


def test_fit_dist_equals_isodist_at_fitted_percent():
    measured = forward("C6H12O6", "2H", 0.2, length=16)
    result = calc_heavy_percent(measured, "C6H12O6", "2H")
    expected = forward("C6H12O6", "2H", result.percent / 100, length=16)
    np.testing.assert_allclose(result.fit_dist, expected, atol=1e-7)


def test_scale_invariance_and_two_column_input():
    measured = forward("C6H12O6", "13C", 0.3, length=12)
    axis = isogen.isodist("C6H12O6", type="ATOM", isolen=12)[:, 0]
    a = calc_heavy_percent(measured, "C6H12O6", "13C")
    b = calc_heavy_percent(measured * 12345.0, "C6H12O6", "13C")
    c = calc_heavy_percent(np.column_stack([axis, measured * 7.0]), "C6H12O6", "13C")
    assert b.percent == pytest.approx(a.percent, abs=1e-8)
    assert c.percent == pytest.approx(a.percent, abs=1e-8)
    assert b.scale == pytest.approx(12345.0 * a.scale, rel=1e-4)
    assert calc_heavy_percent(list(measured), "C6H12O6", "13C").percent == pytest.approx(a.percent, abs=1e-8)


def test_truncated_window_still_recovers():
    measured = forward("C40H80NO8P", "2H", 0.1, length=40)[:12]
    assert calc_heavy_percent(measured, "C40H80NO8P", "2H").percent == pytest.approx(10.0, abs=ABS_TOL)


def test_aliases_and_bare_symbols():
    measured = forward("C6H12O6", "2H", 0.25, length=14)
    values = [calc_heavy_percent(measured, "C6H12O6", name).percent for name in ("2H", "D", "H")]
    assert values == pytest.approx([25.0] * 3, abs=ABS_TOL)
    assert calc_heavy_percent(forward("C6H12O6", "13C", 0.5), "C6H12O6", "C").percent == pytest.approx(50.0, abs=ABS_TOL)


def test_other_elements_stay_natural():
    # A deuterium fit of a 13C-enriched sample is a poor fit, unlike the correct isotope.
    measured = forward("C6H12O6", "13C", 0.5, length=20)
    wrong = calc_heavy_percent(measured, "C6H12O6", "2H")
    right = calc_heavy_percent(measured, "C6H12O6", "13C")
    assert right.rmse < 1e-5 < wrong.rmse


@pytest.mark.parametrize("weights", ["none", "poisson"])
def test_noisy_data_is_within_error_estimate(weights):
    rng = np.random.default_rng(7)
    truth = 0.08
    clean = forward("C40H80NO8P", "2H", truth, length=30)
    errors = []
    for _ in range(20):
        counts = rng.poisson(clean * 2e5) / 2e5
        result = calc_heavy_percent(counts, "C40H80NO8P", "2H", weights=weights)
        errors.append(result.percent - 100 * truth)
        assert np.isfinite(result.stderr) and result.stderr > 0
    assert abs(np.mean(errors)) < 0.5
    assert np.std(errors) < 1.0


def test_stderr_shrinks_with_less_noise():
    rng = np.random.default_rng(3)
    clean = forward("C6H12O6", "13C", 0.2, length=10)
    noisy = calc_heavy_percent(clean + rng.normal(0, 0.01, clean.size), "C6H12O6", "13C")
    cleaner = calc_heavy_percent(clean + rng.normal(0, 0.0001, clean.size), "C6H12O6", "13C")
    assert cleaner.stderr < noisy.stderr
    assert noisy.percent == pytest.approx(20.0, abs=2.0)


def test_result_object_behaviour():
    measured = forward("C6H12O6", "2H", 0.1, length=10)
    result = calc_heavy_percent(measured, "C6H12O6", "2H")
    assert isinstance(result, HeavyPercentResult)
    percent, fit = result
    assert percent == result.percent and fit is result.fit_dist
    assert result.fraction == pytest.approx(result.percent / 100)
    assert result.heavy_atom == "2H" and result.n_atoms == 12
    assert "HeavyPercentResult" in repr(result)


@pytest.mark.parametrize("dist, formula, atom, error", [
    ([1.0], "C6H12O6", "2H", ValueError),
    ([0.0, 0.0, 0.0], "C6H12O6", "2H", ValueError),
    ([1.0, np.nan, 0.1], "C6H12O6", "2H", ValueError),
    ([[1.0, 2.0, 3.0]] * 3, "C6H12O6", "2H", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", "15N", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", "99C", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", "12C", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", "Xx", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", "2", ValueError),
    ([1.0, 0.5, 0.1], "C6H12O6", 2, TypeError),
    ([1.0, 0.5], "C6H12O6", "18O", ValueError),
    ([1.0, 0.5, 0.1], "not a formula", "2H", ValueError),
])
def test_invalid_inputs(dist, formula, atom, error):
    with pytest.raises(error):
        calc_heavy_percent(dist, formula, atom)


def test_invalid_weights():
    with pytest.raises(ValueError):
        calc_heavy_percent([1.0, 0.5, 0.1], "C6H12O6", "2H", weights="bad")
