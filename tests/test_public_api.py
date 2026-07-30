"""Tests for the public ``isogen`` API.

Pyteomics is deliberately used only as a test dependency. It provides an
independent elemental-mass reference without becoming a runtime dependency of
IsoGen.
"""

import numpy as np
import pytest
from pyteomics import mass as pyteomics_mass

import isogen


RNA_RESIDUE_FORMULAS = {
    "A": "C10H12N5O6P",
    "C": "C9H12N3O7P",
    "G": "C10H12N5O7P",
    "U": "C9H11N2O8P",
}

DNA_RESIDUE_FORMULAS = {
    "A": "C10H12N5O5P",
    "C": "C9H12N3O6P",
    "G": "C10H12N5O6P",
    "T": "C10H13N2O7P",
}

LIGHT_ISOTOPES = {
    "H": 1,
    "C": 12,
    "N": 14,
    "O": 16,
    "P": 31,
    "S": 32,
}


def pyteomics_oligo_mass(sequence, residue_formulas, average=False):
    """Calculate the default 3'-OH/5'-monophosphate oligo mass."""
    composition = pyteomics_mass.Composition(formula="H2O")
    for residue in sequence:
        composition += pyteomics_mass.Composition(
            formula=residue_formulas[residue]
        )
    return pyteomics_mass.calculate_mass(
        composition=composition,
        average=average,
    )


def pyteomics_nominal_distribution(sequence, isolen):
    """Bin Pyteomics isotopologues by their nominal neutron shift."""
    distribution = np.zeros(isolen)
    isotopologues = pyteomics_mass.isotopologues(
        sequence=sequence,
        report_abundance=True,
        isotope_threshold=1e-5,
        overall_threshold=1e-12,
    )
    for composition, abundance in isotopologues:
        isotope_number = 0
        for isotope_name, count in composition.items():
            element, mass_number = isotope_name.rstrip("]").split("[")
            isotope_number += (
                int(mass_number) - LIGHT_ISOTOPES[element]
            ) * count
        if isotope_number < isolen:
            distribution[isotope_number] += abundance
    return distribution / distribution.max()


def test_public_exports_are_available():
    """Every documented top-level function should be importable."""
    expected = {
        "isodist",
        "calc_atom_mass_axis",
        "calc_atom_monoisotopic_mass",
        "calc_dna_mass",
        "calc_dna_mass_axis",
        "calc_dna_monoisotopic_mass",
        "calc_mass_axis",
        "calc_pep_mass",
        "calc_pep_mass_axis",
        "calc_pep_monoisotopic_mass",
        "calc_rna_mass",
        "calc_rna_mass_axis",
        "calc_rna_monoisotopic_mass",
        "gen_mass_axis",
    }
    assert expected <= set(isogen.__all__)
    assert all(callable(getattr(isogen, name)) for name in expected)


def test_calc_mass_axis_length_origin_and_spacing():
    """The generic axis should preserve its requested geometry."""
    axis = isogen.calc_mass_axis(1234.5, isolen=7, isotope_spacing=1.01)
    assert axis.shape == (7,)
    assert axis[0] == pytest.approx(1234.5)
    np.testing.assert_allclose(np.diff(axis), 1.01)


def test_gen_mass_axis_dispatch_and_spacing_override():
    """The generic dispatcher should expose sequence-specific mass options."""
    axis = isogen.gen_mass_axis(
        "PEPTID",
        type="PEPTIDE",
        isolen=5,
        isotope_spacing=1.01,
        ion_type="b",
    )
    expected = pyteomics_mass.calculate_mass(
        sequence="PEPTID",
        ion_type="b",
    )
    assert axis[0] == pytest.approx(expected, abs=3e-5)
    np.testing.assert_allclose(np.diff(axis), 1.01)


def test_peptide_masses_match_pyteomics():
    """Intact peptide masses should agree with the Pyteomics reference."""
    sequence = "PEPTIDE"
    reference_mono = pyteomics_mass.calculate_mass(sequence=sequence)
    reference_average = pyteomics_mass.calculate_mass(
        sequence=sequence,
        average=True,
    )

    # IsoGen residue monoisotopic masses are stored to five decimal places.
    assert isogen.calc_pep_monoisotopic_mass(sequence) == pytest.approx(
        reference_mono,
        abs=3e-5,
    )
    # The legacy average residue table is lower precision (four decimals).
    assert isogen.calc_pep_mass(sequence, round_to=6) == pytest.approx(
        reference_average,
        abs=0.02,
    )


@pytest.mark.parametrize("ion_type", ["a", "b", "c", "x", "y", "z"])
def test_peptide_fragment_masses_match_pyteomics(ion_type):
    """Neutral a/b/c/x/y/z fragment masses should match Pyteomics."""
    full_sequence = "PEPTIDE"
    fragment = (
        full_sequence[:-1]
        if ion_type in {"a", "b", "c"}
        else full_sequence[1:]
    )
    expected = pyteomics_mass.calculate_mass(
        sequence=fragment,
        ion_type=ion_type,
    )
    observed = isogen.calc_pep_monoisotopic_mass(
        fragment,
        ion_type=ion_type,
    )
    assert observed == pytest.approx(expected, abs=3e-5)


@pytest.mark.parametrize(
    ("kind", "sequence", "formulas"),
    [
        ("rna", "AUGC", RNA_RESIDUE_FORMULAS),
        ("dna", "ATGC", DNA_RESIDUE_FORMULAS),
    ],
)
def test_nucleic_acid_masses_match_pyteomics(kind, sequence, formulas):
    """Default oligonucleotide masses should match elemental formulas."""
    expected_mono = pyteomics_oligo_mass(sequence, formulas)
    expected_average = pyteomics_oligo_mass(
        sequence,
        formulas,
        average=True,
    )
    observed_mono = getattr(
        isogen,
        f"calc_{kind}_monoisotopic_mass",
    )(sequence)
    observed_average = getattr(isogen, f"calc_{kind}_mass")(sequence)

    assert observed_mono == pytest.approx(expected_mono, abs=1e-5)
    # IsoGen's average nucleotide table is intentionally one-decimal input.
    assert observed_average == pytest.approx(expected_average, abs=0.05)


@pytest.mark.parametrize(
    ("kind", "sequence", "formulas", "spacing"),
    [
        ("pep", "PEPTIDE", None, 1.0033),
        ("rna", "AUGC", RNA_RESIDUE_FORMULAS, 1.0027),
        ("dna", "ATGC", DNA_RESIDUE_FORMULAS, 1.0027),
    ],
)
def test_sequence_mass_axes_start_at_pyteomics_mass(
    kind,
    sequence,
    formulas,
    spacing,
):
    """Each sequence axis should use a reference monoisotopic origin."""
    if kind == "pep":
        expected_mass = pyteomics_mass.calculate_mass(sequence=sequence)
    else:
        expected_mass = pyteomics_oligo_mass(sequence, formulas)

    axis = getattr(isogen, f"calc_{kind}_mass_axis")(sequence, isolen=6)
    assert axis.shape == (6,)
    assert axis[0] == pytest.approx(expected_mass, abs=3e-5)
    np.testing.assert_allclose(np.diff(axis), spacing)


@pytest.mark.parametrize("method", ["FFT", "NN"])
@pytest.mark.parametrize(
    ("analyte_type", "sequence", "expected_mass"),
    [
        (
            "PEPTIDE",
            "PEPTIDE",
            lambda: pyteomics_mass.calculate_mass(sequence="PEPTIDE"),
        ),
        (
            "RNA",
            "AUGC",
            lambda: pyteomics_oligo_mass("AUGC", RNA_RESIDUE_FORMULAS),
        ),
        (
            "DNA",
            "ATGC",
            lambda: pyteomics_oligo_mass("ATGC", DNA_RESIDUE_FORMULAS),
        ),
    ],
)
def test_isodist_sequence_outputs(method, analyte_type, sequence, expected_mass):
    """Top-level sequence distributions should have valid reference axes."""
    distribution = isogen.isodist(
        sequence,
        type=analyte_type,
        isolen=16,
        method=method,
    )

    assert distribution.shape == (16, 2)
    assert np.isfinite(distribution).all()
    assert distribution[0, 0] == pytest.approx(expected_mass(), abs=3e-5)
    assert distribution[:, 1].min() >= 0
    assert distribution[:, 1].max() == pytest.approx(1.0, abs=1e-6)


def test_fft_peptide_intensities_match_pyteomics_isotopologues():
    """FFT intensities should match a Pyteomics isotope enumeration."""
    sequence = "AG"
    observed = isogen.isodist(
        sequence,
        type="PEPTIDE",
        isolen=5,
        method="FFT",
    )[:, 1]
    expected = pyteomics_nominal_distribution(sequence, isolen=5)
    np.testing.assert_allclose(observed, expected, rtol=2e-5, atol=1e-8)


def test_isodist_forwards_fragment_mass_options():
    """Mass kwargs should reach the peptide mass-axis calculation."""
    fragment = "PEPTID"
    distribution = isogen.isodist(
        fragment,
        type="PEPTIDE",
        isolen=8,
        ion_type="b",
    )
    expected = pyteomics_mass.calculate_mass(
        sequence=fragment,
        ion_type="b",
    )
    assert distribution[0, 0] == pytest.approx(expected, abs=3e-5)


def test_numeric_isodist_uses_input_as_axis_origin():
    """Numeric input should be treated as the supplied monoisotopic mass."""
    distribution = isogen.isodist(
        10_000,
        type="PEPTIDE",
        isolen=8,
        isotope_spacing=1.01,
    )
    assert distribution.shape == (8, 2)
    assert distribution[0, 0] == 10_000
    np.testing.assert_allclose(np.diff(distribution[:, 0]), 1.01)


def test_dna_uses_rna_intensities_but_dna_mass_axis():
    """DNA's documented RNA intensity approximation should stay explicit."""
    dna = isogen.isodist("ATGC", type="DNA", isolen=16)
    rna = isogen.isodist("AUGC", type="RNA", isolen=16)

    np.testing.assert_allclose(dna[:, 1], rna[:, 1])
    assert dna[0, 0] != rna[0, 0]
