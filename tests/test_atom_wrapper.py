"""Tests for the elemental-formula native wrappers.

Pyteomics is used only as an independent test reference and remains an
optional test dependency.
"""

import numpy as np
import pytest
from pyteomics import mass as pyteomics_mass

import isogen
from isogen import isogenwrapper


LIGHT_ISOTOPES = {
    "H": 1,
    "C": 12,
    "N": 14,
    "O": 16,
}


def pyteomics_formula_distribution(formula, isolen):
    """Bin formula isotopologues by nominal neutron shift."""
    distribution = np.zeros(isolen)
    composition = pyteomics_mass.Composition(formula=formula)
    isotopologues = pyteomics_mass.isotopologues(
        composition=composition,
        report_abundance=True,
        isotope_threshold=1e-7,
        overall_threshold=1e-14,
    )
    for isotope_composition, abundance in isotopologues:
        isotope_number = 0
        for isotope_name, count in isotope_composition.items():
            element, mass_number = isotope_name.rstrip("]").split("[")
            isotope_number += (
                int(mass_number) - LIGHT_ISOTOPES[element]
            ) * count
        if isotope_number < isolen:
            distribution[isotope_number] += abundance
    return distribution / distribution.max()


def test_atom_formula_to_vector_uses_atomic_number_order():
    """Formula counts should occupy their atomic-number-minus-one indices."""
    counts = isogenwrapper.atom_formula_to_vector("C6H12O6")

    assert counts.dtype == np.int32
    assert counts.shape == (isogenwrapper.ATOM_ELEMENT_COUNT,)
    assert counts[0] == 12
    assert counts[5] == 6
    assert counts[7] == 6
    assert counts.sum() == 24


def test_atom_fft_distribution_matches_pyteomics():
    """Native formula intensities should match independent enumeration."""
    observed = isogenwrapper.fft_gen_atom_isodist("C2H5NO", isolen=8)
    expected = pyteomics_formula_distribution("C2H5NO", isolen=8)

    assert observed.dtype == np.float32
    np.testing.assert_allclose(observed, expected, rtol=2e-5, atol=1e-11)


def test_public_atom_isodist_has_formula_mass_axis():
    """ATOM should propagate through the public mass/intensity API."""
    formula = "C6H12O6"
    observed = isogen.isodist(
        formula,
        type="ATOM",
        isolen=8,
        method="FFT",
    )
    expected_mass = pyteomics_mass.calculate_mass(formula=formula)

    assert observed.shape == (8, 2)
    assert observed[0, 0] == pytest.approx(expected_mass, abs=1e-8)
    np.testing.assert_allclose(np.diff(observed[:, 0]), 1.0033)
    assert observed[:, 1].max() == pytest.approx(1.0, abs=1e-6)


def test_public_atom_type_is_case_insensitive_and_spacing_is_editable():
    """Public ATOM dispatch should normalize names and forward mass options."""
    observed = isogen.isodist(
        "H2O",
        type="atom",
        isolen=4,
        isotope_spacing=1.01,
    )

    np.testing.assert_allclose(np.diff(observed[:, 0]), 1.01)


def test_public_atom_rejects_neural_network_method():
    """Elemental formulas have no neural-network model."""
    with pytest.raises(ValueError, match="only the FFT"):
        isogen.isodist("H2O", type="ATOM", method="NN")


def test_atom_fft_offset_and_generic_dispatch():
    """Offsets and FORMULA dispatch should preserve the native distribution."""
    direct = isogenwrapper.fft_gen_atom_isodist("C2H5NO", isolen=8)
    offset = isogenwrapper.fft_atom_formula_to_dist(
        "C2H5NO", isolen=10, offset=2
    )
    dispatched = isogenwrapper.fft_gen_isodist(
        "C2H5NO", type="FORMULA", isolen=8
    )
    compatibility = isogenwrapper.isogen_atom("C2H5NO", isolen=8)

    np.testing.assert_array_equal(offset[:2], 0)
    np.testing.assert_allclose(offset[2:], direct, rtol=2e-5, atol=1e-12)
    np.testing.assert_allclose(dispatched, direct)
    np.testing.assert_allclose(compatibility, direct)


@pytest.mark.parametrize("formula", ["", "C6H12O6?", "NotAnElement"])
def test_atom_formula_rejects_invalid_input(formula):
    """Malformed or unsupported formulas should raise a Python error."""
    with pytest.raises(ValueError, match="Invalid elemental formula"):
        isogenwrapper.atom_formula_to_vector(formula)
    with pytest.raises(ValueError, match="Invalid elemental formula"):
        isogenwrapper.fft_gen_atom_isodist(formula)


@pytest.mark.parametrize(
    ("isolen", "offset"),
    [(0, 0), (8, -1), (8, 8)],
)
def test_atom_fft_rejects_invalid_output_geometry(isolen, offset):
    """Invalid vector lengths and offsets should fail before entering C."""
    with pytest.raises(ValueError):
        isogenwrapper.fft_gen_atom_isodist(
            "H2O", isolen=isolen, offset=offset
        )
