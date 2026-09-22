"""Tests for the public ``isogen`` API.

Pyteomics is deliberately used only as a test dependency. It provides an
independent elemental-mass reference without becoming a runtime dependency of
IsoGen.
"""

import ctypes

import numpy as np
import pytest
from pyteomics import mass as pyteomics_mass

import isogen
from isogen import isogenwrapper


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
        "calc_pep_fragments",
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


@pytest.mark.parametrize(
    ("mass_function", "sequence", "cleaned_sequence"),
    [
        (isogen.calc_rna_mass, "AUGX", "AUG"),
        (isogen.calc_rna_monoisotopic_mass, "AUGX", "AUG"),
        (isogen.calc_dna_mass, "ATGX", "ATG"),
        (isogen.calc_dna_monoisotopic_mass, "ATGX", "ATG"),
    ],
)
def test_bad_residue_codes_are_zero_and_quiet_by_default(
    mass_function, sequence, cleaned_sequence, capsys
):
    """Invalid residues should contribute zero without stopping processing."""
    assert mass_function(sequence) == mass_function(cleaned_sequence)
    assert capsys.readouterr().out == ""

    mass_function(sequence, verbose=True)
    assert "Bad" in capsys.readouterr().out


def test_rna_t_as_u_notice_is_opt_in(capsys):
    """RNA accepts thymine silently unless verbose output is requested."""
    observed = isogen.calc_rna_monoisotopic_mass("AT")
    expected = isogen.calc_rna_monoisotopic_mass("AU")
    assert observed == expected
    assert capsys.readouterr().out == ""

    isogen.calc_rna_monoisotopic_mass("AT", verbose=True)
    assert "Assuming T means U" in capsys.readouterr().out


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


def test_calc_pep_fragments_defaults_to_monoisotopic_b_and_y_ions():
    sequence = "PEPTIDE"
    observed = isogen.calc_pep_fragments(sequence)

    assert list(observed) == [
        *(f"b{i}" for i in range(1, len(sequence))),
        *(f"y{i}" for i in range(1, len(sequence))),
    ]
    for ion_name, fragment_mass in observed.items():
        ion_type, length = ion_name[0], int(ion_name[1:])
        fragment = (
            sequence[:length]
            if ion_type in {"a", "b", "c"}
            else sequence[-length:]
        )
        expected = pyteomics_mass.calculate_mass(
            sequence=fragment,
            ion_type=ion_type,
        )
        assert fragment_mass == pytest.approx(expected, abs=3e-5)


def test_calc_pep_fragments_accepts_other_ion_combinations_and_average_mass():
    sequence = "PEPTIDE"
    observed = isogen.calc_pep_fragments(
        sequence,
        ion_types="acz",
        monoisotopic=False,
    )

    assert set(observed) == {
        f"{ion_type}{length}"
        for ion_type in "acz"
        for length in range(1, len(sequence))
    }
    assert observed["a3"] == pytest.approx(
        pyteomics_mass.calculate_mass(
            sequence=sequence[:3],
            ion_type="a",
            average=True,
        ),
        abs=0.02,
    )
    assert observed["z2"] == pytest.approx(
        pyteomics_mass.calculate_mass(
            sequence=sequence[-2:],
            ion_type="z",
            average=True,
        ),
        abs=0.02,
    )


@pytest.mark.parametrize("monoisotopic", [True, False])
def test_calc_pep_fragments_matches_pyteomics(monoisotopic):
    sequence = "PEPTIDE"
    observed = isogen.calc_pep_fragments(
        sequence,
        ion_types="abcxyz",
        monoisotopic=monoisotopic,
    )

    for ion_type in "abcxyz":
        for length in range(1, len(sequence)):
            fragment = (
                sequence[:length]
                if ion_type in "abc"
                else sequence[-length:]
            )
            expected = pyteomics_mass.calculate_mass(
                sequence=fragment,
                ion_type=ion_type,
                average=not monoisotopic,
            )
            tolerance = 3e-5 if monoisotopic else 0.02
            assert observed[f"{ion_type}{length}"] == pytest.approx(
                expected,
                abs=tolerance,
            )


@pytest.mark.parametrize("monoisotopic", [True, False])
def test_z_prime_is_one_hydrogen_heavier_than_pyteomics_z(monoisotopic):
    sequence = "PEPTIDE"
    observed = isogen.calc_pep_fragments(
        sequence,
        ion_types="cz'",
        monoisotopic=monoisotopic,
    )
    expected = pyteomics_mass.calculate_mass(
        sequence=sequence[-3:],
        ion_type="z",
        average=not monoisotopic,
    ) + pyteomics_mass.calculate_mass(
        formula="H", average=not monoisotopic
    )
    tolerance = 3e-5 if monoisotopic else 0.02

    assert "c3" in observed
    assert observed["z'3"] == pytest.approx(expected, abs=tolerance)


@pytest.mark.parametrize("alias", ["z+1", "z•", "z·", "z."])
def test_z_prime_aliases_use_the_same_mass_and_canonical_name(alias):
    expected = isogen.calc_pep_fragments("PEPTIDE", ion_types="z'")
    observed = isogen.calc_pep_fragments("PEPTIDE", ion_types=alias)

    assert observed == expected


@pytest.mark.parametrize(
    ("ion_type", "base_ion_type", "hydrogen_multiplier"),
    [("a+1", "a", 1), ("x+1", "x", 1), ("y-1", "y", -1)],
)
def test_uvpd_hydrogen_shifted_ions_match_pyteomics(
    ion_type, base_ion_type, hydrogen_multiplier
):
    sequence = "PEPTIDE"
    fragment = sequence[:3] if ion_type.startswith("a") else sequence[-3:]
    observed = isogen.calc_pep_monoisotopic_mass(
        fragment, ion_type=ion_type
    )
    expected = pyteomics_mass.calculate_mass(
        sequence=fragment, ion_type=base_ion_type
    ) + hydrogen_multiplier * pyteomics_mass.calculate_mass(formula="H")

    assert observed == pytest.approx(expected, abs=3e-5)


@pytest.mark.parametrize(
    ("fragmentation_type", "ion_types"),
    [
        ("CID", ("b", "y")),
        ("HCD", ("b", "y")),
        ("SID", ("b", "y")),
        ("IRMPD", ("b", "y")),
        ("ETD", ("c", "z'")),
        ("ECD", ("c", "z'")),
        ("EThcD", ("b", "y", "c", "z'")),
        ("BYCZ*", ("b", "y", "c", "z'")),
        ("UVPD", ("a", "b", "c", "x", "y", "z'")),
        ("UVPD4", ("a", "a+1", "x+1", "y-1")),
        ("UVPD6", ("a", "a+1", "x+1", "x", "y-1", "z'")),
        (
            "UVPD9",
            ("a", "a+1", "b", "c", "x", "x+1", "y", "y-1", "z'"),
        ),
    ],
)
def test_fragmentation_type_selects_expected_ion_series(
    fragmentation_type, ion_types
):
    observed = isogen.calc_pep_fragments(
        "PEP", fragmentation_type=fragmentation_type
    )

    def ion_name(ion_type, length):
        if ion_type in {"a+1", "x+1", "y-1"}:
            return f"{ion_type[0]}{length}{ion_type[1:]}"
        return f"{ion_type}{length}"

    assert list(observed) == [
        ion_name(ion_type, length)
        for ion_type in ion_types
        for length in (1, 2)
    ]


def test_explicit_ion_types_override_fragmentation_type():
    observed = isogen.calc_pep_fragments(
        "PEP", ion_types="a", fragmentation_type="ETD"
    )
    assert set(observed) == {"a1", "a2"}


def test_unknown_fragmentation_type_is_rejected():
    with pytest.raises(ValueError, match="fragmentation_type"):
        isogen.calc_pep_fragments("PEP", fragmentation_type="unknown")


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


@pytest.mark.parametrize("method", ["FFT", "NN", "BRAIN"])
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


@pytest.mark.parametrize(
    ("reference_method", "minimum_similarity"),
    [("FFT", 0.99999), ("NN", 0.998)],
)
@pytest.mark.parametrize(
    ("analyte_type", "input_value"),
    [
        ("PEPTIDE", "PEPTIDE"),
        ("PEPTIDE", 10_000.0),
        ("RNA", "AUGCAGUACGUA"),
        ("RNA", 10_000.0),
    ],
)
def test_brain_intensities_agree_with_existing_methods(
    reference_method, minimum_similarity, analyte_type, input_value
):
    """BRAIN should reproduce the established distribution shape."""
    brain = isogen.isodist(
        input_value,
        type=analyte_type,
        isolen=128,
        method="BRAIN",
    )[:, 1]
    reference = isogen.isodist(
        input_value,
        type=analyte_type,
        isolen=128,
        method=reference_method,
    )[:, 1]

    cosine_similarity = np.dot(brain, reference) / (
        np.linalg.norm(brain) * np.linalg.norm(reference)
    )
    assert cosine_similarity > minimum_similarity


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


def test_modified_peptide_nn_uses_full_base_sequence():
    """Bracketed peptide modifications should not truncate the NN sequence."""
    observed = isogenwrapper.nn_gen_seq_isodist(
        "AC[O1]DE",
        type="PEPTIDE",
        isolen=64,
    )
    expected = isogenwrapper.nn_gen_seq_isodist(
        "ACDE",
        type="PEPTIDE",
        isolen=64,
    )
    np.testing.assert_allclose(observed, expected)


def test_modified_peptide_fft_applies_bracketed_formula():
    """FFT peptide intensities should reflect the modification formula."""
    unmodified = isogenwrapper.fft_gen_seq_isodist(
        "ACDE",
        type="PEPTIDE",
        isolen=64,
    )
    modified = isogenwrapper.fft_gen_seq_isodist(
        "AC[O1]DE",
        type="PEPTIDE",
        isolen=64,
    )
    assert not np.allclose(modified, unmodified)


def test_native_peptide_mass_boundary_uses_64_bin_at_11000():
    """The 11 kDa boundary should stay in the 64-length peptide NN model."""
    native = ctypes.CDLL(isogenwrapper.dllpath)
    native.nn_pep_mass_to_isolen.argtypes = [ctypes.c_float]
    native.nn_pep_mass_to_isolen.restype = ctypes.c_int

    assert native.nn_pep_mass_to_isolen(ctypes.c_float(11000.0)) == 64


@pytest.mark.parametrize(
    ("function_name", "sequence", "tail_start"),
    [
        ("nn_pep_seq_to_dist", b"AC", 16),
        ("nn_rna_seq_to_dist", b"AUGC", 64),
    ],
)
def test_native_nn_sequence_outputs_zero_unused_tail(
    function_name,
    sequence,
    tail_start,
):
    """Direct C callers should receive zeroed output beyond the model length."""
    native = ctypes.CDLL(isogenwrapper.dllpath)
    function = getattr(native, function_name)
    function.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int,
        ctypes.c_int,
    ]
    function.restype = ctypes.c_float

    output = (ctypes.c_float * 128)(*([7.0] * 128))
    result = function(sequence, output, 128, 0)
    values = np.ctypeslib.as_array(output)

    assert result > 0
    np.testing.assert_allclose(values[tail_start:], 0.0)


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


@pytest.mark.parametrize(
    ("polarity", "proton_sign"),
    [("positive", 1), ("negative", -1)],
)
@pytest.mark.parametrize("charge", [1, 2])
def test_isodist_charge_returns_charge_adjusted_mz_axis(
    charge, polarity, proton_sign
):
    """Public distributions should convert neutral axes to the requested m/z."""
    neutral_mass = 10_000.0
    distribution = isogen.isodist(
        neutral_mass,
        type="PEPTIDE",
        isolen=8,
        isotope_spacing=1.01,
        charge=charge,
        polarity=polarity,
    )

    expected_origin = (
        neutral_mass + proton_sign * charge * 1.00727647
    ) / charge
    assert distribution[0, 0] == pytest.approx(expected_origin)
    np.testing.assert_allclose(np.diff(distribution[:, 0]), 1.01 / charge)


@pytest.mark.parametrize("charge", [None, 0, -1])
def test_isodist_charge_less_than_one_preserves_neutral_axis(charge):
    """None, zero, and negative charges should retain the neutral mass axis."""
    neutral_mass = 10_000.0
    distribution = isogen.isodist(
        neutral_mass,
        type="PEPTIDE",
        isolen=8,
        isotope_spacing=1.01,
        charge=charge,
    )

    assert distribution[0, 0] == neutral_mass
    np.testing.assert_allclose(np.diff(distribution[:, 0]), 1.01)


def test_dna_uses_rna_intensities_but_dna_mass_axis():
    """DNA's documented RNA intensity approximation should stay explicit."""
    dna = isogen.isodist("ATGC", type="DNA", isolen=16)
    rna = isogen.isodist("AUGC", type="RNA", isolen=16)

    np.testing.assert_allclose(dna[:, 1], rna[:, 1])
    assert dna[0, 0] != rna[0, 0]
