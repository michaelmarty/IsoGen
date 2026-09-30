"""Composition-aware prediction compared with independently specified chemistry."""

import ctypes
import numpy as np
import pytest
from pyteomics import mass as reference_mass

import isogen
from isogen import isogenwrapper as wrapper
from isogen.isogen_tools import peptide_to_dist
from isogen.protein_composition import parse_composition_formula, resolve_proforma_composition


@pytest.mark.parametrize("sequence,formula", [
    ("M[Oxidation]", "C5H11NO3S"),
    ("M[UNIMOD:35]", "C5H11NO3S"),
    ("M[MOD:00719]", "C5H11NO3S"),
    ("M[RESID:AA0581]", "C5H11NO3S"),
    ("[Acetyl]-A", "C5H9NO3"),
    ("A-[Amidated]", "C3H8N2O"),
    ("S[Phospho]", "C3H8NO6P"),
    ("N[Deamidated]", "C4H7NO4"),
    ("C[Carbamidomethyl]", "C5H10N2O3S"),
    ("M[Formula:O]", "C5H11NO3S"),
    ("M[Oxidation|+15.9949|INFO:confirmed]", "C5H11NO3S"),
    ("<[Oxidation]@M>MM", "C10H20N2O5S2"),
    ("[Oxidation]^2?MM", "C10H20N2O5S2"),
    ("{Oxidation}M", "C5H11NO3S"),
    ("(MM)[Oxidation]", "C10H20N2O4S2"),
    ("M[#g1]M[Oxidation#g1]", "C10H20N2O4S2"),
    ("A[Formula:Cl]", "C3H7NO2Cl"),
    ("N[Glycan:HexNAc1Hex2]", "C24H41N3O18"),
    ("S[UNIMOD:41]", "C9H17NO8"),
    ("N[UNIMOD:43]", "C12H21N3O8"),
])
def test_resolved_envelope_matches_independent_formula(sequence, formula, capsys):
    observed = isogen.isodist(sequence, isolen=64)
    expected = isogen.isodist(formula, type="ATOM", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], rtol=3e-5, atol=3e-9)
    assert observed[0, 0] == pytest.approx(expected[0, 0], abs=5e-5)
    assert not capsys.readouterr().err


@pytest.mark.parametrize("method", ["FFT", "BRAIN", "NN"])
def test_default_uses_modifications_and_opt_out_preserves_legacy(method):
    modified = isogen.isodist("M[Oxidation]", method=method, isolen=32)
    ignored = isogen.isodist("M[Oxidation]", method=method, isolen=32, use_modifications=False)
    plain = isogen.isodist("M", method=method, isolen=32)
    assert not np.allclose(modified[:, 1], ignored[:, 1])
    np.testing.assert_array_equal(ignored[:, 1], plain[:, 1])
    np.testing.assert_array_equal(modified[:, 0], ignored[:, 0])


@pytest.mark.parametrize("helper", [wrapper.fft_gen_isodist, wrapper.fft_gen_seq_isodist,
                                    wrapper.nn_gen_isodist, wrapper.nn_gen_seq_isodist,
                                    wrapper.brain_gen_isodist, wrapper.brain_gen_seq_isodist,
                                    wrapper.brain_pep_seq_to_dist, peptide_to_dist])
def test_helpers_use_composition_and_offer_opt_out(helper):
    observed = helper("M[Oxidation]", isolen=32)
    ignored = helper("M[Oxidation]", isolen=32, use_modifications=False)
    assert not np.allclose(observed, ignored)


@pytest.mark.parametrize("tag", ["+15.9949", "Not a modification", "GNO:unknown", "Formula:C-99"])
def test_unsupported_modification_does_not_discard_supported_one(tag, capsys):
    observed = isogen.isodist("M[Oxidation][{}]".format(tag), isolen=32)
    expected = isogen.isodist("M[Oxidation]", isolen=32)
    np.testing.assert_array_equal(observed[:, 1], expected[:, 1])
    assert "Ignoring modification" in capsys.readouterr().err
    if tag == "+15.9949":
        assert observed[0, 0] - expected[0, 0] == pytest.approx(15.9949)


def test_unknown_annotation_is_ignored_on_mass_axis_too(capsys):
    observed = isogen.isodist("M[Unknown]")
    expected = isogen.isodist("M")
    np.testing.assert_array_equal(observed, expected)
    assert capsys.readouterr().err.count("Warning:") == 1
    with pytest.raises(ValueError, match="Unknown"):
        isogen.calc_pep_monoisotopic_mass("M[Unknown]")


def test_opt_out_does_not_warn_for_unresolved_annotation(capsys):
    isogen.isodist("M[+12]", use_modifications=False)
    assert not capsys.readouterr().err


@pytest.mark.parametrize("formula", ["H(-2) 2H(2)", "H -2 (2)H 2", "H-2[2H]2", "H-2[2H2]", "H-2D2"])
def test_fixed_isotope_formula_replaces_natural_atoms(formula):
    sequence = "A[Formula:{}]".format(formula)
    observed = isogen.isodist(sequence, isolen=64)
    # Two fixed deuteriums contribute mass, but no stochastic isotope width.
    expected = isogen.isodist("C3H5NO2", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], atol=3e-9)
    expected_origin = reference_mass.calculate_mass(formula="C3H5NO2") + 2 * reference_mass.nist_mass["H"][2][0]
    assert observed[0, 0] == pytest.approx(expected_origin, abs=5e-5)


def test_global_carbon_label_updates_mass_and_envelope():
    observed = isogen.isodist("<13C>A", isolen=64)
    expected = isogen.isodist("H7NO2", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], atol=3e-9)
    assert observed[0, 0] == pytest.approx(expected[0, 0] + 3 * reference_mass.nist_mass["C"][13][0], abs=5e-5)
    assert isogen.calc_pep_monoisotopic_mass("<13C>A") == pytest.approx(observed[0, 0])
    ignored = isogen.isodist("<13C>A", isolen=64, use_modifications=False)
    np.testing.assert_array_equal(ignored[:, 1], isogen.isodist("A", isolen=64)[:, 1])
    np.testing.assert_array_equal(ignored[:, 0], observed[:, 0])


def test_global_label_includes_modification_atoms():
    observed = isogen.isodist("<13C>[Acetyl]-A", isolen=64)
    expected = isogen.isodist("H9NO3", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], atol=3e-9)
    assert observed[0, 0] == pytest.approx(expected[0, 0] + 5 * reference_mass.nist_mass["C"][13][0], abs=5e-5)


def test_unimod_heavy_label_replaces_carbon_and_nitrogen():
    observed = isogen.isodist("K[UNIMOD:259]", isolen=64)
    expected = isogen.isodist("H14O2", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], atol=3e-9)
    expected_origin = (expected[0, 0] + 6 * reference_mass.nist_mass["C"][13][0]
                       + 2 * reference_mass.nist_mass["N"][15][0])
    assert observed[0, 0] == pytest.approx(expected_origin, abs=5e-5)


def test_deuterium_global_alias():
    np.testing.assert_array_equal(isogen.isodist("<D>A"), isogen.isodist("<2H>A"))


@pytest.mark.parametrize("ion_type", ["a", "b", "c", "y", "z"])
def test_global_hydrogen_label_uses_mass_axis_terminal_composition(ion_type):
    composition = reference_mass.Composition(sequence="AA", ion_type=ion_type)
    hydrogens = composition.pop("H")
    expected = reference_mass.calculate_mass(composition=composition)
    expected += hydrogens * reference_mass.nist_mass["H"][2][0]
    observed = isogen.isodist("<D>AA", ion_type=ion_type)
    assert observed[0, 0] == pytest.approx(expected, abs=5e-5)
    assert isogen.calc_proforma_mass("<D>AA", ion_type=ion_type) == pytest.approx(expected, abs=5e-5)
    expected_average = reference_mass.calculate_mass(composition=composition, average=True)
    expected_average += hydrogens * reference_mass.nist_mass["H"][2][0]
    assert isogen.calc_proforma_mass("<D>AA", ion_type=ion_type, monoisotopic=False) == pytest.approx(expected_average, abs=0.002)


def test_unsupported_isotope_prints_warning_and_keeps_supported_modification(capsys):
    observed = isogen.isodist("M[Oxidation][Formula:[99C2]]")
    expected = isogen.isodist("M[Oxidation]")
    np.testing.assert_array_equal(observed, expected)
    assert "No isotope data" in capsys.readouterr().err
    with pytest.raises(ValueError, match="No isotope data"):
        isogen.calc_pep_monoisotopic_mass("<99C>A")


def test_selenium_mass_origin_matches_elemental_fft():
    observed = isogen.isodist("U", isolen=64)
    expected = isogen.isodist("C3H7NO2Se", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed, expected, atol=5e-5)


def test_cv_selenium_substitution_uses_light_isotope_origin():
    observed = isogen.isodist("C[MOD:00007]", isolen=64)
    expected = isogen.isodist("C3H7NO2Se", type="ATOM", isolen=64)
    np.testing.assert_allclose(observed, expected, atol=5e-5)


def test_mass_options_remain_independent_of_envelope():
    observed = isogen.isodist("EC[Oxidation]", all_cyst_ox=True, pyroglu=True)
    expected_mass = isogen.calc_pep_monoisotopic_mass("EC[Oxidation]", all_cyst_ox=True, pyroglu=True)
    assert observed[0, 0] == pytest.approx(expected_mass, abs=5e-5)


def test_cli_opt_out_preserves_csv_and_warns_only_on_stderr(capsys):
    from isogen.__main__ import main
    assert main(["dist", "M[Unknown]", "--isolen", "4"]) == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("mass,intensity\n")
    assert "Warning:" not in captured.out
    assert "Warning:" in captured.err
    assert main(["dist", "M[Oxidation]", "--isolen", "4", "--ignore-modifications"]) == 0
    captured = capsys.readouterr()
    assert not captured.err


def test_formula_dialects_produce_same_signed_composition():
    assert parse_composition_formula("H(-1) 2H(3) C(2) O") == parse_composition_formula("C 2 H -1 (2)H 3 O 1")


def test_custom_model_keeps_requested_nn_and_prints_warning(capsys):
    from pathlib import Path
    model = Path(wrapper.current_path) / "models" / "isogenpep_model_16.bin"
    observed = isogen.isodist_custom("M[Oxidation]", model, isolen=16)
    expected = isogen.isodist_custom("M", model, isolen=16)
    np.testing.assert_array_equal(observed[:, 1], expected[:, 1])
    assert observed[0, 0] - expected[0, 0] == pytest.approx(15.99491462)
    assert "Custom residue-count NN" in capsys.readouterr().err
    isogen.isodist_custom("M[Oxidation]", model, isolen=16, use_modifications=False)
    assert not capsys.readouterr().err


def test_dist_only_still_prints_warning_once(capsys):
    isogen.isodist("M[+12]", dist_only=True)
    assert capsys.readouterr().err.count("Warning:") == 1


def test_mass_only_tag_follows_successful_composition_resolution(capsys):
    isogen.isodist("M[Oxidation|+12]", isolen=8)
    assert "Conflicting masses" in capsys.readouterr().err


@pytest.mark.parametrize("isolen,offset", [(0, 0), (8, -1), (8, 8)])
def test_modified_output_geometry_is_validated(isolen, offset):
    with pytest.raises(ValueError):
        wrapper.fft_gen_seq_isodist("M[Oxidation]", isolen=isolen, offset=offset)


def test_small_output_truncates_without_fft_wraparound():
    sequence = "A" * 400 + "M[Oxidation]"
    short = isogen.isodist(sequence, isolen=8, dist_only=True)
    long = isogen.isodist(sequence, isolen=128, dist_only=True)
    np.testing.assert_allclose(short, long[:8], atol=1e-8)
    assert short.max() < 1


def test_modified_offset_shifts_distribution():
    direct = wrapper.fft_gen_seq_isodist("M[Oxidation]", isolen=64)
    shifted = wrapper.fft_gen_seq_isodist("M[Oxidation]", isolen=66, offset=2)
    np.testing.assert_array_equal(shifted[:2], 0)
    np.testing.assert_allclose(shifted[2:], direct, atol=1e-8)


def test_heavy_labels_leave_original_natural_fft_unchanged():
    before = wrapper.fft_gen_atom_isodist("C6H12O6", isolen=64)
    isogen.isodist("<13C><15N>A[Formula:H-2[2H]2]", isolen=64)
    np.testing.assert_array_equal(wrapper.fft_gen_atom_isodist("C6H12O6", isolen=64), before)


def test_brain_handles_phosphate_and_falls_back_for_heavy_atoms(capsys):
    observed = isogen.isodist("S[Phospho]", method="BRAIN", isolen=32)
    expected = isogen.isodist("S[Phospho]", method="FFT", isolen=32)
    # The existing BRAIN coefficients use slightly different isotope
    # constants from the general-element FFT table.
    np.testing.assert_allclose(observed[:, 1], expected[:, 1], rtol=0.002, atol=3e-6)
    isogen.isodist("<13C>A", method="BRAIN")
    assert "using FFT" in capsys.readouterr().err


def test_native_counts_reject_negative_atoms():
    counts = np.zeros(109, dtype=np.int32)
    counts[0] = -1
    output = np.zeros(8, dtype=np.float32)
    result = wrapper._fft_atom_counts_to_dist_c(
        counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        None, None, None, 0, output.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 8, 0)
    assert result < 0


def test_residue_composition_does_not_copy_legacy_table_errors():
    resolved = resolve_proforma_composition("ENQS[Formula:O]")
    assert resolved.counts == parse_composition_formula("C17H28N6O11")


def test_charge_spacing_and_terminal_mass_options_still_apply():
    neutral = isogen.isodist("M[Oxidation]", isolen=8, ion_type="b", isotope_spacing=1.01)
    charged = isogen.isodist("M[Oxidation]", isolen=8, ion_type="b", isotope_spacing=1.01, charge=2)
    np.testing.assert_allclose(np.diff(neutral[:, 0]), 1.01)
    np.testing.assert_allclose(charged[:, 0], (neutral[:, 0] + 2 * isogen.mass.mass_proton) / 2)


@pytest.mark.parametrize("sequence", ["M[Oxidation", "M+OTHER", "K[XLMOD:02001#XL1]"])
def test_malformed_or_structurally_unsupported_input_still_raises(sequence):
    with pytest.raises(ValueError):
        isogen.isodist(sequence)
