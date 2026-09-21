"""Tests for ProForma modified-protein mass calculations."""

from pathlib import Path

import numpy as np
import pytest
from pyteomics import mass as pyteomics_mass

import isogen


MONO_TOLERANCE = 4e-5


@pytest.mark.parametrize(
    "proforma",
    [
        "EM[Oxidation]E",
        "EM[U:Oxidation]E",
        "EM[UNIMOD:35]E",
        "EM[MOD:00719]E",
        "EM[M:L-methionine sulfoxide]E",
    ],
)
def test_unimod_and_psi_mod_match_pyteomics(proforma):
    # Pyteomics' ProForma resolvers require optional lxml/psims packages. Use
    # Pyteomics for the unmodified reference mass and add the CV-defined delta.
    expected = pyteomics_mass.calculate_mass(sequence="EME") + 15.994915
    assert isogen.calc_pep_monoisotopic_mass(proforma) == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )


@pytest.mark.parametrize(
    "proforma",
    [
        "EM[RESID:AA0581]E",
        "EM[R:L-methionine (R)-sulfoxide]E",
    ],
)
def test_resid_origin_specific_correction(proforma):
    expected = pyteomics_mass.calculate_mass(sequence="EME") + 15.994915
    assert isogen.calc_pep_monoisotopic_mass(proforma) == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )


def test_resid_rejects_the_wrong_origin_residue():
    with pytest.raises(ValueError, match="not defined for site"):
        isogen.calc_pep_monoisotopic_mass("ES[RESID:AA0581]E")


@pytest.mark.parametrize("residue", ["J", "O", "U"])
def test_unusual_residues_match_pyteomics(residue):
    expected = pyteomics_mass.calculate_mass(sequence="A" + residue)
    observed = isogen.calc_pep_monoisotopic_mass("A" + residue)
    assert observed == pytest.approx(expected, abs=MONO_TOLERANCE)


@pytest.mark.parametrize("residue", ["J", "O", "U"])
def test_unusual_residue_average_masses_match_pyteomics(residue):
    expected = pyteomics_mass.calculate_mass(sequence="A" + residue, average=True)
    observed = isogen.calc_pep_mass("A" + residue, round_to=6)
    assert observed == pytest.approx(expected, abs=0.02)


def test_b_and_z_use_the_documented_ambiguity_midpoint():
    alternatives = [
        pyteomics_mass.calculate_mass(sequence="A" + first + second)
        for first in "DN"
        for second in "EQ"
    ]
    assert isogen.calc_pep_monoisotopic_mass("ABZ") == pytest.approx(
        np.mean(alternatives), abs=MONO_TOLERANCE
    )


def test_x_is_a_zero_base_known_mass_gap():
    expected = pyteomics_mass.calculate_mass(sequence="RTAAWT") + 367.0537
    observed = isogen.calc_pep_monoisotopic_mass("RTAAX[+367.0537]WT")
    assert observed == pytest.approx(expected, abs=MONO_TOLERANCE)


def test_bare_x_is_rejected():
    with pytest.raises(ValueError, match="explicit signed mass gap"):
        isogen.calc_pep_monoisotopic_mass("PEPXIDE")


@pytest.mark.parametrize(
    ("proforma", "plain", "delta"),
    [
        ("PEP[+15.9949]TIDE", "PEPTIDE", 15.9949),
        ("SEQUEN[Formula:C2H2O]CE", "SEQUENCE", 42.010564684),
        ("[Acetyl]-PEPTIDE", "PEPTIDE", 42.010565),
        ("PEPTIDE-[Amidated]", "PEPTIDE", -0.984016),
        ("{Oxidation}PEPTIDE", "PEPTIDE", 15.994915),
        ("[Phospho]?PEPTIDE", "PEPTIDE", 79.966331),
        ("[Phospho]^2?PEPTIDE", "PEPTIDE", 2 * 79.966331),
        ("PROT(EOSFORMS)[+19.0523]ISK", "PROTEOSFORMSISK", 19.0523),
    ],
)
def test_supported_proforma_locations_and_explicit_deltas(
    proforma, plain, delta
):
    expected = pyteomics_mass.calculate_mass(sequence=plain) + delta
    observed = isogen.calc_pep_monoisotopic_mass(proforma)
    assert observed == pytest.approx(expected, abs=MONO_TOLERANCE)


def test_global_fixed_modification_applies_to_each_selected_residue():
    proforma = "<[Carbamidomethyl]@C>ACDC"
    expected = pyteomics_mass.calculate_mass(sequence="ACDC") + 2 * 57.021464
    assert isogen.calc_pep_monoisotopic_mass(proforma) == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )


def test_modified_average_mass_uses_the_ontology_average_delta():
    # Pyteomics treats ProForma deltas as monoisotopic even with average=True,
    # so validate the average base sequence independently and add UniMod's
    # average oxidation delta.
    expected = pyteomics_mass.calculate_mass(sequence="EME", average=True) + 15.9994
    observed = isogen.calc_pep_mass("EM[Oxidation]E", round_to=6)
    assert observed == pytest.approx(expected, abs=0.02)


def test_ambiguous_localization_group_contributes_one_modification():
    proforma = "EMEVT[#g1]S[#g1]ES[Phospho#g1]PEK"
    plain = "EMEVTSESpeK".upper()
    expected = pyteomics_mass.calculate_mass(sequence=plain) + 79.966331
    assert isogen.calc_pep_monoisotopic_mass(proforma) == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )


def test_pipe_descriptors_are_checked_and_info_is_ignored():
    proforma = "EM[Oxidation|+15.9949|INFO:confirmed]E"
    expected = pyteomics_mass.calculate_mass(sequence="EME") + 15.994915
    assert isogen.calc_pep_monoisotopic_mass(proforma) == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )
    with pytest.raises(ValueError, match="Conflicting masses"):
        isogen.calc_pep_monoisotopic_mass("EM[Oxidation|+12.0]E")


def test_neutral_charge_suffix_does_not_change_mass():
    assert isogen.calc_pep_monoisotopic_mass(
        "VAEINPSNGGTT/2"
    ) == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="VAEINPSNGGTT"),
        abs=MONO_TOLERANCE,
    )


@pytest.mark.parametrize(
    ("proforma", "message"),
    [
        ("PEP[Not a modification]TIDE", "Unknown"),
        ("PEP[GNO:G62765YT]TIDE", "not yet supported"),
        ("PEPK[XLMOD:02001#XL1]TIDE", "not supported"),
        ("<13C>PEPTIDE", "not yet supported"),
        ("PEPTIDE+OTHER", "Chimeric"),
        ("PEP[Oxidation", "Unbalanced"),
    ],
)
def test_unsupported_or_malformed_constructs_raise(proforma, message):
    with pytest.raises(ValueError, match=message):
        isogen.calc_pep_monoisotopic_mass(proforma)


def test_strip_proforma_handles_all_supported_annotation_locations():
    value = "<[Carbamidomethyl]@C>{Oxidation}[Phospho]?-AC[+1]D-[Amidated]"
    # The leading '-' after an unlocalized modification is malformed, while a
    # conventional N-terminal tag is accepted.
    with pytest.raises(ValueError):
        isogen.strip_proforma(value)
    assert isogen.strip_proforma(
        "<[Carbamidomethyl]@C>{Oxidation}[Acetyl]-AC[+1]D-[Amidated]/2"
    ) == "ACD"


@pytest.mark.parametrize("method", ["FFT", "NN", "BRAIN"])
def test_isodist_strips_modifications_but_uses_modified_mass_axis(method):
    modified = isogen.isodist(
        "EM[Oxidation]E", type="PEPTIDE", method=method, isolen=16
    )
    plain = isogen.isodist("EME", type="PEPTIDE", method=method, isolen=16)
    np.testing.assert_allclose(modified[:, 1], plain[:, 1])
    assert modified[0, 0] == pytest.approx(
        plain[0, 0] + 15.994915, abs=MONO_TOLERANCE
    )


def test_distribution_only_still_validates_unsupported_modifications():
    with pytest.raises(ValueError, match="not yet supported"):
        isogen.isodist(
            "PEP[GNO:G62765YT]TIDE",
            type="PEPTIDE",
            method="FFT",
            isolen=16,
            dist_only=True,
        )


def test_plain_sequence_regression_and_public_exports():
    expected = pyteomics_mass.calculate_mass(sequence="PEPTIDE")
    assert isogen.calc_pep_monoisotopic_mass("PEPTIDE") == pytest.approx(
        expected, abs=MONO_TOLERANCE
    )
    assert "calc_proforma_mass" in isogen.__all__
    assert "strip_proforma" in isogen.__all__


def test_modification_resource_is_present_and_compact():
    resource = Path(isogen.__file__).with_name("resources") / "protein_modifications.json.gz"
    assert resource.is_file()
    assert resource.stat().st_size < 250_000
