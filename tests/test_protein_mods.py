"""Tests for ProForma modified-protein mass calculations."""

from pathlib import Path

import numpy as np
import pytest
from pyteomics import mass as pyteomics_mass

import isogen


MONO_TOLERANCE = 4e-5
CA_SEQUENCE = (
    "S[Acetylation]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEAT"
    "SRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLV"
    "HWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNV"
    "LDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNR"
    "QVRGFPK"
)
CA_PLAIN_SEQUENCE_LENGTH = 259
CA_EXPERIMENTAL_MASSES = tuple(
    float(value)
    for value in (
        Path(__file__).with_name("data") / "ca_etd_masses.txt"
    ).read_text(encoding="utf-8-sig").split()
)


def test_ca_sequence_fragment_coverage_is_60_to_61_percent():
    predicted_fragments = isogen.calc_pep_fragments(
        CA_SEQUENCE, ion_types=("c", "z'")
    )
    match_candidates = sorted(
        (
            abs(experimental - predicted) / predicted * 1e6,
            ion,
            predicted,
            experimental,
        )
        for ion, predicted in predicted_fragments.items()
        for experimental in CA_EXPERIMENTAL_MASSES
        if abs(experimental - predicted) / predicted * 1e6 <= 20
    )

    matches = []
    matched_ions = set()
    matched_features = set()
    for _, ion, predicted, experimental in match_candidates:
        if ion in matched_ions or experimental in matched_features:
            continue
        matched_ions.add(ion)
        matched_features.add(experimental)
        matches.append((ion, predicted, experimental))

    cleavage_sites = {
        int(ion.lstrip("abcxyz'").partition("#")[0])
        if ion.startswith("c")
        else CA_PLAIN_SEQUENCE_LENGTH
        - int(ion.lstrip("abcxyz'").partition("#")[0])
        for ion, _, _ in matches
    }
    coverage = len(cleavage_sites) / (CA_PLAIN_SEQUENCE_LENGTH - 1) * 100

    assert len(matches) == 159
    assert len(cleavage_sites) == 155
    assert 60 <= coverage <= 61


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


def test_fragments_include_only_modifications_on_their_side_of_cleavage():
    observed = isogen.calc_pep_fragments(
        "EM[Oxidation]E", ion_types="by"
    )
    oxidation = 15.994915

    assert observed["b1"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="E", ion_type="b"),
        abs=MONO_TOLERANCE,
    )
    assert observed["b2"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="EM", ion_type="b")
        + oxidation,
        abs=MONO_TOLERANCE,
    )
    assert observed["y2"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="ME", ion_type="y")
        + oxidation,
        abs=MONO_TOLERANCE,
    )


def test_modified_z_prime_is_one_hydrogen_heavier_than_z():
    observed = isogen.calc_pep_fragments(
        "PEP[Oxidation]TIDE", ion_types=("z", "z'")
    )

    assert observed["z'4"] - observed["z4"] == pytest.approx(
        pyteomics_mass.calculate_mass(formula="H"), abs=1e-9
    )


def test_fragmentation_type_applies_to_proforma_sequences():
    observed = isogen.calc_pep_fragments(
        "PEP[Oxidation]TIDE", fragmentation_type="ECD"
    )
    expected = isogen.calc_pep_fragments(
        "PEP[Oxidation]TIDE", ion_types=("c", "z'")
    )

    assert observed == expected


def test_fragment_terminal_modifications_follow_the_retained_terminus():
    observed = isogen.calc_pep_fragments(
        "[Acetyl]-PEPTIDE-[Amidated]", ion_types="by"
    )

    assert observed["b1"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="P", ion_type="b")
        + 42.010565,
        abs=MONO_TOLERANCE,
    )
    assert observed["y1"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="E", ion_type="y")
        - 0.984016,
        abs=MONO_TOLERANCE,
    )


def test_global_fixed_modifications_are_applied_at_each_retained_site():
    observed = isogen.calc_pep_fragments(
        "<[Carbamidomethyl]@C>ACDC", ion_types="by"
    )
    carbamidomethyl = 57.021464

    assert observed["b2"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="AC", ion_type="b")
        + carbamidomethyl,
        abs=MONO_TOLERANCE,
    )
    assert observed["y3"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="CDC", ion_type="y")
        + 2 * carbamidomethyl,
        abs=MONO_TOLERANCE,
    )


def test_ambiguous_localization_fragments_default_to_reject():
    observed = isogen.calc_pep_fragments(
        "AS[#g1]T[Phospho#g1]K", ion_types="b"
    )

    assert set(observed) == {"b1", "b3"}
    assert observed["b3"] == pytest.approx(
        pyteomics_mass.calculate_mass(sequence="AST", ion_type="b")
        + 79.966331,
        abs=MONO_TOLERANCE,
    )


def test_ambiguous_rule_both_returns_all_distinct_fragment_masses():
    observed = isogen.calc_pep_fragments(
        "AS[#g1]T[Phospho#g1]K",
        ion_types="b",
        ambiguous_rule="both",
    )
    unmodified = pyteomics_mass.calculate_mass(sequence="AS", ion_type="b")

    assert "b2" not in observed
    assert observed["b2#1"] == pytest.approx(unmodified, abs=MONO_TOLERANCE)
    assert observed["b2#2"] == pytest.approx(
        unmodified + 79.966331, abs=MONO_TOLERANCE
    )


def test_multiple_unlocalized_modifications_respect_available_sites():
    sequence = "[Phospho]^2?PEP"
    assert isogen.calc_pep_fragments(sequence, ion_types="b") == {}

    observed = isogen.calc_pep_fragments(
        sequence, ion_types="b", ambiguous_rule="both"
    )
    b1 = pyteomics_mass.calculate_mass(sequence="P", ion_type="b")
    b2 = pyteomics_mass.calculate_mass(sequence="PE", ion_type="b")
    phospho = 79.966331
    assert observed["b1#1"] == pytest.approx(b1, abs=MONO_TOLERANCE)
    assert observed["b1#2"] == pytest.approx(
        b1 + phospho, abs=MONO_TOLERANCE
    )
    assert observed["b2#1"] == pytest.approx(
        b2 + phospho, abs=MONO_TOLERANCE
    )
    assert observed["b2#2"] == pytest.approx(
        b2 + 2 * phospho, abs=MONO_TOLERANCE
    )


def test_ambiguous_residue_fragments_follow_ambiguity_rule():
    assert "b2" not in isogen.calc_pep_fragments("ABK", ion_types="b")

    observed = isogen.calc_pep_fragments(
        "ABK", ion_types="b", ambiguous_rule="both"
    )
    expected = sorted((
        pyteomics_mass.calculate_mass(sequence="AD", ion_type="b"),
        pyteomics_mass.calculate_mass(sequence="AN", ion_type="b"),
    ))
    assert observed["b2#1"] == pytest.approx(expected[0], abs=MONO_TOLERANCE)
    assert observed["b2#2"] == pytest.approx(expected[1], abs=MONO_TOLERANCE)


def test_fragment_ambiguity_rule_is_validated():
    with pytest.raises(ValueError, match="ambiguous_rule"):
        isogen.calc_pep_fragments("PEPTIDE", ambiguous_rule="midpoint")
