"""Independent checks for heavily modified, labeled peptides and ion termini."""

import numpy as np
import pytest
from pyteomics import mass

import isogen
from benchmarks.proforma_reference import independently_modified, pyteomics_envelope


CASES = [
    ("[Acetyl]-AC[Carbamidomethyl]DM[Oxidation]ST[Phospho]NK[UNIMOD:259]",
     "ACDMSTNK", ("C2H2O", "C2H3NO", "O", "HPO3"), (("C", 13, 6), ("N", 15, 2)), {}),
    ("<[Carbamidomethyl]@C><[Oxidation]@M>ACMCMK[UNIMOD:259]",
     "ACMCMK", ("C2H3NO", "C2H3NO", "O", "O"), (("C", 13, 6), ("N", 15, 2)), {}),
    ("<13C><15N>[Acetyl]-AC[Carbamidomethyl]DM[Oxidation]NK",
     "ACDMNK", ("C2H2O", "C2H3NO", "O"), (), {"C": 13, "N": 15}),
    ("N[Glycan:HexNAc2Hex3NeuAc1]M[Oxidation]C[Carbamidomethyl]K[UNIMOD:259]",
     "NMCK", ("C45H73N3O33", "O", "C2H3NO"), (("C", 13, 6), ("N", 15, 2)), {}),
    ("[Acetyl]-M[Oxidation]S[Phospho]N[Deamidated]K[UNIMOD:259]-[Amidated]",
     "MSNK", ("C2H2O", "O", "HPO3", "H-1N-1O", "HNO-1"), (("C", 13, 6), ("N", 15, 2)), {}),
    ("A[Formula:H-4[2H4]]C[Formula:S-1Se]M[Oxidation]K[UNIMOD:259]",
     "ACMK", ("S-1Se", "O"), (("H", 2, 4), ("C", 13, 6), ("N", 15, 2)), {}),
]


@pytest.mark.parametrize("proforma,plain,deltas,replacements,labels", CASES)
def test_heavily_modified_peptides_match_pyteomics(proforma, plain, deltas, replacements, labels):
    natural, fixed = independently_modified(plain, deltas, replacements, labels)
    origin, expected = pyteomics_envelope(natural, fixed, isolen=128)
    observed = isogen.isodist(proforma, isolen=128)
    # Native natural-abundance data and Pyteomics NIST data differ slightly.
    np.testing.assert_allclose(observed[:, 1], expected, atol=0.0015, rtol=0.006)
    assert observed[0, 0] == pytest.approx(origin, abs=5e-5)


@pytest.mark.parametrize("ion_type", ["a", "a+1", "b", "c", "x", "x+1", "y", "y-1", "z", "z'"])
@pytest.mark.parametrize("sequence", ["ENQS", "ACM", "PEPTIDE"])
def test_terminal_envelopes_match_pyteomics(sequence, ion_type):
    reference_type = "z+1" if ion_type == "z'" else ion_type
    natural = dict(mass.Composition(sequence=sequence, ion_type="y" if ion_type == "y-1" else reference_type))
    if ion_type == "y-1":  # Pyteomics has no entry for this radical y ion.
        natural["H"] -= 1
    origin, expected = pyteomics_envelope(natural, isolen=64)
    observed = isogen.isodist(sequence, ion_type=ion_type, isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected, atol=0.0002, rtol=0.005)
    assert observed[0, 0] == pytest.approx(origin, abs=5e-5)


def test_label_and_modification_counts_include_correct_fragment_termini():
    natural, fixed = independently_modified("AMK", ("O",), (("C", 13, 6), ("N", 15, 2)), {"H": 2})
    fixed[("H", 2)] -= 2  # b-ion has no intact terminal H2O
    natural["O"] -= 1
    origin, expected = pyteomics_envelope(natural, fixed, isolen=64)
    observed = isogen.isodist("<D>AM[Oxidation]K[UNIMOD:259]", ion_type="b", isolen=64)
    np.testing.assert_allclose(observed[:, 1], expected, atol=0.0002, rtol=0.005)
    assert observed[0, 0] == pytest.approx(origin, abs=5e-5)


def test_merged_batch_matches_pyteomics_for_every_fragment():
    sequence = "ENQSMK"
    batch = isogen.calc_pep_fragment_isodists(sequence, ion_types=("a", "b", "c", "x", "y", "z'"))
    for label, origin, observed in zip(batch.labels, batch.masses, batch.intensities):
        ion = "z+1" if label.startswith("z'") else label[0]
        length = int(label[2:] if label.startswith("z'") else label[1:])
        fragment = sequence[:length] if ion[0] in "abc" else sequence[-length:]
        natural = dict(mass.Composition(sequence=fragment, ion_type=ion))
        expected_origin, expected = pyteomics_envelope(natural)
        assert origin == pytest.approx(expected_origin, abs=5e-5)
        np.testing.assert_allclose(observed, expected, atol=0.0002, rtol=0.005)
