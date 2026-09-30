import numpy as np
import pytest

import isogen


def test_batch_uses_fragment_composition_and_ion_termini():
    batch = isogen.calc_pep_fragment_isodists("PEPTIDE", ion_types="by")
    assert batch.labels == tuple(isogen.calc_pep_fragments("PEPTIDE", ion_types="by"))
    assert batch.intensities.shape == (len(batch.labels), 128)
    b1 = batch.intensities[batch.labels.index("b1")]
    y1 = batch.intensities[batch.labels.index("y1")]
    np.testing.assert_allclose(
        b1[:6], isogen.isodist("C5H7NO", type="ATOM", dist_only=True)[:6],
        rtol=2e-5, atol=1e-8,
    )
    assert not np.array_equal(b1, y1)


def test_known_modification_formula_changes_fragment_envelope():
    plain = isogen.calc_pep_fragment_isodists("SHHS")
    modified = isogen.calc_pep_fragment_isodists("S[Acetylation]HHS")
    assert modified.masses[modified.labels.index("b1")] > plain.masses[plain.labels.index("b1")]
    assert not np.array_equal(modified.intensities[modified.labels.index("b1")],
                              plain.intensities[plain.labels.index("b1")])
    np.testing.assert_array_equal(modified.intensities[modified.labels.index("y1")],
                                  plain.intensities[plain.labels.index("y1")])


def test_mass_only_modification_requires_composition():
    with pytest.raises(ValueError, match="composition is unknown"):
        isogen.calc_pep_fragment_isodists("S[+42.0106]HHS")


def test_batch_handles_empty_sequence_and_terminal_ion_series():
    empty = isogen.calc_pep_fragment_isodists("")
    assert empty.labels == ()
    assert empty.intensities.shape == (0, 128)
    batch = isogen.calc_pep_fragment_isodists("PEPTIDE", ion_types=("a+1", "c", "x+1", "z'"))
    assert set(label[0] for label in batch.labels) == {"a", "c", "x", "z"}
    assert np.all(batch.intensities.max(axis=1) > 0)
