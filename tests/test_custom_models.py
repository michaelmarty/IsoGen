"""Tests for loading neural-network weights from caller-provided files."""

from pathlib import Path

import numpy as np
import pytest

import isogen
from isogen import isogenwrapper


MODEL_DIRECTORY = Path(isogenwrapper.current_path) / "models"


@pytest.mark.parametrize(
    ("input_value", "analyte_type", "isolen", "model_name"),
    [
        ("A" * 10, "PEPTIDE", 16, "isogenpep_model_16.bin"),
        ("A" * 60, "PEPTIDE", 64, "isogenpep_model_64.bin"),
        ("A" * 301, "PEPTIDE", 128, "isogenpep_model_128.bin"),
        ("A" * 10, "RNA", 64, "isogenrna_model_64.bin"),
        ("A" * 201, "RNA", 128, "isogenrna_model_128.bin"),
        (1_000.0, "PEPTIDE", 8, "isogenmass_model_8.bin"),
        (5_000.0, "PEPTIDE", 32, "isogenmass_model_32.bin"),
        (20_000.0, "PEPTIDE", 64, "isogenmass_model_64.bin"),
        (60_000.0, "PEPTIDE", 128, "isogenmass_model_128.bin"),
        (10_000.0, "RNA", 32, "isogen_rnaveragine_model32.bin"),
        (30_000.0, "RNA", 64, "isogen_rnaveragine_model64.bin"),
        (70_000.0, "RNA", 128, "isogen_rnaveragine_model128.bin"),
    ],
)
def test_provided_model_matches_default(
    input_value, analyte_type, isolen, model_name
):
    """Provided model files should reproduce their compiled-in defaults."""
    expected = isogenwrapper.nn_gen_isodist(
        input_value,
        type=analyte_type,
        isolen=isolen,
    )
    observed = isogenwrapper.nn_gen_isodist(
        input_value,
        type=analyte_type,
        isolen=isolen,
        model_path=MODEL_DIRECTORY / model_name,
    )

    np.testing.assert_allclose(observed, expected)


def test_gen_isodist_forwards_custom_model_path():
    """The general dispatcher should pass custom model paths to NN calls."""
    model_path = MODEL_DIRECTORY / "isogenmass_model_32.bin"
    expected = isogenwrapper.nn_gen_isodist(
        5_000.0,
        type="PEPTIDE",
        isolen=32,
        model_path=model_path,
    )

    observed = isogenwrapper.gen_isodist(
        5_000.0,
        type="PEPTIDE",
        isolen=32,
        method="NN",
        model_path=model_path,
    )

    np.testing.assert_allclose(observed, expected)


def test_invalid_custom_model_path_raises_value_error(tmp_path):
    """Native custom-model loading failures should be visible to callers."""
    with pytest.raises(ValueError, match="Unable to use custom model file"):
        isogenwrapper.nn_gen_isodist(
            5_000.0,
            type="PEPTIDE",
            isolen=32,
            model_path=tmp_path / "missing.bin",
        )


def test_public_isodist_custom_combines_mass_axis_and_custom_intensities():
    """The public custom API should expose native output with a mass axis."""
    model_path = MODEL_DIRECTORY / "isogenrna_model_64.bin"
    expected_intensities = isogenwrapper.nn_gen_isodist(
        "AUGC",
        type="RNA",
        isolen=64,
        model_path=model_path,
    )

    observed = isogen.isodist_custom(
        "AUGC",
        model_file=model_path,
        isolen=64,
        type="rna",
    )

    assert observed.shape == (64, 2)
    np.testing.assert_allclose(observed[:, 1], expected_intensities)
    assert np.all(np.diff(observed[:, 0]) > 0)


def test_public_isodist_custom_rejects_atom_type():
    """Elemental formulas have no compatible custom neural-network model."""
    with pytest.raises(ValueError, match="Custom models support"):
        isogen.isodist_custom(
            "H2O",
            model_file="unused.bin",
            isolen=16,
            type="ATOM",
        )
