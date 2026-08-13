"""Regression tests for the native BRAIN isotope implementation."""

import ctypes

import numpy as np
import pytest

from isogen.isogenwrapper import isogen_c_lib


FLOAT_PTR = ctypes.POINTER(ctypes.c_float)
INT_PTR = ctypes.POINTER(ctypes.c_int)


def _configure_native_functions():
    isogen_c_lib.brain_list_to_dist.argtypes = [INT_PTR, ctypes.c_int, FLOAT_PTR]
    isogen_c_lib.brain_list_to_dist.restype = ctypes.c_float

    for name in ("brain_rna_mass_to_dist", "brain_pep_mass_to_dist"):
        function = getattr(isogen_c_lib, name)
        function.argtypes = [ctypes.c_float, FLOAT_PTR, ctypes.c_int, ctypes.c_int]
        function.restype = ctypes.c_float

    for name in ("brain_rna_seq_to_dist", "brain_pep_seq_to_dist"):
        function = getattr(isogen_c_lib, name)
        function.argtypes = [ctypes.c_char_p, FLOAT_PTR, ctypes.c_int, ctypes.c_int]
        function.restype = ctypes.c_float

    isogen_c_lib.fft_list_to_dist.argtypes = [INT_PTR, ctypes.c_int, FLOAT_PTR]
    isogen_c_lib.fft_list_to_dist.restype = ctypes.c_float

    isogen_c_lib.normalize_isodist.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]
    isogen_c_lib.normalize_isodist.restype = ctypes.c_double


_configure_native_functions()


def _float_pointer(values):
    return values.ctypes.data_as(FLOAT_PTR)


def test_brain_recurrence_remains_finite_for_large_formula():
    formula = (ctypes.c_int * 5)(10_000, 16_000, 2_000, 3_000, 200)
    observed = np.zeros(1024, dtype=np.float32)

    result = isogen_c_lib.brain_list_to_dist(
        formula, len(observed), _float_pointer(observed)
    )

    assert np.isfinite(result)
    assert np.isfinite(observed).all()
    assert observed.sum() == pytest.approx(1.0, abs=1e-6)


def test_brain_formula_entry_point_rejects_invalid_arguments():
    formula = (ctypes.c_int * 5)(10, 20, 2, 3, 0)
    observed = np.zeros(16, dtype=np.float32)

    assert (
        isogen_c_lib.brain_list_to_dist(
            None, len(observed), _float_pointer(observed)
        )
        == -1.0
    )
    assert isogen_c_lib.brain_list_to_dist(formula, 0, _float_pointer(observed)) == -1.0
    assert isogen_c_lib.brain_list_to_dist(formula, len(observed), None) == -1.0


def test_normalization_excludes_clamped_negative_values_from_sum():
    values = (ctypes.c_double * 3)(0.8, 0.3, -0.1)

    maximum = isogen_c_lib.normalize_isodist(values, len(values))

    assert list(values) == pytest.approx([0.8 / 1.1, 0.3 / 1.1, 0.0])
    assert sum(values) == pytest.approx(1.0)
    assert maximum == pytest.approx(0.8 / 1.1)


@pytest.mark.parametrize(
    ("function_name", "input_value"),
    [
        ("brain_rna_mass_to_dist", 320_000.0),
        ("brain_pep_mass_to_dist", 250_000.0),
        ("brain_rna_seq_to_dist", b"A" * 1_000),
        ("brain_pep_seq_to_dist", b"A" * 2_000),
    ],
)
def test_brain_entry_points_handle_large_inputs(function_name, input_value):
    observed = np.zeros(128, dtype=np.float32)

    result = getattr(isogen_c_lib, function_name)(
        input_value, _float_pointer(observed), len(observed), 0
    )

    assert result > 0
    assert np.isfinite(observed).all()


@pytest.mark.parametrize(
    ("function_name", "input_value"),
    [
        ("brain_rna_mass_to_dist", 10_000.0),
        ("brain_pep_mass_to_dist", 10_000.0),
        ("brain_rna_seq_to_dist", b"ACGU"),
        ("brain_pep_seq_to_dist", b"PEPTIDE"),
    ],
)
def test_brain_entry_points_reject_invalid_output_geometry(function_name, input_value):
    observed = np.zeros(16, dtype=np.float32)
    function = getattr(isogen_c_lib, function_name)

    assert function(input_value, _float_pointer(observed), len(observed), -1) == -1.0
    assert function(input_value, None, len(observed), 0) == -1.0


def test_rna_sequence_counts_each_nucleotide_phosphate_once():
    observed = np.zeros(32, dtype=np.float32)
    expected = np.zeros_like(observed)
    # Sum of the A, C, G, and U vectors; phosphorus has only one stable isotope.
    formula = (ctypes.c_int * 5)(38, 43, 15, 28, 0)

    isogen_c_lib.fft_rna_seq_to_dist(
        b"ACGU", _float_pointer(observed), len(observed), 0
    )
    isogen_c_lib.fft_list_to_dist(
        formula, len(expected), _float_pointer(expected)
    )
    expected /= expected.max()

    np.testing.assert_allclose(observed, expected, rtol=1e-6, atol=1e-8)
