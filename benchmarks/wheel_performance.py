"""Benchmark an installed IsoGen wheel and write machine-readable results."""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from pathlib import Path
from typing import Callable

import isogen
import numpy as np


PEPTIDE_SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 20
RNA_SEQUENCE = "AUGC" * 100


def cases() -> dict[str, Callable[[], np.ndarray]]:
    """Return representative native workloads for every calculation method."""
    return {
        "fft_peptide_mass": lambda: isogen.isodist(
            50_000.0, type="PEPTIDE", isolen=128, method="FFT", dist_only=True
        ),
        "fft_peptide_sequence": lambda: isogen.isodist(
            PEPTIDE_SEQUENCE,
            type="PEPTIDE",
            isolen=128,
            method="FFT",
            dist_only=True,
        ),
        "nn_peptide_mass": lambda: isogen.isodist(
            50_000.0, type="PEPTIDE", isolen=128, method="NN", dist_only=True
        ),
        "nn_peptide_sequence": lambda: isogen.isodist(
            PEPTIDE_SEQUENCE,
            type="PEPTIDE",
            isolen=128,
            method="NN",
            dist_only=True,
        ),
        "brain_peptide_mass": lambda: isogen.isodist(
            50_000.0, type="PEPTIDE", isolen=128, method="BRAIN", dist_only=True
        ),
        "brain_peptide_sequence": lambda: isogen.isodist(
            PEPTIDE_SEQUENCE,
            type="PEPTIDE",
            isolen=128,
            method="BRAIN",
            dist_only=True,
        ),
        "fft_rna_sequence": lambda: isogen.isodist(
            RNA_SEQUENCE, type="RNA", isolen=128, method="FFT", dist_only=True
        ),
        "nn_rna_sequence": lambda: isogen.isodist(
            RNA_SEQUENCE, type="RNA", isolen=128, method="NN", dist_only=True
        ),
        "brain_rna_sequence": lambda: isogen.isodist(
            RNA_SEQUENCE, type="RNA", isolen=128, method="BRAIN", dist_only=True
        ),
    }


def time_case(
    function: Callable[[], np.ndarray], target_seconds: float, repeats: int
) -> tuple[int, list[float]]:
    """Measure a call with enough iterations to reduce timer noise."""
    function()
    iterations = 1
    while True:
        start = time.perf_counter()
        for _ in range(iterations):
            function()
        elapsed = time.perf_counter() - start
        if elapsed >= target_seconds or iterations >= 1_048_576:
            break
        iterations *= 2

    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        for _ in range(iterations):
            function()
        samples.append((time.perf_counter() - start) / iterations)
    return iterations, samples


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-seconds", type=float, default=0.15)
    parser.add_argument("--repeats", type=int, default=5)
    arguments = parser.parse_args()

    benchmark_results = {}
    for name, function in cases().items():
        expected = np.asarray(function(), dtype=np.float64)
        if not np.all(np.isfinite(expected)):
            raise RuntimeError(f"{name} returned a non-finite value")
        iterations, samples = time_case(
            function, arguments.target_seconds, arguments.repeats
        )
        benchmark_results[name] = {
            "iterations": iterations,
            "median_seconds": statistics.median(samples),
            "samples_seconds": samples,
            "result": expected.tolist(),
        }

    output = {
        "label": arguments.label,
        "python": sys.version,
        "platform": platform.platform(),
        "native_library": isogen.isogenwrapper.dllpath,
        "benchmarks": benchmark_results,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(output, indent=2), encoding="utf-8")

    for name, result in benchmark_results.items():
        microseconds = result["median_seconds"] * 1_000_000
        print(f"{arguments.label:>12} {name:<28} {microseconds:10.2f} us")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
