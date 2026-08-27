"""Benchmark RNA isotope-distribution methods with generated sequences.

Run with ``python -m isogen.timing_test_rna``. Unlike ``timing_test.py``,
this script does not require external sequence files.
"""

import argparse
import random
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import rcParams

if __package__:
    from .isogen import isodist
    from .mass import calc_rna_mass
else:
    from isogen import isodist
    from mass import calc_rna_mass


RNA_ALPHABET = "ACGU"
MAX_NN_RNA_LENGTH = 475
METHODS = ("FFT", "NN", "BRAIN")


def remove_outliers(data, multiplier=2):
    """Remove IQR outliers from a sequence of timing values."""
    if len(data) < 4:
        return data

    sorted_data = sorted(data)
    q1 = sorted_data[len(sorted_data) // 4]
    q3 = sorted_data[(len(sorted_data) * 3) // 4]
    iqr = q3 - q1
    lower_bound = q1 - multiplier * iqr
    upper_bound = q3 + multiplier * iqr
    return [value for value in data if lower_bound <= value <= upper_bound]


def calculate_cosine_similarity(left, right):
    """Return the cosine similarity of two isotope-intensity vectors."""
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0:
        return 0.0
    return float(np.dot(left, right) / denominator)


def average_metric_map(metric_map, remove_outlier_values=False):
    """Average each sorted measurement bucket."""
    averaged = {}
    for key, values in sorted(metric_map.items()):
        if remove_outlier_values:
            values = remove_outliers(values)
        averaged[key] = np.average(values)
    return averaged


def append_value(values_by_key, key, value):
    values_by_key.setdefault(key, []).append(value)


def random_rna_sequences(max_length, samples_per_length, random_generator):
    """Yield random RNA sequences for every requested nucleotide length."""
    for length in range(1, max_length + 1):
        for _ in range(samples_per_length):
            yield "".join(random_generator.choices(RNA_ALPHABET, k=length))


def get_rna_distribution(target, method):
    return isodist(target, type="RNA", isolen=128, method=method, dist_only=True)


def benchmark_rna_sequences(sequences, mass_bin_width):
    """Time all RNA engines and collect agreement with sequence FFT output."""
    timings_by_length = {}
    timings_by_mass = {}
    similarities_by_length = {}
    similarities_by_mass = {}

    for sequence in sequences:
        length = len(sequence)
        mass = calc_rna_mass(sequence)
        mass_key = int(mass // mass_bin_width) * mass_bin_width

        start = time.perf_counter()
        fft_sequence = get_rna_distribution(sequence, "FFT")
        fft_sequence_time = time.perf_counter() - start
        append_value(timings_by_length, ("FFT-Seq", length), fft_sequence_time)
        append_value(timings_by_mass, ("FFT-Seq", mass_key), fft_sequence_time)

        for method in METHODS:
            target_name = f"{method}-Mass"
            start = time.perf_counter()
            distribution = get_rna_distribution(mass, method)
            elapsed = time.perf_counter() - start
            similarity = calculate_cosine_similarity(distribution, fft_sequence)

            append_value(timings_by_length, (target_name, length), elapsed)
            append_value(timings_by_mass, (target_name, mass_key), elapsed)
            append_value(similarities_by_length, (target_name, length), similarity)
            append_value(similarities_by_mass, (target_name, mass_key), similarity)

        for method in ("NN", "BRAIN"):
            target_name = f"{method}-Seq"
            start = time.perf_counter()
            distribution = get_rna_distribution(sequence, method)
            elapsed = time.perf_counter() - start
            similarity = calculate_cosine_similarity(distribution, fft_sequence)

            append_value(timings_by_length, (target_name, length), elapsed)
            append_value(timings_by_mass, (target_name, mass_key), elapsed)
            append_value(similarities_by_length, (target_name, length), similarity)
            append_value(similarities_by_mass, (target_name, mass_key), similarity)

    return (
        average_metric_map(timings_by_length, remove_outlier_values=True),
        average_metric_map(timings_by_mass, remove_outlier_values=True),
        average_metric_map(similarities_by_length),
        average_metric_map(similarities_by_mass),
    )


def plot_metric(axis, metrics, title, xlabel, ylabel):
    """Plot a metric map whose keys are ``(method, measurement)`` pairs."""
    methods = sorted({method for method, _ in metrics})
    for method in methods:
        points = [(key, value) for (name, key), value in metrics.items() if name == method]
        if points:
            x_values, y_values = zip(*points)
            axis.scatter(x_values, y_values, label=method)

    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.legend()


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-length",
        type=int,
        default=MAX_NN_RNA_LENGTH,
        help=f"Maximum RNA length (1-{MAX_NN_RNA_LENGTH}; default: %(default)s).",
    )
    parser.add_argument(
        "--samples-per-length",
        type=int,
        default=10,
        help="Random sequences generated at each length (default: %(default)s).",
    )
    parser.add_argument(
        "--mass-bin-width",
        type=int,
        default=1000,
        help="Mass-bin width in Da (default: %(default)s).",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed (default: %(default)s).")
    parser.add_argument("--output", help="Optional image path instead of displaying the plot.")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if not 1 <= arguments.max_length <= MAX_NN_RNA_LENGTH:
        raise ValueError(
            f"--max-length must be between 1 and {MAX_NN_RNA_LENGTH} to support all models."
        )
    if arguments.samples_per_length < 1:
        raise ValueError("--samples-per-length must be positive.")
    if arguments.mass_bin_width < 1:
        raise ValueError("--mass-bin-width must be positive.")

    random_generator = random.Random(arguments.seed)
    sequences = random_rna_sequences(
        arguments.max_length, arguments.samples_per_length, random_generator
    )
    timing_by_length, timing_by_mass, similarity_by_length, similarity_by_mass = (
        benchmark_rna_sequences(sequences, arguments.mass_bin_width)
    )

    figure, axes = plt.subplots(2, 2, figsize=(16, 12))
    plot_metric(axes[0, 0], timing_by_length, "RNA speed by length", "Nucleotide count", "Average time (sec)")
    plot_metric(axes[0, 1], similarity_by_length, "RNA agreement by length", "Nucleotide count", "Cosine similarity")
    plot_metric(axes[1, 0], timing_by_mass, "RNA speed by mass", "Mass (Da)", "Average time (sec)")
    plot_metric(axes[1, 1], similarity_by_mass, "RNA agreement by mass", "Mass (Da)", "Cosine similarity")
    figure.tight_layout()

    if arguments.output:
        figure.savefig(arguments.output, dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    rcParams["ps.useafm"] = True
    rcParams["ps.fonttype"] = 42
    rcParams["pdf.fonttype"] = 42
    rcParams["font.family"] = "Arial"
    main()
