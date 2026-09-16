"""Compare sequence- and mass-based IsoGen predictions for FASTA proteins."""

import argparse
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np

if __package__:
    from .. import isodist
    from ..mass import read_fasta
else:
    from isogen import isodist
    from isogen.mass import read_fasta

import matplotlib

matplotlib.use("TkAgg")

TOPDIR = Path(r"C:\Data\IsoNN\fastas")
FASTA_SUFFIXES = {".fa", ".faa", ".fas", ".fasta"}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
BIN_WIDTH = 2_000


def cosine_similarity(left, right):
    """Return the cosine similarity between two intensity vectors."""
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0:
        return np.nan
    return float(np.dot(left, right) / denominator)


def maximum_percent_deviation(left, right):
    """Return the largest pointwise error as a percentage of peak intensity."""
    maximum = max(np.max(np.abs(left)), np.max(np.abs(right)))
    if maximum == 0:
        return np.nan
    return float(100 * np.max(np.abs(left - right)) / maximum)


def absolute_apex_mass_deviation(sequence_prediction, mass_prediction):
    """Return the absolute difference between distribution apex masses."""
    sequence_apex = sequence_prediction[np.argmax(sequence_prediction[:, 1]), 0]
    mass_apex = mass_prediction[np.argmax(mass_prediction[:, 1]), 0]
    return float(abs(sequence_apex - mass_apex))


def plot_butterfly(sequence_prediction, mass_prediction, sequence, css):
    """Plot mirrored sequence/mass distributions beside the protein sequence."""
    monoisotopic_mass = sequence_prediction[0, 0]
    isotope_number = np.arange(len(sequence_prediction))
    figure, axes = plt.subplots(
        1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": (2, 1)}
    )

    axes[0].vlines(
        isotope_number,
        0,
        sequence_prediction[:, 1],
        color="C0",
        label="Sequence",
    )
    axes[0].vlines(
        isotope_number,
        0,
        -mass_prediction[:, 1],
        color="C1",
        label="Mass",
    )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_title(f"Monoisotopic mass: {monoisotopic_mass:,.2f} Da")
    axes[0].set_xlabel("Isotope number")
    axes[0].set_ylabel("Relative intensity")
    axes[0].legend()

    axes[1].axis("off")
    axes[1].set_title(f"Sequence (CSS: {css:.4f})")
    axes[1].text(
        0,
        1,
        textwrap.fill(sequence, width=50),
        family="monospace",
        va="top",
        transform=axes[1].transAxes,
    )
    figure.tight_layout()
    return figure, axes


def fasta_files(directory):
    """Return all FASTA files below a directory in stable order."""
    return sorted(
        path
        for path in directory.rglob("*")
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def analyze_fastas(directory, method="BRAIN", isolen=128, butterfly_cutoff=0.7):
    """Calculate CSS, MPD, and apex deviation by mass bin and M+C percentage."""
    files = fasta_files(directory)
    if not files:
        raise FileNotFoundError(f"No FASTA files found below {directory}")

    css_by_bin = {}
    mpd_by_bin = {}
    apex_deviation_by_bin = {}
    mc_css = []
    mc_mpd = []
    mc_apex_deviation = []
    processed = 0
    skipped = 0

    for path in files:
        for sequence in read_fasta(path).values():
            sequence = sequence.upper().rstrip("*")
            if not sequence or not set(sequence) <= AMINO_ACIDS:
                skipped += 1
                continue

            sequence_prediction = isodist(
                sequence,
                type="PEPTIDE",
                isolen=isolen,
                method=method,
            )
            monoisotopic_mass = float(sequence_prediction[0, 0])
            mass_prediction = isodist(
                monoisotopic_mass,
                type="PEPTIDE",
                isolen=isolen,
                method=method,
            )
            sequence_intensities = sequence_prediction[:, 1]
            mass_intensities = mass_prediction[:, 1]
            css = cosine_similarity(sequence_intensities, mass_intensities)
            mpd = maximum_percent_deviation(sequence_intensities, mass_intensities)
            apex_deviation = absolute_apex_mass_deviation(
                sequence_prediction, mass_prediction
            )
            if np.isfinite(css) and np.isfinite(mpd) and np.isfinite(apex_deviation):
                if butterfly_cutoff is not None and css < butterfly_cutoff:
                    plot_butterfly(
                        sequence_prediction, mass_prediction, sequence, css
                    )
                bin_start = int(monoisotopic_mass // BIN_WIDTH) * BIN_WIDTH
                css_by_bin.setdefault(bin_start, []).append(css)
                mpd_by_bin.setdefault(bin_start, []).append(mpd)
                apex_deviation_by_bin.setdefault(bin_start, []).append(
                    apex_deviation
                )
                mc_percentage = 100 * (sequence.count("M") + sequence.count("C")) / len(sequence)
                mc_css.append((mc_percentage, css))
                mc_mpd.append((mc_percentage, mpd))
                mc_apex_deviation.append((mc_percentage, apex_deviation))
                processed += 1
            else:
                skipped += 1

    if not css_by_bin:
        raise ValueError("No valid protein sequences were found in the FASTA files")

    print(
        f"Processed {processed} sequences from {len(files)} FASTA files; "
        f"skipped {skipped}."
    )
    return (
        css_by_bin,
        mpd_by_bin,
        apex_deviation_by_bin,
        mc_css,
        mc_mpd,
        mc_apex_deviation,
    )


def plot_deviations(
    css_by_bin,
    mpd_by_bin,
    apex_deviation_by_bin,
    mc_css,
    mc_mpd,
    mc_apex_deviation,
):
    """Plot all comparison metrics by mass bin and M+C percentage."""
    bin_starts = np.array(sorted(css_by_bin))
    mass_kda = (bin_starts + BIN_WIDTH / 2) / 1_000
    mean_css = [np.mean(css_by_bin[start]) for start in bin_starts]
    min_css = [np.min(css_by_bin[start]) for start in bin_starts]
    mean_mpd = [np.mean(mpd_by_bin[start]) for start in bin_starts]
    max_mpd = [np.max(mpd_by_bin[start]) for start in bin_starts]
    mean_apex_deviation = [
        np.mean(apex_deviation_by_bin[start]) for start in bin_starts
    ]
    min_apex_deviation = [
        np.min(apex_deviation_by_bin[start]) for start in bin_starts
    ]
    max_apex_deviation = [
        np.max(apex_deviation_by_bin[start]) for start in bin_starts
    ]

    figure, axes = plt.subplots(3, 2, figsize=(15, 18))

    axes[0, 0].plot(mass_kda, mean_css, marker="o", label="Mean CSS")
    axes[0, 0].plot(mass_kda, min_css, marker="o", label="Minimum CSS")
    axes[0, 0].set_xlabel("Monoisotopic mass (kDa)")
    axes[0, 0].set_ylabel("Cosine similarity score")
    axes[0, 0].set_ylim(0, 1.01)
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    mc_percentage, protein_css = zip(*mc_css)
    axes[0, 1].scatter(mc_percentage, protein_css, s=8, alpha=0.25)
    axes[0, 1].set_xlabel("M + C residues (%)")
    axes[0, 1].set_ylabel("Cosine similarity score")
    axes[0, 1].set_ylim(0, 1.01)
    axes[0, 1].grid(alpha=0.25)

    axes[1, 0].plot(mass_kda, mean_mpd, marker="o", label="Mean MPD")
    axes[1, 0].plot(mass_kda, max_mpd, marker="o", label="Maximum MPD")
    axes[1, 0].set_xlabel("Monoisotopic mass (kDa)")
    axes[1, 0].set_ylabel("Maximum percent deviation (%)")
    axes[1, 0].set_ylim(bottom=0)
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    mc_percentage, protein_mpd = zip(*mc_mpd)
    axes[1, 1].scatter(mc_percentage, protein_mpd, s=8, alpha=0.25)
    axes[1, 1].set_xlabel("M + C residues (%)")
    axes[1, 1].set_ylabel("Maximum percent deviation (%)")
    axes[1, 1].set_ylim(bottom=0)
    axes[1, 1].grid(alpha=0.25)

    axes[2, 0].plot(
        mass_kda, mean_apex_deviation, marker="o", label="Mean deviation"
    )
    axes[2, 0].plot(
        mass_kda, min_apex_deviation, marker="o", label="Minimum deviation"
    )
    axes[2, 0].plot(
        mass_kda, max_apex_deviation, marker="o", label="Maximum deviation"
    )
    axes[2, 0].set_xlabel("Monoisotopic mass (kDa)")
    axes[2, 0].set_ylabel("Absolute apex mass deviation (Da)")
    axes[2, 0].set_ylim(bottom=0)
    axes[2, 0].grid(alpha=0.25)
    axes[2, 0].legend()

    mc_percentage, protein_apex_deviation = zip(*mc_apex_deviation)
    axes[2, 1].scatter(
        mc_percentage, protein_apex_deviation, s=8, alpha=0.25
    )
    axes[2, 1].set_xlabel("M + C residues (%)")
    axes[2, 1].set_ylabel("Absolute apex mass deviation (Da)")
    axes[2, 1].set_ylim(bottom=0)
    axes[2, 1].grid(alpha=0.25)

    figure.tight_layout()
    return figure, axes


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        default=TOPDIR,
        help=f"FASTA directory (default: {TOPDIR})",
    )
    parser.add_argument(
        "--method",
        choices=("FFT", "NN", "BRAIN"),
        default="FFT",
        help="IsoGen prediction method (default: %(default)s)",
    )
    parser.add_argument(
        "--isolen",
        type=int,
        default=128,
        help="Number of isotope positions (default: %(default)s)",
    )
    parser.add_argument(
        "--butterfly-cutoff",
        type=float,
        default=0.0,
        help="Create a butterfly plot below this CSS (default: %(default)s)",
    )
    parser.add_argument("--output", type=Path, help="Save the plot instead of displaying it")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if arguments.isolen < 1:
        raise ValueError("--isolen must be positive")
    if not 0 <= arguments.butterfly_cutoff <= 1:
        raise ValueError("--butterfly-cutoff must be between 0 and 1")

    results = analyze_fastas(
        arguments.directory,
        arguments.method,
        arguments.isolen,
        arguments.butterfly_cutoff,
    )
    figure, _ = plot_deviations(*results)
    if arguments.output:
        figure.savefig(arguments.output, dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    main()
