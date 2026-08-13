import matplotlib.pyplot as plt
import numpy as np
import time
import random
from matplotlib.pyplot import rcParams
import math
import os
import pandas as pd

if __package__:
    from .isogen import isodist
    from .mass import calc_pep_mass
    from .isogen_tools import (
        pep_makemass,
        peptide_to_aacount,
    )
else:
    from isogen import isodist
    from mass import calc_pep_mass
    from isogen_tools import (
        pep_makemass,
        peptide_to_aacount,
    )

def remove_outliers(data, multiplier=2):
    """
    Remove outliers from a list using the IQR method.

    Args:
        data (list or iterable of float): The input list of numbers.
        multiplier (float): The IQR multiplier (default 1.5). Larger = less strict.

    Returns:
        list: A list with outliers removed.
    """
    if len(data) < 4:
        # Too few data points to determine outliers reliably
        return data

    sorted_data = sorted(data)
    q1 = sorted_data[len(sorted_data) // 4]
    q3 = sorted_data[(len(sorted_data) * 3) // 4]
    iqr = q3 - q1

    lower_bound = q1 - multiplier * iqr
    upper_bound = q3 + multiplier * iqr

    return [x for x in data if lower_bound <= x <= upper_bound]


def calculate_cosine_similarity(left, right):
    """Return the cosine similarity of two isotope-intensity vectors."""
    left = np.asarray(left)
    right = np.asarray(right)
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0:
        return 0.0
    return float(np.dot(left, right) / denominator)

def get_s_count_from_pepmass(mass):
    _, _, isolist = pep_makemass(mass)
    return isolist[4]

def get_mass_s_count_array():
    masses = np.arange(0, 100000, 10)
    transition_masses = []
    previous_s_count = 0
    for mass in masses:
        s_count = get_s_count_from_pepmass(mass)
        if s_count > previous_s_count:
            transition_masses.append(mass)
            previous_s_count = s_count
    return transition_masses


def get_key_from_mass(mass, bin_width):
    key = int(math.floor(mass / bin_width) * bin_width)
    return key

def style_axis(axis, xlabel, ylabel, xlim=None, ylim=None):
    """Apply the shared plot styling used by all benchmark panels."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(2)
    axis.spines["left"].set_linewidth(2)
    axis.tick_params(length=5, width=2, labelsize=tick_size)
    axis.set_xlabel(xlabel, fontsize=label_size)
    axis.set_ylabel(ylabel, fontsize=label_size)
    if xlim is not None:
        axis.set_xlim(*xlim)
    if ylim is not None:
        axis.set_ylim(*ylim)
    axis.legend()

def open_tsv(filename):
    """Open a TSV file and return a list of rows, each as a list of strings."""
    df = pd.read_csv(filename, sep="\t")
    seqs = df["Sequence"].to_list()

    # Drop seqs longer than 1000
    seqs = [seq for seq in seqs if len(seq) <= 1000]
    masses = [calc_pep_mass(seq) for seq in seqs]

    return seqs, masses


def append_value(values_by_key, key, value):
    values_by_key.setdefault(key, []).append(value)


def append_value_for_count_and_mass(
    values_by_count, aa_count, values_by_mass, mass_key, value
):
    append_value(values_by_count, aa_count, value)
    append_value(values_by_mass, mass_key, value)


def average_metric_map(metric_map, remove_outlier_values=False):
    averaged = {}
    for key, values in sorted(metric_map.items()):
        if remove_outlier_values:
            values = remove_outliers(values)
        averaged[key] = np.average(values)
    return averaged


def get_peptide_distribution(target, method):
    return isodist(target, type="PEPTIDE", isolen=128, method=method, dist_only=True)

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

tick_size = 16
label_size = 20

if __name__ == "__main__":
    topdir = r"Z:\Group Share\JGP\PeptideTraining\IntactProtein\Training"
    os.chdir(topdir)
    assessment_files = ["human_protein_seqs.tsv", "mouse_protein_seqs.tsv", "ecoli_protein_seqs.tsv", "yeast_protein_seqs.tsv"]

    seqs = []
    masses = []

    for file in assessment_files:
        if not os.path.exists(file):
            print("File not found: ", file)
            continue
        if file.endswith(".tsv"):
            seqs_tsv, masses_tsv = open_tsv(file)
            seqs.extend(seqs_tsv)
            masses.extend(masses_tsv)
        elif file.endswith(".npz"):
            data = np.load(file)
            seqs.extend(data["seqs"])
            masses.extend(data["masses"])
        else:
            print("Unsupported file format: ", file)
            continue

    fft_seq_times = {}
    fft_mass_times = {}
    nn_seq_times = {}
    nn_mass_times = {}
    brain_mass_times = {}
    brain_seq_times = {}

    fft_mass_css_vals = {}
    nn_seq_css_vals = {}
    nn_mass_css_vals = {}
    brain_mass_css_vals = {}
    brain_seq_css_vals = {}


    #These are the values separated into a mass dictionary (100 Da / bin, and the key is the bottom)
    fft_seq_times_bymass = {}
    fft_mass_times_bymass = {}
    nn_seq_times_bymass = {}
    nn_mass_times_bymass = {}
    brain_mass_times_bymass = {}
    brain_seq_times_bymass = {}

    fft_mass_css_vals_bymass = {}
    nn_seq_css_vals_bymass = {}
    nn_mass_css_vals_bymass = {}
    brain_mass_css_vals_bymass = {}
    brain_seq_css_vals_bymass = {}

    shuffled_indices = np.arange(0, len(masses))
    random.shuffle(shuffled_indices)


    num_processed = 0
    mass_bin_width = 1000
    for i in shuffled_indices:
        num_processed += 1
        percent_processed = round((num_processed / len(masses)) * 100, 2)
        if num_processed % 1000 == 0:
            print("Processed ", percent_processed, "% of sequences")

        aa_count = peptide_to_aacount(seqs[i])

        mass_key = get_key_from_mass(masses[i], mass_bin_width)

        timestart = time.perf_counter()
        fft_seq_dist = get_peptide_distribution(seqs[i], "FFT")
        fft_seq_time = time.perf_counter() - timestart

        append_value_for_count_and_mass(
            fft_seq_times, aa_count, fft_seq_times_bymass, mass_key, fft_seq_time
        )


        timestart = time.perf_counter()
        fft_mass_dist = get_peptide_distribution(masses[i], "FFT")
        fft_mass_time = time.perf_counter() - timestart
        fft_mass_css = calculate_cosine_similarity(fft_mass_dist, fft_seq_dist)

        append_value_for_count_and_mass(
            fft_mass_times,
            aa_count,
            fft_mass_times_bymass,
            mass_key,
            fft_mass_time,
        )
        append_value_for_count_and_mass(
            fft_mass_css_vals,
            aa_count,
            fft_mass_css_vals_bymass,
            mass_key,
            fft_mass_css,
        )


        timestart = time.perf_counter()
        nn_seq_dist = get_peptide_distribution(seqs[i], "NN")
        nn_seq_time = time.perf_counter() - timestart
        nn_seq_css = calculate_cosine_similarity(nn_seq_dist, fft_seq_dist)

        append_value_for_count_and_mass(
            nn_seq_times, aa_count, nn_seq_times_bymass, mass_key, nn_seq_time
        )
        append_value_for_count_and_mass(
            nn_seq_css_vals, aa_count, nn_seq_css_vals_bymass, mass_key, nn_seq_css
        )

        timestart = time.perf_counter()
        nn_mass_dist = get_peptide_distribution(masses[i], "NN")
        nn_mass_time = time.perf_counter() - timestart
        nn_mass_css = calculate_cosine_similarity(nn_mass_dist, fft_seq_dist)

        append_value_for_count_and_mass(
            nn_mass_times, aa_count, nn_mass_times_bymass, mass_key, nn_mass_time
        )
        append_value_for_count_and_mass(
            nn_mass_css_vals, aa_count, nn_mass_css_vals_bymass, mass_key, nn_mass_css
        )

        timestart = time.perf_counter()
        brain_mass_dist = get_peptide_distribution(masses[i], "BRAIN")
        brain_mass_time = time.perf_counter() - timestart
        brain_mass_css = calculate_cosine_similarity(brain_mass_dist, fft_seq_dist)

        append_value_for_count_and_mass(
            brain_mass_times,
            aa_count,
            brain_mass_times_bymass,
            mass_key,
            brain_mass_time,
        )
        append_value_for_count_and_mass(
            brain_mass_css_vals,
            aa_count,
            brain_mass_css_vals_bymass,
            mass_key,
            brain_mass_css,
        )


        timestart = time.perf_counter()
        brain_seq_dist = get_peptide_distribution(seqs[i], "BRAIN")
        brain_seq_time = time.perf_counter() - timestart
        brain_seq_css = calculate_cosine_similarity(brain_seq_dist, fft_seq_dist)

        append_value_for_count_and_mass(
            brain_seq_times,
            aa_count,
            brain_seq_times_bymass,
            mass_key,
            brain_seq_time,
        )
        append_value_for_count_and_mass(
            brain_seq_css_vals,
            aa_count,
            brain_seq_css_vals_bymass,
            mass_key,
            brain_seq_css,
        )

    fft_seq_times = average_metric_map(fft_seq_times, remove_outlier_values=True)
    fft_seq_times_bymass = average_metric_map(
        fft_seq_times_bymass, remove_outlier_values=True
    )
    fft_mass_times = average_metric_map(fft_mass_times, remove_outlier_values=True)
    fft_mass_times_bymass = average_metric_map(
        fft_mass_times_bymass, remove_outlier_values=True
    )

    nn_seq_times = average_metric_map(nn_seq_times, remove_outlier_values=True)
    nn_seq_times_bymass = average_metric_map(
        nn_seq_times_bymass, remove_outlier_values=True
    )
    nn_mass_times = average_metric_map(nn_mass_times, remove_outlier_values=True)
    nn_mass_times_bymass = average_metric_map(
        nn_mass_times_bymass, remove_outlier_values=True
    )
    fft_mass_css_vals = average_metric_map(fft_mass_css_vals)
    fft_mass_css_vals_bymass = average_metric_map(fft_mass_css_vals_bymass)

    nn_seq_css_vals = average_metric_map(nn_seq_css_vals)
    nn_seq_css_vals_bymass = average_metric_map(nn_seq_css_vals_bymass)
    nn_mass_css_vals = average_metric_map(nn_mass_css_vals)
    nn_mass_css_vals_bymass = average_metric_map(nn_mass_css_vals_bymass)

    brain_mass_times = average_metric_map(
        brain_mass_times, remove_outlier_values=True
    )
    brain_mass_times_bymass = average_metric_map(
        brain_mass_times_bymass, remove_outlier_values=True
    )
    brain_seq_times = average_metric_map(
        brain_seq_times, remove_outlier_values=True
    )
    brain_seq_times_bymass = average_metric_map(
        brain_seq_times_bymass, remove_outlier_values=True
    )
    brain_mass_css_vals = average_metric_map(brain_mass_css_vals)
    brain_mass_css_vals_bymass = average_metric_map(brain_mass_css_vals_bymass)
    brain_seq_css_vals = average_metric_map(brain_seq_css_vals)
    brain_seq_css_vals_bymass = average_metric_map(brain_seq_css_vals_bymass)

    aa_counts = list(fft_seq_times.keys())
    fft_seq_avgtimes = list(fft_seq_times.values())
    fft_mass_avgtimes = list(fft_mass_times.values())
    nn_seq_avgtimes = list(nn_seq_times.values())
    nn_mass_avgtimes = list(nn_mass_times.values())
    brain_mass_avgtimes = list(brain_mass_times.values())
    brain_seq_avgtimes = list(brain_seq_times.values())

    fft_mass_avgcss = list(fft_mass_css_vals.values())
    nn_mass_avgcss = list(nn_mass_css_vals.values())
    nn_seq_avgcss = list(nn_seq_css_vals.values())
    brain_mass_avgcss = list(brain_mass_css_vals.values())
    brain_seq_avgcss = list(brain_seq_css_vals.values())

    masses = list(fft_seq_times_bymass.keys())
    fft_seq_avgtimes_bymass = list(fft_seq_times_bymass.values())
    fft_mass_avgtimes_bymass = list(fft_mass_times_bymass.values())
    nn_seq_avgtimes_bymass = list(nn_seq_times_bymass.values())
    nn_mass_avgtimes_bymass = list(nn_mass_times_bymass.values())
    brain_mass_avgtimes_bymass = list(brain_mass_times_bymass.values())
    brain_seq_avgtimes_bymass = list(brain_seq_times_bymass.values())

    fft_mass_avgcss_bymass = list(fft_mass_css_vals_bymass.values())
    nn_mass_avgcss_bymass = list(nn_mass_css_vals_bymass.values())
    nn_seq_avgcss_bymass = list(nn_seq_css_vals_bymass.values())
    brain_mass_avgcss_bymass = list(brain_mass_css_vals_bymass.values())
    brain_seq_avgcss_bymass = list(brain_seq_css_vals_bymass.values())
    figure, axes = plt.subplots(2, 2, figsize=(16, 12))

    aa_speed_axis = axes[0, 0]
    aa_speed_axis.scatter(aa_counts, fft_seq_avgtimes, label="FFT-Seq")
    aa_speed_axis.scatter(aa_counts, fft_mass_avgtimes, label="FFT-Mass")
    aa_speed_axis.scatter(aa_counts, nn_seq_avgtimes, label="NN-Seq")
    aa_speed_axis.scatter(aa_counts, nn_mass_avgtimes, label="NN-Mass")
    aa_speed_axis.scatter(aa_counts, brain_mass_avgtimes, label="BRAIN-Mass")
    aa_speed_axis.scatter(aa_counts, brain_seq_avgtimes, label="BRAIN-Seq")
    style_axis(
        aa_speed_axis,
        "AA Count",
        "Average Time (sec)",
        xlim=(0, 1000),
        ylim=(0, None),
    )

    aa_css_axis = axes[0, 1]
    aa_css_axis.scatter(aa_counts, fft_mass_avgcss, label="FFT-Mass")
    aa_css_axis.scatter(aa_counts, nn_seq_avgcss, label="NN-Seq")
    aa_css_axis.scatter(aa_counts, nn_mass_avgcss, label="NN-Mass")
    aa_css_axis.scatter(aa_counts, brain_mass_avgcss, label="BRAIN-Mass")
    aa_css_axis.scatter(aa_counts, brain_seq_avgcss, label="BRAIN-Seq")
    style_axis(
        aa_css_axis,
        "AA Count",
        "CSS",
        xlim=(0, 1000),
        ylim=(None, 1),
    )

    mass_speed_axis = axes[1, 0]
    mass_speed_axis.scatter(masses, fft_seq_avgtimes_bymass, label="FFT-Seq")
    mass_speed_axis.scatter(masses, fft_mass_avgtimes_bymass, label="FFT-Mass")
    mass_speed_axis.scatter(masses, nn_seq_avgtimes_bymass, label="NN-Seq")
    mass_speed_axis.scatter(masses, nn_mass_avgtimes_bymass, label="NN-Mass")
    mass_speed_axis.scatter(
        masses, brain_mass_avgtimes_bymass, label="BRAIN-Mass"
    )
    mass_speed_axis.scatter(
        masses, brain_seq_avgtimes_bymass, label="BRAIN-Seq"
    )
    style_axis(
        mass_speed_axis,
        "Mass (Da)",
        "Average Time (sec)",
        ylim=(0, None),
    )

    mass_s_transitions = get_mass_s_count_array()

    mass_css_axis = axes[1, 1]
    mass_css_axis.plot(masses, fft_mass_avgcss_bymass, label="FFT-Mass")
    mass_css_axis.plot(masses, nn_seq_avgcss_bymass, label="NN-Seq")
    mass_css_axis.plot(masses, nn_mass_avgcss_bymass, label="NN-Mass")
    mass_css_axis.plot(
        masses, brain_mass_avgcss_bymass, label="BRAIN-Mass"
    )
    mass_css_axis.plot(
        masses, brain_seq_avgcss_bymass, label="BRAIN-Seq"
    )
    for mass in mass_s_transitions:
        mass_css_axis.vlines(mass, 0.975, 1, color="r", linestyle="-")
    style_axis(mass_css_axis, "Mass (Da)", "CSS", ylim=(None, 1))

    figure.tight_layout()
    plt.show()
