import matplotlib.pyplot as plt

if __package__:
    from . import isodist
else:
    from isogen import isodist

test_mass = 10000
test_pep_seq = "MKTVVLAVAVLFLTGSQARHFWQRDDPQTPWDRVKDFATVYVDAVKDSGREYVSQFETSALGKQLNLNLLENWDTLGSTVGRLQEQLGPVTQEFWDNLEKETEW" \
     "LRREMNKDLEEVKAKVQPYLDQFQTKWQEEVALYRQKMEPLGAELRDGARQKLQELQEKLTPLGEDLRDRMRHHVDALRTKMTPYSDQMRDRLAERLAQLKDSPTL" \
     "AEYHTKAADHLKAFGEKAKPALEDLRQGLMPVFESFKTRIMSMVEEASKKLNAQ"
test_rna_seq = "AUGCAGUACGUA"
test_dna_seq = "ATGCAGTACGTAATGCAGTACGTAATGCAGTACGTAATGCAGTACGTA"
test_atomic_formula = "Pd2H4O2"


def stick_plot(ax, distribution, title):
    """Draw a mass/intensity array as a stick plot.

    Args:
        ax: Matplotlib axes on which to draw.
        distribution: Two-column mass/intensity NumPy array.
        title: Axes title.
    """
    ax.stem(
        distribution[:, 0],
        distribution[:, 1],
        linefmt="C0-",
        markerfmt=" ",
        basefmt=" ",
    )
    ax.set_title(title)
    ax.set_xlabel("Mass (Da)")
    ax.set_ylabel("Relative intensity")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.25)


def plot_isodist_examples(isolen=64, method="FFT"):
    """Create protein, RNA, DNA, and formula example plots.

    Args:
        isolen: Number of isotope positions in each subplot.
        method: Distribution engine, either ``FFT`` or ``NN``, for the
            protein and RNA example panels.

    Returns:
        ``(figure, axes)`` from a 3-by-2 Matplotlib subplot layout.
    """
    examples = [
        ("Protein — mass input ({})".format(method), test_mass, "PEPTIDE", method),
        ("Protein — sequence input ({})".format(method), test_pep_seq, "PEPTIDE", method),
        ("RNA — mass input ({})".format(method), test_mass, "RNA", method),
        ("RNA — sequence input ({})".format(method), test_rna_seq, "RNA", method),
        ("DNA — sequence input (NN)", test_dna_seq, "DNA", "NN"),
        ("Formula — atomic input (FFT)", test_atomic_formula, "ATOM", "FFT"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    for ax, (title, input_value, analyte_type, example_method) in zip(axes.flat, examples):
        distribution = isodist(
            input=input_value,
            type=analyte_type,
            isolen=isolen,
            method=example_method,
        )
        stick_plot(ax, distribution, title)

    fig.suptitle("IsoGen isotope distributions")
    fig.tight_layout()
    return fig, axes


if __name__ == "__main__":
    plot_isodist_examples()
    plt.show()
