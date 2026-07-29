import matplotlib.pyplot as plt
from isogen import isodist

test_mass = 10000
test_pep_seq = "MKTVVLAVAVLFLTGSQARHFWQRDDPQTPWDRVKDFATVYVDAVKDSGREYVSQFETSALGKQLNLNLLENWDTLGSTVGRLQEQLGPVTQEFWDNLEKETEW" \
     "LRREMNKDLEEVKAKVQPYLDQFQTKWQEEVALYRQKMEPLGAELRDGARQKLQELQEKLTPLGEDLRDRMRHHVDALRTKMTPYSDQMRDRLAERLAQLKDSPTL" \
     "AEYHTKAADHLKAFGEKAKPALEDLRQGLMPVFESFKTRIMSMVEEASKKLNAQ"
test_rna_seq = "AUGCAGUACGUA"


def stick_plot(ax, distribution, title):
    """Plot a two-column mass/intensity distribution as sticks."""
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
    """Plot protein and RNA isotope distributions from mass and sequence inputs."""
    examples = [
        ("Protein — mass input", test_mass, "PEPTIDE"),
        ("Protein — sequence input", test_pep_seq, "PEPTIDE"),
        ("RNA — mass input", test_mass, "RNA"),
        ("RNA — sequence input", test_rna_seq, "RNA"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (title, input_value, analyte_type) in zip(axes.flat, examples):
        distribution = isodist(
            input=input_value,
            type=analyte_type,
            isolen=isolen,
            method=method,
        )
        stick_plot(ax, distribution, title)

    fig.suptitle("{} IsoGen isotope distributions".format(method))
    fig.tight_layout()
    return fig, axes


if __name__ == "__main__":
    plot_isodist_examples()
    plt.show()
