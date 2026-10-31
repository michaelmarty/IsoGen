"""Plot the POPC (C42H82NO8P) isotope distribution with 0% and 5% 2H."""

import matplotlib.pyplot as plt
import numpy as np

import isogen

FORMULA = "C42H82NO8P"
ISOLEN = 10
PERCENT = 0.01


def main():
    natural = isogen.isodist(FORMULA, type="ATOM", isolen=ISOLEN)
    labeled = isogen.isodist(FORMULA, type="ATOM", isolen=ISOLEN, heavy_isotopes={"2H": PERCENT})

    fig, ax = plt.subplots(figsize=(8, 4.5))
    width = 0.35
    index = np.arange(ISOLEN)
    ax.bar(index - width / 2, natural[:, 1], width, label="0% 2H (natural)")
    ax.bar(index + width / 2, labeled[:, 1], width, label=f"{PERCENT*100:.0f}% 2H")
    ax.set_xlabel(f"Mass offset from monoisotopic ({natural[0, 0]:.4f} Da)")
    ax.set_ylabel("Relative intensity")
    ax.set_title(f"POPC ({FORMULA}) isotope distribution")
    ax.legend()
    fig.tight_layout()
    # fig.savefig("popc_deuterium.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
