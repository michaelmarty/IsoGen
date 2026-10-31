"""Plot the POPC (C42H82NO8P) isotope distribution with 0% and 1% 2H.

The labeled distribution is then fitted with ``calc_heavy_percent`` and the
recovered enrichment and fit time are shown in the plot title.
"""

import time

import matplotlib.pyplot as plt
import numpy as np

import isogen
import matplotlib as mpl

mpl.use("WxAgg")

FORMULA = "C42H82NO8P"
ISOLEN = 10
PERCENT = 0.01


def fit_enrichment(labeled):
    """Fit the 2H percent of a labeled distribution; returns (result, seconds)."""
    start = time.perf_counter()
    result = isogen.calc_heavy_percent(labeled[:, 1], FORMULA, "2H")
    return result, time.perf_counter() - start


def main(show=True):
    natural = isogen.isodist(FORMULA, type="ATOM", isolen=ISOLEN)
    labeled = isogen.isodist(FORMULA, type="ATOM", isolen=ISOLEN, heavy_isotopes={"2H": PERCENT})
    result, seconds = fit_enrichment(labeled)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    width = 0.35
    index = np.arange(ISOLEN)
    ax.bar(index - width / 2, natural[:, 1], width, label="0% 2H (natural)")
    ax.bar(index + width / 2, labeled[:, 1], width, label=f"{PERCENT*100:.0f}% 2H")
    ax.plot(index + width / 2, result.fit_dist, "k.", label="fit")
    ax.set_xlabel(f"Mass offset from monoisotopic ({natural[0, 0]:.4f} Da)")
    ax.set_ylabel("Relative intensity")
    ax.set_title(f"POPC ({FORMULA}) fitted 2H = {result.percent:.4f}% "
                 f"\u00b1 {result.stderr:.2g} (true {PERCENT*100:.4g}%), {seconds*1000:.1f} ms")
    ax.legend()
    fig.tight_layout()
    # fig.savefig("popc_deuterium.png", dpi=150)
    if show:
        plt.show()
    return fig, result, seconds


if __name__ == "__main__":
    main()
