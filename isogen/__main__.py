"""Command-line entry points for ``python -m isogen``."""

import argparse
import sys

import numpy as np

from .isogen import isodist


def parse_input(value):
    """Return numeric mass inputs as floats and leave sequences as strings."""
    try:
        return float(value)
    except ValueError:
        return value


def build_parser():
    parser = argparse.ArgumentParser(
        prog="python -m isogen",
        description="Generate and plot protein, RNA, or DNA isotope distributions.",
        epilog=(
            "Examples:\n"
            "  python -m isogen dist 10000 --type PEPTIDE\n"
            "  python -m isogen dist ACDEFGHIK --type PEPTIDE --method NN\n"
            "  python -m isogen dist AUGCAGUACGUA --type RNA --output rna.csv\n"
            "  python -m isogen plot --save isodist_examples.png --no-show"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")

    dist_parser = subparsers.add_parser(
        "dist",
        help="generate a two-column mass/intensity distribution",
    )
    dist_parser.add_argument("input_value", metavar="INPUT", help="numeric mass or sequence")
    dist_parser.add_argument(
        "--type",
        default="PEPTIDE",
        choices=("PEPTIDE", "RNA", "DNA"),
        type=str.upper,
        help="analyte type (default: PEPTIDE)",
    )
    dist_parser.add_argument(
        "--method",
        default="FFT",
        choices=("FFT", "NN"),
        type=str.upper,
        help="distribution method (default: FFT)",
    )
    dist_parser.add_argument(
        "--isolen",
        default=128,
        type=int,
        help="number of isotope values (default: 128)",
    )
    dist_parser.add_argument(
        "-o",
        "--output",
        help="write CSV output to this path instead of standard output",
    )
    dist_parser.set_defaults(handler=run_distribution)

    plot_parser = subparsers.add_parser(
        "plot",
        help="show protein/RNA mass and sequence stick-plot examples",
    )
    plot_parser.add_argument(
        "--method",
        default="FFT",
        choices=("FFT", "NN"),
        type=str.upper,
        help="distribution method (default: FFT)",
    )
    plot_parser.add_argument(
        "--isolen",
        default=64,
        type=int,
        help="number of isotope values per plot (default: 64)",
    )
    plot_parser.add_argument("--save", metavar="PATH", help="save the figure to an image file")
    plot_parser.add_argument(
        "--dpi",
        default=150,
        type=int,
        help="saved image resolution (default: 150)",
    )
    plot_parser.add_argument(
        "--no-show",
        action="store_true",
        help="do not open an interactive plot window",
    )
    plot_parser.set_defaults(handler=run_plots)

    return parser


def run_distribution(args):
    if args.isolen < 0:
        raise ValueError("--isolen cannot be negative")

    distribution = isodist(
        input=parse_input(args.input_value),
        type=args.type,
        isolen=args.isolen,
        method=args.method,
    )
    destination = args.output if args.output else sys.stdout
    np.savetxt(
        destination,
        distribution,
        delimiter=",",
        header="mass,intensity",
        comments="",
    )
    return 0


def run_plots(args):
    if args.isolen < 0:
        raise ValueError("--isolen cannot be negative")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")

    import matplotlib.pyplot as plt
    from .plot_tests import plot_isodist_examples

    fig, _ = plot_isodist_examples(isolen=args.isolen, method=args.method)
    if args.save:
        fig.savefig(args.save, dpi=args.dpi, bbox_inches="tight")
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        return 0

    try:
        return args.handler(args)
    except (OSError, TypeError, ValueError) as exception:
        parser.error(str(exception))


if __name__ == "__main__":
    raise SystemExit(main())
