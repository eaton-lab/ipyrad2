"""Command-line parser for the experimental locus exporter."""

from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 seqex -d assembly.hdf5 -o OUT/
$ ipyrad2 seqex -d assembly.hdf5 -w Chr1 -N 100 -s 123 -c 4 -O nex
$ ipyrad2 seqex -d assembly.hdf5 -L 150 -C -O fa
$ ipyrad2 seqex -d assembly.hdf5 -i imap.tsv -a -X
"""


def _setup_seqex_subparser(subparsers, header=None):
    """Add the top-level ipyrad2 seqex parser."""
    tool = subparsers.add_parser(
        "seqex",
        description=header,
        help="Extract filtered, delimited loci from an HDF5 database.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d",
        "--data",
        metavar="Path",
        type=Path,
        required=True,
        help="Assembly HDF5 file containing sequence and locus-map datasets.",
    )
    core.add_argument(
        "-n",
        "--name",
        metavar="str",
        default="alignment",
        help="Prefix for alignment and stats files. [default=alignment]",
    )
    core.add_argument(
        "-o",
        "--out",
        metavar="Path",
        type=Path,
        default="output-seqex",
        help="Output directory. [default=output-seqex]",
    )
    core.add_argument(
        "-O",
        "--out-format",
        metavar="str",
        choices=["phy", "nex", "fa"],
        default="phy",
        help="Output format: PHYLIP, NEXUS, or FASTA. [default=phy]",
    )

    sampling = tool.add_argument_group("Locus sampling")
    sampling.add_argument(
        "-w",
        "--windows",
        metavar="str",
        nargs="*",
        help="Select complete loci overlapping scaffold names, regexes, regions, or one BED file.",
    )
    sampling.add_argument(
        "-N",
        "--max-loci",
        metavar="int",
        type=int,
        help="Randomly retain at most this many loci after filtering. [default=all]",
    )
    sampling.add_argument(
        "-s",
        "--random-seed",
        metavar="int",
        type=int,
        help="Non-negative seed for reproducible --max-loci sampling.",
    )
    sampling.add_argument(
        "-L",
        "--min-length",
        metavar="int",
        type=int,
        help="Minimum locus length after site filtering. [default=disabled]",
    )

    filtering = tool.add_argument_group("Filtering and samples")
    filtering.add_argument(
        "-m",
        "--min-sample-coverage",
        metavar="int",
        type=int,
        default=4,
        help="Minimum samples represented in each locus and retained site. [default=4]",
    )
    filtering.add_argument(
        "-r",
        "--max-sample-missing",
        metavar="float",
        type=float,
        default=1.0,
        help="Maximum missing fraction per locus and, with -C, across the final matrix. [default=1.0]",
    )
    filtering.add_argument(
        "-e",
        "--exclude",
        metavar="str",
        nargs="*",
        help="Exclude samples by name; takes precedence over IMAP membership and -R.",
    )
    filtering.add_argument(
        "-R",
        "--include-reference",
        action="store_true",
        help="Include assembly_reference_sequence.",
    )
    filtering.add_argument(
        "-i",
        "--imap",
        metavar="Path",
        type=Path,
        help="Sample-to-population mapping file.",
    )
    filtering.add_argument(
        "-g",
        "--minmap",
        metavar="Path",
        type=Path,
        help="Population-to-minimum-coverage mapping; overrides -m.",
    )
    filtering.add_argument(
        "-c",
        "--cores",
        metavar="int",
        type=int,
        default=1,
        help="Number of processes used to filter locus batches. [default=1]",
    )

    output = tool.add_argument_group("Output control")
    layout = output.add_mutually_exclusive_group()
    layout.add_argument(
        "-C",
        "--concatenate",
        action="store_true",
        help="Concatenate accepted loci into one alignment matrix.",
    )
    layout.add_argument(
        "-X",
        "--split",
        action="store_true",
        help="Write every accepted locus to a separate alignment file.",
    )
    output.add_argument(
        "-a",
        "--append-population",
        action="store_true",
        help="Write names as sample^population; requires --imap.",
    )
    output.add_argument(
        "-P",
        "--print-scaffold-table",
        action="store_true",
        help="Print the scaffold table to stdout and exit.",
    )
    output.add_argument(
        "-x",
        "--stdout",
        action="store_true",
        help="Write the alignment to stdout instead of a file.",
    )
    output.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing output files with identical names.",
    )

    logging = tool.add_argument_group("Logging")
    logging.add_argument(
        "-l",
        "--log-level",
        metavar="str",
        default="INFO",
        help="Log level. [default=INFO]",
    )
    logging.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )
