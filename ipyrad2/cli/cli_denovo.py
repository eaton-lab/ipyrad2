#!/usr/bin/env python


import argparse
from pathlib import Path
from .make_wide import make_wide


EPILOG = """\
Examples
--------
$ ipyrad denovo --fastqs DATA/*.fastq.gz --out OUT/ -c 0.95 -C 0.85
"""


def _setup_denovo_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad assemble` subcommand parser.

    """
    tool = subparsers.add_parser(
        "denovo",
        description=header,
        help="construct a denovo reference library by clustering with 'vsearch'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="One or more paths to fastq data files (regex allowed; e.g., './data/*.fastq.gz')",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./CLUSTERS",
        help="Directory to write trimmed read files. Will be created if it doesn't exist. [default=CLUSTERS]",
    )
    tool.add_argument(
        "-s", "--similarity-threshold-within", metavar="float", type=float, default=0.95,
        help="Sequence similarity threshold for clustering within samples. [default=0.95]",
    )
    tool.add_argument(
        "-S", "--similarity-threshold-across", metavar="float", type=float, default=0.85,
        help="Sequence similarity threshold for clustering across samples. [default=0.85]",
    )
    tool.add_argument(
        "-m", "--min-dereplication-size", metavar="int", type=int, default=2,
        help="Min replication of sequences for inclusion. [default=2]",
    )
    tool.add_argument(
        "-i", "--min-length", metavar="int", type=int, default=35,
        help="Min length of merged paired sequences. [default=35]",
    )
    tool.add_argument(
        "-g", "--min-merge-overlap", metavar="int", type=int, default=20,
        help="Min overlap to merge paired reads. [default=20]",
    )
    tool.add_argument(
        "-e", "--max-merge-diffs", metavar="int", type=int, default=4,
        help="Max difference in merged region. [default=4]",
    )
    tool.add_argument(
        "-b", "--strand-both", action="store_true",
        help="Match both strands. Use for single-enzyme GBS.",
    )
    tool.add_argument(
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Max number of cores to use. [default=6]",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Run c/t multi-threaded jobs concurrently. Larger -t reduces RAM and I/O. [default=3]",
    )
    tool.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite if out dir contains fastq file with identical name.",
    )
    tool.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Set delim substring 'dx' to override name parsing from files. [default=None]"
    )
    tool.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Set delim index. Extracts substring left of the 'di'-th 'dx' in filename. [default=1]",
    )
    tool.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    tool.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
