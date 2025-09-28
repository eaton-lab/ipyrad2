#!/usr/bin/env python


import argparse
from pathlib import Path
from .make_wide import make_wide, intlike


EPILOG = """\
Examples
--------
$ ipyrad map --fastqs DATA/*.fastq.gz --ref REF --out OUT/
"""


def _setup_map_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad assemble` subcommand parser.

    """
    tool = subparsers.add_parser(
        "map",
        description=header,
        help="map, filter, and sort reads with 'bwa-mem2' and 'samtools'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="One or more paths to fastq data files (regex allowed; e.g., './data/*.fastq.gz')",
    )
    tool.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Directory to write trimmed read files. Will be created if it doesn't exist.",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./trimmed",
        help="Directory to write trimmed read files. Will be created if it doesn't exist.",
    )
    # tool.add_argument(
    #     "-i", "--imap", metavar="Path", type=Path,
    #     help="Optional file with sample\tgroup to assign samples to read groups.",
    # )
    tool.add_argument(
        "-c", "--cores", metavar="int", type=int, default=4,
        help="Number of cores available for processing.",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=2,
        help="Number of threads (e.g., -c 4 -t 2 will run 2 2-threaded jobs).",
    )
    tool.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite if out dir contains fastq file with identical name.",
    )
    tool.add_argument(
        "-nx", "--name-delim", metavar="str", type=str, default=None,
        help="Set delim substring 'nx' to override auto name parsing from files. [default=None]"
    )
    tool.add_argument(
        "-ni", "--name-index", metavar="int", type=int, default=1,
        help="Set delim index to split file names to keep substring left of the "
        "'ni'-th occurrence of substring 'nx', if valid. [default=1]",
    )
    tool.add_argument(
        "--logger", type=str, nargs="*",
        help=(
            "Logging info entered as one value for LOGLEVEL, or two values "
            "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    )
