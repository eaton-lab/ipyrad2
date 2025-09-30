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
        "-o", "--out", metavar="Path", type=Path, default="./MAPPED",
        help="Directory to write trimmed read files. Will be created if it doesn't exist. [default=MAPPED]",
    )
    # tool.add_argument(
    #     "-i", "--imap", metavar="Path", type=Path,
    #     help="Optional file with sample\tgroup to assign samples to read groups.",
    # )
    tool.add_argument(
        "-m", "--mark-dups-by-coords", action="store_true",
        help="Mark PCR duplicates by coordinates. Only use for WGS data.",
    )
    tool.add_argument(
        "-u", "--mark-dups-by-umis", action="store_true",
        help="Mark PCR duplicates by UMIs. Only use with i5 tags (see ipyrad trim -u).",
    )
    tool.add_argument(
        "-w", "--workers", metavar="int", type=int, default=2,
        help="N concurrent workers (jobs) to parallelize. [default=4]",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=4,
        help="N threads per worker (e.g., -w 2 -t 4 uses up to 8 threads). [default=2]",
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
