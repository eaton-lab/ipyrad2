#!/usr/bin/env python


import argparse
from pathlib import Path
from .make_wide import make_wide, intlike


EPILOG = """\
Examples
--------
$ ipyrad trim --fastqs DATA/*.fastq.gz --out OUT/
$ ipyrad trim -d DATA/*.gz -o OUT -q 20 -n 5 -m 30 -x 1e5 -c 12 -t 4
$ ipyrad trim -d DATA/*.gz -o OUT -nx _ -ni 1
$ ipyrad trim -d DATA/*.gz -o OUT -nx=-R -ni -1
"""


def _setup_trim_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad assemble` subcommand parser.

    """
    tool = subparsers.add_parser(
        "trim",
        description=header,
        help="trim reads for quality, adapters, and restriction overhangs using 'fastp'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="One or more paths to fastq data files (regex allowed; e.g., './data/*.fastq.gz')",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./trimmed",
        help="Directory to write results. Created if it doesn't exist. [default=./trimmed]",
    )
    tool.add_argument(
        "-q", "--min-quality", metavar="int", type=int, default=20,
        help="Minimum base quality score. [default=20]",
    )
    tool.add_argument(
        "-n", "--max-low-quality-bases", metavar="int", type=int, default=5,
        help="Maximum number of low quality bases. [default=5]",
    )
    tool.add_argument(
        "-m", "--min-trimmed-length", metavar="int", type=int, default=35,
        help="Minimum length of retained trimmed reads. [default=35]",
    )
    tool.add_argument(
        "-x", "--max-reads", metavar="int", type=intlike, default=35,
        help="Maximum number of reads per file (useful for testing or to normalize inputs). [default=None]"
    )
    tool.add_argument(
        "-r", "--restriction-overhangs", metavar="str", nargs=2, type=str,
        help="Restriction overhangs (overrides inferred REs). [default=(None, None)]",
    )
    tool.add_argument(
        "-p", "--phred-qscore-offset", metavar="int", type=int, default=33,
        help="Q score offset (to accommodate legacy data). [default=33]",
    )
    tool.add_argument(
        "-k", "--max-reads-kmer", metavar="int", type=intlike, default=500_000,
        help="Maximum number of reads sampled across files to infer REs from kmers. [default=5e5]",
    )
    tool.add_argument(
        "-R", "--disable-infer-re-overhangs", action="store_true",
        help="Disable infer restriction overhangs using kmer analysis.",
    )
    tool.add_argument(
        "-Q", "--disable-quality-filtering", action="store_true",
        help="Disable quality filtering.",
    )
    tool.add_argument(
        "-A", "--disable-adapter-trimming", action="store_true",
        help="Disable adapter trimming.",
    )
    tool.add_argument(
        "-c", "--cores", metavar="int", type=int, default=4,
        help="Number of cores available for processing. [default=4]",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=2,
        help="Number of threads (e.g., -c 4 -t 2 will run 2 2-threaded jobs). [default=2]",
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
