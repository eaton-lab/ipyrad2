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
        help="Trim reads for quality, adapters, and restriction overhangs using 'fastp'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="One or more paths to fastq data files (or glob patterns; e.g., './data/*.fastq.gz')",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./TRIMMED",
        help="Directory to write results. Created if it doesn't exist. [%(default)s]",
    )
    tool.add_argument(
        "-u", "--max-unqualified-percent", metavar="int", type=int, default=15,
        help="Max percent of unqualified bases in a read. [%(default)s]",
    )
    tool.add_argument(
        "-q", "--min-quality", metavar="int", type=int, default=20,
        help="Min base quality. The threshold Q score below which a site is unqualified (see -). [%(default)s]",
    )
    tool.add_argument(
        "-M", "--min-mean-window-quality", metavar="int", type=int, default=30,
        help="Min mean base quality when trimming in sliding windows from 5' and 3'. [%(default)s]",
    )
    tool.add_argument(
        "-W", "--cut-window-size", metavar="int", type=int, default=5,
        help="Size of sliding windows when trimming if mean quality <= -m. [%(default)s]",
    )
    tool.add_argument(
        "-n", "--max-ns", metavar="int", type=int, default=5,
        help="Maximum number of uncalled (N) bases. [%(default)s]",
    )
    tool.add_argument(
        "-e", "--min-trimmed-length", metavar="int", type=int, default=35,
        help="Minimum length of retained trimmed reads. [%(default)s]",
    )
    tool.add_argument(
        "-x", "--max-reads", metavar="int", type=intlike, default=None,
        help="Maximum number of reads per file (useful for testing or to normalize inputs). [%(default)s]",
    )
    tool.add_argument(
        "-r", "--restriction-overhangs", metavar="str", nargs=2, type=str,
        help="Restriction overhangs (overrides inferred REs). [default=None None]",
    )
    tool.add_argument(
        "-p", "--phred-qscore-offset", metavar="int", type=int, default=33,
        help="Q score offset (to accommodate legacy data). [%(default)s]",
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
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Max number of cores to use. [%(default)s]",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Run c/t multi-threaded jobs concurrently. Larger -t reduces RAM and I/O. [%(default)s]",
    )
    tool.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite .trimmed.fastq.gz files if they exist in outdir. Does not clear outdir.",
    )
    tool.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Set delim substring 'dx' to override name parsing from files. [%(default)s]"
    )
    tool.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Set delim index. Extracts substring left of the 'di'-th 'dx' in filename. [%(default)s]",
    )
    # TODO
    tool.add_argument(
        "-s", "--suffix", metavar="str", type=str,
        help=r"Add a suffix to sample names (output: {name}{suffix}.trimmed._R*.fastq.gz). [%(default)s]",
    )
    tool.add_argument(
        "-U", "--umi-tag-in-i5", action="store_true",
        help="Move i5 index into read name. Use if i5 contains UMIs for marking PCR duplicates.",
    )
    tool.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    tool.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
