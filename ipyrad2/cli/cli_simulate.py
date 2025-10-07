#!/usr/bin/env python

import argparse
from pathlib import Path
from .make_wide import make_wide, intlike


EPILOG = """\
Examples
--------
$ ipyrad simulate -g 10000 -r1 A^TCGG -r2 T^TCAG -o READS/ -sg 123 -m 123 --tree
$ ipyrad simulate -g REF   -r1 ... -r2 ... -o READS         -sm 123 --tree 10
$ ipyrad simulate -g REF   -r1 ... -r2 ... -o READS         -sm 123 --tree NWK
$ ipyrad simulate -g 10000 -r1 ...         -o READS -sg 123 -sm 123
$ ipyrad simulate
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
        "-g", "--genome-size", metavar="int", type=int, required=True,
        help="Simulate a reference genome composed of random {A,C,T,G} of length '-g'."
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./TRIMMED",
        help="Directory to write results. Created if it doesn't exist. [default=SIM]",
    )
    tool.add_argument(
        "-n", "--name", metavar="str", type=str, default="simulated",
        help="Name prefix for output files [default='simulated']",
    )
    tool.add_argument(
        "-r1", "--re-pattern-1", metavar="str", type=str, default="AT^CGAT",
        help="Restriction cut pattern associated with index 1. Type II (symmetric cut) enzymes only. Mark the forward strand cut '^'. [default='AT^CGAT']",
    )
    tool.add_argument(
        "-r2", "--re-pattern-2", metavar="str", type=str, default="G^GATCC",
        help="Restriction cut pattern associated with index 2. Type II (symmetric cut) enzymes only. Mark the forward strand cut '^'. [default='G^GATCC']",
    )
    tool.add_argument(
        "-w", "--size-window", metavar="int", type=int, nargs=2, default=(300, 500),
        help="Digested fragment size selection window [default=(300, 500)]"
    )
    tool.add_argument(
        "-t", "--tree", metavar="str", type=str,
        help="Newick tree of samples to simulate data for [default=None=(,);]",
    )
    tool.add_argument(
        "-u", "--umi-tag-in-i5", action="store_true",
        help="Assign UMIs to i7 index attached to all reads (library indicator) [default='ATCGATCG']",
    )
    tool.add_argument(
        "--i7", metavar="str", type=str, default="ATCGATCG",
        help="i7 index attached to all reads (library indicator) [default='ATCGATCG']",
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
        "-x", "--max-reads", metavar="int", type=intlike, default=None,
        help="Maximum number of reads per file (useful for testing or to normalize inputs). [default=None]"
    )
    tool.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite .trimmed.fastq.gz files if they exist in outdir. Does not clear outdir.",
    )
    tool.add_argument(
        "-u", "--umi-tag-in-i5", action="store_true",
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
