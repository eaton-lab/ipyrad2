#!/usr/bin/env python


import argparse
from pathlib import Path
from .make_wide import make_wide


EPILOG = """\
Examples
--------
$ ipyrad assemble --bams BAMs/*.bam --ref REF --out OUT -m 4 -d 5 -q 20
$ ipyrad assemble -b BAMs/a*.bam -B BAMS/b*.bam --ref REF --out OUT -m 4 -d 5 -q 20
"""


def _setup_assemble_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad assemble` subcommand parser.

    """
    tool = subparsers.add_parser(
        "assemble",
        description=header,
        help="assemble loci and call variants in shared mapping beds using 'bedtools' and 'bcftools'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-b", "--rad-bams", metavar="Path", type=Path, required=True, nargs="*",
        help="Bam files from RAD-type data. These samples are used to delimit locus beds. (regex allowed; e.g., './bam/{a,b}*.bam')",
    )
    tool.add_argument(
        "-B", "--wgs-bams", metavar="Path", type=Path, nargs="*",
        help="Optional bam files from WGS-type data. These samples are not used to delimit locus beds, but will have variants called within the RAD locus beds. (regex allowed; e.g., './bam/{a,b}*.bam')",
    )
    tool.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Path to the reference fasta used in the mapping step.",
    )
    tool.add_argument(
        "-n", "--name", metavar="str", type=str, default="assembly",
        help="Prefix name for output files. [default=assembly]",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./OUT",
        help="Directory for results and stat files. Created if it doesn't exist. [default=./OUT]",
    )
    tool.add_argument(
        "-q", "--min-gq", metavar="int", type=int, default=20,
        help="Min per-sample genotype quality score. [default=20]",
    )
    tool.add_argument(
        "-Q", "--min-qual", metavar="int", type=int, default=20,
        help="Min across-sample genotype quality score. [default=20]",
    )
    tool.add_argument(
        "-d", "--min-sample-depth", metavar="int", type=int, default=1,
        help="Min read depth within a sample to make variant calls. [default=1]",
    )
    tool.add_argument(
        "-m", "--min-locus-sample-coverage", metavar="int", type=int, default=4,
        help="Min num samples that must be present to retain a locus. [default=4]",
    )
    # This isn't super necessary. It reduces the size of the seqs h5 a bit,
    # but otherwise this filter is applied when you use wex.
    tool.add_argument(
        "-a", "--min-locus-trim-sample-coverage", metavar="int", type=int, default=4,
        help="Min num samples with non-N calls for trimming locus edges. Must be <= '-m'. [default=4]",
    )
    tool.add_argument(
        "-z", "--min-locus-length", metavar="int", type=int, default=25,
        help="Min length of locus after trimming. [default=25]",
    )
    # tool.add_argument(
    #     "-L", "--max-locus-length", metavar="int", type=int, default=None,
    #     help="Max length of locus (to prevent overlapping locus beds). [default=None]",
    # )
    tool.add_argument(
        "-g", "--min-locus-merge-distance", metavar="int", type=int, default=300,
        help="Merge locus beds if they overlap within nbases. [default=300]",
    )
    tool.add_argument(
        "-u", "--max-locus-hetero-frequency", metavar="float", type=float, default=0.3,
        help="Max frequency of samples heterozygous *at the same site* in a locus. [default=0.3]",
    )
    tool.add_argument(
        "-s", "--max-locus-variant-frequency", metavar="float", type=float, default=1.0,
        help="Max frequency of sites that are variant in a locus. [default=1.0]",
    )
    tool.add_argument(
        "-p", "--populations", metavar="Path", type=Path,
        help="Path to a population file where each line lists 'sample-name\tgroup-name'. [default=None]"
    )
    tool.add_argument(
        "-x", "--masks", metavar="str", nargs="*", type=str,
        help="Site patterns to mask (e.g., restriction overhangs). [default=None]",
    )
    tool.add_argument(
        "-e", "--exclude-reference", action="store_true",
        help="Do not include the reference sequence as a sample in outputs",
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
        help="Overwrite if out dir already contains result files with identical name.",
    )
    tool.add_argument(
        "-nx", "--name-delim", metavar="str", type=str, default=None,
        help="Set name delim substring 'nx' to override auto name parsing from files. [default=None]"
    )
    tool.add_argument(
        "-ni", "--name-index", metavar="int", type=int, default=1,
        help="Set name delim index to keep substring left of the 'ni'-th substring 'nx', if valid [default=1]",
    )
    tool.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    tool.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
