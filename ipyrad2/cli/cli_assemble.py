#!/usr/bin/env python


import argparse
from pathlib import Path
from .make_wide import make_wide


EPILOG = """\
Examples
--------
$ ipyrad assemble -d BAMs/r*.bam --ref REF --out OUT -m 4 -s 5 -q 20
$ ipyrad assemble -d BAMs/r*.bam -w BAMS/w*.bam --ref REF --out OUT -m 4 -s 5 -q 20
$ ipyrad assemble -w BAMS/w*.bam --ref REF -b loci.bed --out OUT -m 4 -q 20
"""


def _setup_assemble_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad assemble` subcommand parser.

    """
    tool = subparsers.add_parser(
        "assemble",
        description=header,
        help="Assemble loci and call variants using 'bedtools' and 'bcftools'.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    tool.add_argument(
        "-d", "--rad-bams", metavar="Path", type=Path, required=True, nargs="*",
        help="Bam files from RAD-type samples. (glob supported; e.g., './bam/{a,b}*.bam'). "
             "These data are used to delimit loci regions (unless overruled by -b), and assembled",
    )
    tool.add_argument(
        "-w", "--wgs-bams", metavar="Path", type=Path, nargs="*",
        help="Optional bam files from WGS-type data. (glob supported; e.g., './bam/{a,b}*.bam') "
             "These data are only assembled within loci regions delimited by RAD samples (or set using -b)"
    )
    tool.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Path to the reference fasta used in the mapping step",
    )
    tool.add_argument(
        "-b", "--loci-bed", metavar="Path", type=Path,
        help="Optional bed file delimiting loci on the reference genome.",
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
        "-qs", "--min-site-q", metavar="int", type=int, default=20,
        help="Min variant site quality (QUAL: 'confidence a site is variant'). [default=20]."
    )
    tool.add_argument(
        "-qg", "--min-geno-q", metavar="int", type=int, default=13,
        help="Min genotype quality (GQ: 'confidence in a sample's genotype). [default=13]",
    )
    tool.add_argument(
        "-qb", "--min-base-q", metavar="int", type=int, default=13,
        help="Min base quality score (BQ: 'confidence in base call'). [default=13]",
    )
    # tool.add_argument(
    #     "-q", "--min-map-q", metavar="int", type=int, default=10,
    #     help="Min alignment quality... allow user only to apply this in mapper step.
    # )
    tool.add_argument(
        "-s", "--min-sample-depth", metavar="int", type=int, default=1,
        help="Min read depth within a sample to make variant calls. [default=1]",
    )
    tool.add_argument(
        "-m", "--min-locus-sample-coverage", metavar="int", type=int, default=4,
        help="Min num samples that must be present to retain a locus. [default=4]",
    )
    # This isn't super necessary. It reduces the size of the seqs h5 a bit,
    # but otherwise this filter is applied when you use wex. Meh, let's keep it,
    # it makes loci edges looks nicer.
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
        "-y", "--max-locus-variant-frequency", metavar="float", type=float, default=1.0,
        help="Max frequency of sites that are variant in a locus. [default=1.0]",
    )
    tool.add_argument(
        "-p", "--populations", metavar="Path", type=Path,
        help=r"Pop file ('name\tpop' lines) to group samples for joint variant calls. [default=None]"
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
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Max number of cores to use. [default=6]",
    )
    tool.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Run c/t multi-threaded jobs concurrently. Larger -t reduces RAM and I/O. [default=3]",
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
