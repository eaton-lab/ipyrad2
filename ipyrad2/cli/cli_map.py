#!/usr/bin/env python

"""Map command-line parser."""

import argparse
from pathlib import Path
from .common import RAW_HELP_FORMATTER


EPILOG = """\
Examples
--------
$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS
$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS -i IMAP.tsv
$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS --reindex-reference
$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS -m
$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS -u
"""


def _setup_map_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 map` subcommand parser."""
    tool = subparsers.add_parser(
        "map",
        description=header,
        help="Map reads with 'bwa-mem2' and 'samtools' and write coordinate-sorted BAMs.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")
    duplicates = tool.add_argument_group("Duplicate removal")
    naming = tool.add_argument_group("Sample naming and grouping")
    performance = tool.add_argument_group("Performance and overwrite")
    logging = tool.add_argument_group("Logging")

    core.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="Input FASTQ files or glob patterns.",
    )
    core.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Reference FASTA to map against. Build bwa-mem2 indexes if missing, otherwise reuse existing indexes.",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./MAPPED",
        help="Output directory for coordinate-sorted BAMs and map stats. Created if needed. [default=%(default)s]",
    )
    core.add_argument(
        "--unmate", action="store_true",
        help="Treat paired FASTQs as single-end by concatenating R1 and R2 per sample before mapping.",
    )

    duplicates.add_argument(
        "-m", "--mark-dups-by-coords", action="store_true",
        help="Remove PCR duplicates by coordinates; intended for WGS data.",
    )
    duplicates.add_argument(
        "-u", "--mark-dups-by-umis", action="store_true",
        help="Remove PCR duplicates by UMI tags from `ipyrad2 trim -U`.",
    )

    naming.add_argument(
        "-i", "--imap", metavar="Path", type=Path,
        help="Sample-to-group table for subsetting, renaming, or merging samples.",
    )
    naming.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Delimiter substring used to parse sample names from filenames.",
    )
    naming.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Delimiter index: positive from left, negative from right. [default=%(default)s]",
    )

    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Maximum total cores to use. [default=%(default)s]",
    )
    performance.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Threads per mapping job; larger values trade concurrency for lower I/O overhead. [default=%(default)s]",
    )
    performance.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite existing BAM outputs for matching sample names. Does not re-index the reference.",
    )
    performance.add_argument(
        "--reindex-reference", action="store_true",
        help="Rebuild bwa-mem2 reference indexes even when matching sidecar files already exist.",
    )

    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
