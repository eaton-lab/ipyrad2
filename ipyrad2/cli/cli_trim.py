#!/usr/bin/env python

"""Trim command-line parser."""

import argparse
from argparse import ArgumentParser, Namespace
from pathlib import Path
from .common import RAW_HELP_FORMATTER, intlike


EPILOG = """\
Examples
--------
$ ipyrad2 trim -d DATA/*.fastq.gz -o TRIMMED/
$ ipyrad2 trim -d DATA/*.gz -o OUT -q 20 -u 15 -M 30 -W 5 -n 5 -e 35 -c 12 -t 3
$ ipyrad2 trim -d DATA/*.gz -o OUT -dx _R -di 1
$ ipyrad2 trim -d DATA/*.gz -o OUT --phred64 -U
"""


def _parser_error_if(parser: ArgumentParser, condition: bool, message: str) -> None:
    """Raise an argparse error when a trim CLI constraint is violated."""
    if condition:
        parser.error(message)


def validate_trim_args(args: Namespace, parser: ArgumentParser) -> None:
    """Validate numeric trim CLI arguments after parsing."""
    _parser_error_if(parser, args.cores < 1, "--cores must be >= 1")
    _parser_error_if(parser, args.threads < 1, "--threads must be >= 1")
    _parser_error_if(parser, args.threads > args.cores, "--threads cannot exceed --cores")
    _parser_error_if(
        parser,
        not 0 <= args.max_unqualified_percent <= 100,
        "--max-unqualified-percent must be between 0 and 100",
    )
    _parser_error_if(parser, args.min_quality < 0, "--min-quality must be >= 0")
    _parser_error_if(
        parser,
        not 1 <= args.min_mean_window_quality <= 36,
        "--min-mean-window-quality must be between 1 and 36",
    )
    _parser_error_if(
        parser,
        not 1 <= args.cut_window_size <= 1000,
        "--cut-window-size must be between 1 and 1000",
    )
    _parser_error_if(parser, args.max_ns < 0, "--max-ns must be >= 0")
    _parser_error_if(parser, args.min_trimmed_length < 1, "--min-trimmed-length must be >= 1")
    _parser_error_if(
        parser,
        args.max_reads is not None and args.max_reads < 1,
        "--max-reads must be >= 1 when set",
    )
    _parser_error_if(parser, args.max_reads_kmer < 1, "--max-reads-kmer must be >= 1")


def _setup_trim_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 trim` subcommand parser."""
    tool = subparsers.add_parser(
        "trim",
        description=header,
        help="Trim reads for quality, adapters, and cutsite motifs using 'fastp'.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")
    filtering = tool.add_argument_group("Filtering and trimming")
    cutsites = tool.add_argument_group("Cutsite motifs and adapters")
    performance = tool.add_argument_group("Performance and compatibility")
    naming = tool.add_argument_group("Sample naming and UMI")
    logging = tool.add_argument_group("Logging")

    core.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="Input FASTQ files or glob patterns.",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./TRIMMED",
        help="Output directory for trimmed FASTQs. Created if needed. [default=%(default)s]",
    )
    core.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite trim outputs from this run.",
    )

    filtering.add_argument(
        "-q", "--min-quality", metavar="int", type=int, default=20,
        help="Minimum base quality used to mark bases as unqualified. [default=%(default)s]",
    )
    filtering.add_argument(
        "-u", "--max-unqualified-percent", metavar="int", type=int, default=15,
        help="Maximum percent of unqualified bases allowed per read. [default=%(default)s]",
    )
    filtering.add_argument(
        "-M", "--min-mean-window-quality", metavar="int", type=int, default=30,
        help="Minimum sliding-window mean quality. [default=%(default)s]",
    )
    filtering.add_argument(
        "-W", "--cut-window-size", metavar="int", type=int, default=5,
        help="Sliding-window size for front and tail trimming. [default=%(default)s]",
    )
    filtering.add_argument(
        "-n", "--max-ns", metavar="int", type=int, default=5,
        help="Maximum number of N bases allowed per read. [default=%(default)s]",
    )
    filtering.add_argument(
        "-e", "--min-trimmed-length", metavar="int", type=int, default=35,
        help="Minimum retained read length after trimming. [default=%(default)s]",
    )
    filtering.add_argument(
        "-Q", "--disable-quality-filtering", action="store_true",
        help="Skip fastp quality filtering.",
    )

    cutsites.add_argument(
        "-e1", "--cutsite-1", metavar="str", type=str,
        help="5' restriction-site remnant / cutsite motif at the start of R1. Use commas for multiple motifs; overrides inference.",
    )
    cutsites.add_argument(
        "-e2", "--cutsite-2", metavar="str", type=str,
        help="5' restriction-site remnant / cutsite motif at the start of R2. Use commas for multiple motifs; overrides inference.",
    )
    cutsites.add_argument(
        "-k", "--max-reads-kmer", metavar="int", type=intlike, default=500_000,
        help="Total reads sampled across files for cutsite motif inference. [default=500000]",
    )
    cutsites.add_argument(
        "-E", "--disable-infer-cutsite-motifs", action="store_true",
        help="Skip cutsite motif inference.",
    )
    cutsites.add_argument(
        "-A", "--disable-adapter-trimming", action="store_true",
        help="Skip adapter trimming.",
    )

    performance.add_argument(
        "-x", "--max-reads", metavar="int", type=intlike, default=None,
        help="Read up to N reads per file for testing or normalization.",
    )
    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Maximum total cores to use. [default=%(default)s]",
    )
    performance.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Threads per fastp job. [default=%(default)s]",
    )
    performance.add_argument(
        "--phred64", action="store_true",
        help="Treat input qualities as legacy phred64 and convert to phred33.",
    )

    naming.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Delimiter substring used to parse sample names from filenames.",
    )
    naming.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Keep text left of the Nth delimiter when parsing sample names. [default=%(default)s]",
    )
    naming.add_argument(
        "-s", "--suffix", metavar="str", type=str,
        help="Suffix appended to parsed sample names before writing outputs.",
    )
    naming.add_argument(
        "-U", "--umi-tag-in-i5", action="store_true",
        help="Move the i5 index into the read name as a UMI tag.",
    )

    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
