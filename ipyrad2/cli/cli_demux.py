#!/usr/bin/env python

"""Demux command-line parser."""

import argparse
from argparse import ArgumentParser, Namespace
from pathlib import Path
from .common import RAW_HELP_FORMATTER, intlike


EPILOG = """\
Examples
--------
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv -o DEMUX -c 10
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv --log-level DEBUG
$ ipyrad2 demux -d RUN1/*.fastq.gz RUN2/*.fastq.gz -b BARCODES.tsv -M
"""


def _parser_error_if(parser: ArgumentParser, condition: bool, message: str) -> None:
    """Raise an argparse error when a demux CLI constraint is violated."""
    if condition:
        parser.error(message)


def validate_demux_args(args: Namespace, parser: ArgumentParser) -> None:
    """Validate demux CLI arguments after parsing."""
    _parser_error_if(parser, args.cores < 1, "--cores must be >= 1")
    _parser_error_if(parser, args.chunksize < 1, "--chunksize must be >= 1")
    _parser_error_if(
        parser,
        args.max_reads is not None and args.max_reads < 1,
        "--max-reads must be >= 1 when set",
    )
    _parser_error_if(parser, args.max_reads_kmer < 1, "--max-reads-kmer must be >= 1")
    _parser_error_if(
        parser,
        not 0 <= args.max_mismatch <= 2,
        "--max-mismatch must be between 0 and 2",
    )


def _setup_demux_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 demux` subcommand parser."""
    tool = subparsers.add_parser(
        "demux",
        description=header,
        help="Demultiplex pooled data to samples by index or barcode.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")
    mode = tool.add_argument_group("Demultiplexing mode")
    cutsites = tool.add_argument_group("Cutsite motifs")
    performance = tool.add_argument_group("Performance and sampling")
    logging = tool.add_argument_group("Logging")

    core.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="Input FASTQ files or glob patterns.",
    )
    core.add_argument(
        "-b", "--barcodes", metavar="Path", type=Path, required=True,
        help="Barcode/index table with sample barcode1 [barcode2] columns.",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./DEMUX",
        help="Output directory for demultiplexed FASTQs. Created if needed. [default=%(default)s]",
    )
    core.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite demux outputs from this run.",
    )

    mode.add_argument(
        "--i7", action="store_true",
        help="Demultiplex by i7 index instead of inline barcodes.",
    )
    mode.add_argument(
        "-m", "--max_mismatch", metavar="int", type=int, default=0,
        help="Allow up to N barcode/index mismatches. [default=%(default)s]",
    )
    mode.add_argument(
        "-M", "--merge-technical-replicates", action="store_true",
        help="Merge technical replicates that share a sample name.",
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
        "-E", "--disable-infer-cutsite-motifs", action="store_true",
        help="Skip cutsite motif inference; use with explicit cutsite motifs.",
    )

    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=4,
        help="Maximum parallel workers. [default=%(default)s]",
    )
    performance.add_argument(
        "-k", "--chunksize", metavar="int", type=intlike, default=int(1e7),
        help="Reads per write batch; larger values use more RAM. [default=%(default)s]",
    )
    performance.add_argument(
        "-x", "--max_reads", metavar="int", type=intlike,
        help="Read up to N reads per file for testing.",
    )
    performance.add_argument(
        "--max-reads-kmer", metavar="int", type=intlike, default=100_000,
        help="Total reads sampled across files for cutsite motif inference. [default=%(default)s]",
    )
    performance.add_argument(
        "--pigz", action="store_true",
        help="Use pigz for final demux FASTQ compression.",
    )

    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
    # tool.add_argument(
    #     "--logger", type=str, nargs="*", default=("INFO", None),
    #     help=(
    #         "Logging info entered as one value for LOGLEVEL, or two values "
    #         "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad2.txt.'")
    # )
