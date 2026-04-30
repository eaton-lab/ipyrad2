#!/usr/bin/env python

"""Denovo command-line parser."""

import argparse
from argparse import ArgumentParser, Namespace
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = """\
Examples
--------
$ ipyrad2 denovo -d DATA/*.fastq.gz -o output-denovo
$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT -s 0.95 -S 0.85 -c 12 -t 3
$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT --imap denovo.imap.tsv
$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT --no-alignment
$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT -dx _R -di 1 --keep-intermediates
"""


def _parser_error_if(parser: ArgumentParser, condition: bool, message: str) -> None:
    """Raise an argparse error when a denovo CLI constraint is violated."""
    if condition:
        parser.error(message)


def validate_denovo_args(args: Namespace, parser: ArgumentParser) -> None:
    """Validate denovo CLI arguments after parsing."""
    _parser_error_if(parser, args.cores < 1, "--cores must be >= 1")
    _parser_error_if(parser, args.threads < 1, "--threads must be >= 1")
    _parser_error_if(parser, args.threads > args.cores, "--threads cannot exceed --cores")
    _parser_error_if(
        parser,
        not 0 < args.within_similarity <= 1,
        "--within-similarity must be > 0 and <= 1",
    )
    _parser_error_if(
        parser,
        not 0 < args.across_similarity <= 1,
        "--across-similarity must be > 0 and <= 1",
    )
    _parser_error_if(
        parser,
        not 0 < args.query_cov <= 1,
        "--query-cov must be > 0 and <= 1",
    )
    _parser_error_if(parser, args.min_derep_size < 1, "--min-derep-size must be >= 1")
    _parser_error_if(parser, args.min_length < 1, "--min-length must be >= 1")
    _parser_error_if(parser, args.min_merge_overlap < 1, "--min-merge-overlap must be >= 1")
    _parser_error_if(parser, args.max_merge_diffs < 0, "--max-merge-diffs must be >= 0")
    _parser_error_if(
        parser,
        args.delim_idx == 0,
        "--delim-idx cannot be 0",
    )


def _setup_denovo_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 denovo` subcommand parser."""
    tool = subparsers.add_parser(
        "denovo",
        description=header,
        help="Construct a denovo reference library by clustering reads and building locus consensuses.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    selection = tool.add_argument_group("Sample selection")
    clustering = tool.add_argument_group("Clustering and consensus")
    naming = tool.add_argument_group("Sample naming and library type")
    runtime = tool.add_argument_group("Runtime")
    logging = tool.add_argument_group("Logging")

    core.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="+",
        help="Input FASTQ files or glob patterns.",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./output-denovo",
        help="Output directory for the denovo reference library. [default=%(default)s]",
    )
    core.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite denovo outputs created by this command.",
    )

    selection.add_argument(
        "--imap", metavar="Path", type=Path,
        help="Optional IMAP file selecting denovo samples by `sample<TAB>group` or `glob<TAB>group`.",
    )
    selection.add_argument(
        "--use-all-samples", action="store_true",
        help="Disable automatic sample downselection and use every parsed input sample.",
    )

    clustering.add_argument(
        "-s", "--within-similarity", metavar="float", type=float, default=0.95,
        help="Sequence similarity threshold for clustering within samples. [default=%(default)s]",
    )
    clustering.add_argument(
        "-S", "--across-similarity", metavar="float", type=float, default=0.85,
        help="Sequence similarity threshold for clustering across samples. [default=%(default)s]",
    )
    clustering.add_argument(
        "-m", "--min-derep-size", metavar="int", type=int, default=5,
        help="Minimum duplicate count retained during dereplication. [default=%(default)s]",
    )
    clustering.add_argument(
        "-i", "--min-length", metavar="int", type=int, default=35,
        help="Minimum retained sequence length after pair merge or join. [default=%(default)s]",
    )
    clustering.add_argument(
        "-g", "--min-merge-overlap", metavar="int", type=int, default=20,
        help="Minimum overlap required to merge paired reads. [default=%(default)s]",
    )
    clustering.add_argument(
        "-e", "--max-merge-diffs", metavar="int", type=int, default=4,
        help="Maximum mismatches allowed in the merged region. [default=%(default)s]",
    )
    clustering.add_argument(
        "--query-cov", metavar="float", type=float, default=0.75,
        help="Minimum VSEARCH query coverage; lower for variable read lengths. [default=%(default)s]",
    )
    clustering.add_argument(
        "-b", "--allow-reverse-complement", action="store_true",
        help="Cluster both strands rather than plus strand only.",
    )
    clustering.add_argument(
        "--no-alignment", action="store_true",
        help="Skip MAFFT in the final locus step and use the longest stripped sequence per locus.",
    )

    naming.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Delimiter substring used to parse sample names from filenames.",
    )
    naming.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Delimiter index: positive from left, negative from right. [default=%(default)s]",
    )

    runtime.add_argument(
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Maximum total cores to use. [default=%(default)s]",
    )
    runtime.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Threads per vsearch or MAFFT job. [default=%(default)s]",
    )
    runtime.add_argument(
        "--keep-intermediates", action="store_true",
        help="Retain the denovo working directory instead of cleaning it on success.",
    )

    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
