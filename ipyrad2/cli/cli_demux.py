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


# The validator and schema use the same constraints beyond argparse choices.
DEMUX_DEPENDENCIES = (
    ("tmpdir", "requires", "pigz", "--tmpdir requires --pigz"),
    (
        "allow_leading_barcode_deletion", "conflicts_with", "i7",
        "--allow-leading-barcode-deletion applies only to inline barcode demux",
    ),
)
DEMUX_VALUE_CONSTRAINTS = {
    "cores": ({"minimum": 1}, "--cores must be >= 1"),
    "chunksize": ({"minimum": 1}, "--chunksize must be >= 1"),
    "max_reads": ({"minimum": 1}, "--max-reads must be >= 1 when set"),
    "max_reads_kmer": ({"minimum": 1}, "--max-reads-kmer must be >= 1"),
    "max_mismatch": ({"minimum": 0, "maximum": 2}, "--max-mismatch must be between 0 and 2"),
    "delim_idx": ({"excluded_values": [0]}, "--delim-idx cannot be 0"),
}


def validate_demux_args(args: Namespace, parser: ArgumentParser) -> None:
    """Validate demux CLI arguments after parsing."""
    for name, relation, other, message in DEMUX_DEPENDENCIES:
        if getattr(args, name) and bool(getattr(args, other)) != (relation == "requires"):
            parser.error(message)

    for name, (limits, message) in DEMUX_VALUE_CONSTRAINTS.items():
        value = getattr(args, name)
        if value is None:
            continue
        if (
            value < limits.get("minimum", value)
            or value > limits.get("maximum", value)
            or value in limits.get("excluded_values", ())
        ):
            parser.error(message)


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
    naming = tool.add_argument_group("Sample naming and pairing")
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
    mode.add_argument(
        "--barcode-boundary-slack", metavar="int", type=int, choices=(0, 1), default=1,
        help="Position offset allowed for inline barcode matching [default=%(default)s]"
    )
    mode.add_argument(
        "--allow-leading-barcode-deletion", action="store_true",
        help=(
            "Recover inline barcode observations missing the first expected barcode base. "
            "Incompatible with --i7. Disabled by default."
        ),
    )

    naming.add_argument(
        "-dx", "--delim-str", metavar="str", type=str, default=None,
        help="Delimiter substring used to parse sample names from filenames.",
    )
    naming.add_argument(
        "-di", "--delim-idx", metavar="int", type=int, default=1,
        help="Delimiter index: positive from left, negative from right. [default=%(default)s]",
    )

    cutsites.add_argument(
        "-e1", "--cutsite-1", metavar="str", type=str,
        help="5' cutsite motif at the start of R1. Use commas for multiple motifs",
    )
    cutsites.add_argument(
        "-e2", "--cutsite-2", metavar="str", type=str,
        help="5' cutsite motif at the start of R2. Use commas for multiple motifs",
    )
    cutsites.add_argument(
        "-E", "--disable-infer-cutsite-motifs", action="store_true",
        help="Skip cutsite motif inference; explicit motifs are compared to inferred motifs unless -E is set.",
    )

    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=4,
        help=(
            "Core budget for demultiplexing and, with --pigz, for compression "
            "thread allocation. [default=%(default)s]"
        ),
    )
    performance.add_argument(
        "-k", "--chunksize", metavar="int", type=intlike, default=int(1e7),
        help="Maximum matched reads per batch; batches are also capped at 4 MiB. [default=%(default)s]",
    )
    performance.add_argument(
        "-x", "--max_reads", metavar="int", type=intlike,
        help=(
            "Process the first N input read pairs per paired-end FASTQ pair, or the "
            "first N reads per single-end FASTQ file. Counts input before barcode "
            "matching. Defaults to no limit."
        ),
    )
    performance.add_argument(
        "--max-reads-kmer", metavar="int", type=intlike, default=100_000,
        help="Total reads sampled across files for cutsite motif inference. [default=%(default)s]",
    )
    performance.add_argument(
        "--pigz", action="store_true",
        help="Compress temporary FASTQs with pigz (required on PATH); local scratch (--tmpdir) is recommended.",
    )

    performance.add_argument(
        "--tmpdir", metavar="Path", type=Path,
        help="Existing writable scratch directory. Requires --pigz; defaults to private output staging.",
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
