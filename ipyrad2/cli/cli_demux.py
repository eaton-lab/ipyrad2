#!/usr/bin/env python

"""
"""

import argparse
from pathlib import Path
from .make_wide import make_wide, intlike


EPILOG = """\
Examples
--------
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv -o DEMUX -c 10
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv --log-level DEBUG
$ ipyrad2 demux -d RUN1/*.fastq.gz RUN2/*.fastq.gz -b BARCODES.tsv -M
"""


def _setup_demux_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """parser for `ipyrad demux` cmd."""
    tool = subparsers.add_parser(
        "demux",
        description=header,
        help="Demultiplex pooled data to samples by index or barcode.",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140),
    )
    core = tool.add_argument_group("Core inputs")
    mode = tool.add_argument_group("Demultiplexing mode")
    overhangs = tool.add_argument_group("Restriction overhangs")
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

    overhangs.add_argument(
        "-re1", "--restriction-overhang-1", metavar="str", type=str,
        help="Restriction overhang on R1. Overrides kmer inference.",
    )
    overhangs.add_argument(
        "-re2", "--restriction-overhang-2", metavar="str", type=str,
        help="Restriction overhang on R2. Overrides kmer inference.",
    )
    overhangs.add_argument(
        "-R", "--disable-infer-re-overhangs", action="store_true",
        help="Skip kmer inference; use with explicit restriction overhangs.",
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
    # performance.add_argument(
    #     "-k", "--max-reads-kmer", metavar="int", type=intlike, default=500_000,
    #     help="Maximum number of reads sampled across files to infer REs from kmers. [default=5e5]",
    # )

    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Append logs to this file as well as stdout.",
    )
    # tool.add_argument(
    #     "--logger", type=str, nargs="*", default=("INFO", None),
    #     help=(
    #         "Logging info entered as one value for LOGLEVEL, or two values "
    #         "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    # )
    # TOO RISKY perhaps, make the user remove existing dir themselves?
    # tool.add_argument(
    #     "--force", "-f", action="store_true",
    #     help=(
    #         "Force overwrite. Allows overwriting an existing directory of "
    #         "demultiplexed fastq files. (Be careful).")
    # )
