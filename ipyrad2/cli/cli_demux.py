#!/usr/bin/env python

"""
"""

import argparse
from pathlib import Path
from .make_wide import make_wide, intlike


EPILOG = """\
Examples
--------
$ ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv -c 10 -o ./demux
$ ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --i7
$ ipyrad demux -d RAW/*.fastq.gz -b BARCODES.csv --logger DEBUG
$ ipyrad demux -d RAWS1/*.fastq.gz RAWS2/*.fastqs.gz -b BARCODES.tsv
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
    tool.add_argument(
        "-d", "--fastqs", metavar="Path", type=Path, required=True, nargs="*",
        help="One or more paths to fastq data files (regex allowed; e.g., './data/*.fastq.gz')",
    )
    tool.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./DEMUX",
        help="Directory to write results. Created if it doesn't exist. [default=DEMUX]",
    )
    tool.add_argument(
        "-b", "--barcodes", metavar="Path", type=Path, required=True,
        help="Barcode/index file (CSV, TSV, or whitespace delimited) containing name\tindex or name\tindex1\tindex2"
    )
    tool.add_argument(
        "-re1", "--restriction-overhang-1", metavar="str", type=str,
        help="Restriction overhang (junction motif) on read1s. If None it is inferred by kmer analysis."
    )
    tool.add_argument(
        "-re2", "--restriction-overhang-2", metavar="str", type=str,
        help="Restriction overhang (junction motif) on read2s. If None it is inferred by kmer analysis."
    )
    tool.add_argument(
        "-m", "--max_mismatch", metavar="int", type=int, default=0,
        help=("The max number of allowed mismatches between true and "
        "oberved index/barcode [default=0].")
    )
    tool.add_argument(
        "-x", "--max_reads", metavar="int", type=intlike,
        help="Sample only the first N reads per file. Useful for testing. [default=None]."
    )
    tool.add_argument(
        "-k", "--chunksize", metavar="int", type=intlike, default=int(1e7),
        help=("N reads to process between writing to disk. "
        "Larger values = faster, but uses more RAM. [default=1e7]")
    )
    # tool.add_argument(
    #     "-k", "--max-reads-kmer", metavar="int", type=intlike, default=500_000,
    #     help="Maximum number of reads sampled across files to infer REs from kmers. [default=5e5]",
    # )
    tool.add_argument(
        "-R", "--disable-infer-re-overhangs", action="store_true",
        help="Disable automated inference of restriction overhangs using kmer analysis.",
    )
    tool.add_argument(
        "-M", "--merge-technical-replicates", action="store_true",
        help=(
        "If the same sample name is associated with multiple indices in the barcode file "
        "the default behavior is to append '-technical-replicate-x' to each name to make "
        "them unique, unless this option is turned on, in which case they are merged.")
    )
    tool.add_argument(
        "--i7",action="store_true",
        help="Demultiplex on i7 index instead of inline barcodes."
    )
    tool.add_argument(
        "-c", "--cores", metavar="int", type=int, default=4,
        help="Max number of cores. [default=4]",
    )
    tool.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    tool.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
    # tool.add_argument(
    #     "--logger", type=str, nargs="*", default=("INFO", None),
    #     help=(
    #         "Logging info entered as one value for LOGLEVEL, or two values "
    #         "for LOGLEVEL LOGFILE; e.g., 'DEBUG' or 'DEBUG ipyrad.txt.'")
    # )

    # TOO RISKY perhaps, make the user remove existing dir themselves?
    tool.add_argument(
        "--force", "-f", action="count", default=0,
        help=(
            "Force overwrite. Allows overwriting demultiplexed fastq files "
            "in the output directory. (Be careful).")
    )
