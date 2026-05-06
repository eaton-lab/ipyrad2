#!/usr/bin/env python

"""ABBA/BABA admixture analysis command line."""

import argparse
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 baba -d assembly.hdf5 -o OUT/ --tests quartets.tsv
$ ipyrad2 baba -d assembly.hdf5 -o OUT/ --tree species.nwk --f-branch
$ ipyrad2 baba -d assembly.hdf5 -o OUT/ -i IMAP.txt -g MINMAP.txt --tests quartets.tsv
$ ipyrad2 baba -d snps.hdf5 -o OUT/ --tree species.nwk --resampling bootstrap --bootstrap-replicates 2000 --seed 7
"""


def _setup_baba_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 baba` subcommand parser."""
    tool = subparsers.add_parser(
        "baba",
        description=header,
        help="Compute ABBA/BABA quartet statistics and related admixture summaries from SNP HDF5 data.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d", "--data", metavar="Path", type=Path, required=True,
        help="Input SNP HDF5 file containing `genos` and `snpsmap`.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="baba",
        help="Output file prefix. [default=baba]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-baba",
        help="Output directory. [default=output-baba]",
    )
    quartet_source = core.add_mutually_exclusive_group(required=True)
    quartet_source.add_argument(
        "--tests", metavar="Path", type=Path,
        help="Whitespace-delimited quartet file with one `P1 P2 P3 P4` test per line.",
    )
    quartet_source.add_argument(
        "--tree", metavar="Path", type=Path,
        help="Rooted binary Newick tree used to expand valid quartets automatically.",
    )

    filtering = tool.add_argument_group("Filtering and samples")
    filtering.add_argument(
        "-i", "--imap", metavar="Path", type=Path,
        help="Sample-to-population mapping file with `sample<TAB>population` entries.",
    )
    filtering.add_argument(
        "-g", "--minmap", metavar="Path", type=Path,
        help="Population-to-minimum-coverage file; applied per quartet and requires `--imap`.",
    )
    filtering.add_argument(
        "-m", "--min-sample-coverage", metavar="int", type=int, default=1,
        help="Minimum selected samples with data before quartet-specific filters. [default=1]",
    )
    filtering.add_argument(
        "-e", "--exclude", metavar="str", type=str, nargs="*",
        help="Exclude one or more samples by name. This overrides IMAP membership and `-R`.",
    )
    filtering.add_argument(
        "-R", "--include-reference", action="store_true",
        help="Include `assembly_reference_sequence` in the selected sample set.",
    )

    statistics = tool.add_argument_group("Statistics and resampling")
    statistics.add_argument(
        "--resampling",
        metavar="str",
        choices=("auto", "jackknife", "bootstrap", "none"),
        default="auto",
        help="Significance method; `auto` prefers physical-block jackknife. [default=auto]",
    )
    statistics.add_argument(
        "--bootstrap-replicates", metavar="int", type=int, default=1000,
        help="Bootstrap replicate count. [default=1000]",
    )
    statistics.add_argument(
        "--jackknife-block-bp", metavar="int", type=int, default=5_000_000,
        help="Physical block size in bp for jackknife sampling. [default=5000000]",
    )
    statistics.add_argument(
        "--jackknife-block-loci", metavar="int", type=int, default=100,
        help="Fallback loci per block for sparse or denovo-like jackknife runs. [default=100]",
    )
    statistics.add_argument(
        "--seed", metavar="int", type=int,
        help="Random seed for resampling.",
    )

    outputs = tool.add_argument_group("Optional outputs")
    outputs.add_argument(
        "--f-branch", action="store_true",
        help="Write tree-derived branch summary tables. Requires `--tree`.",
    )
    outputs.add_argument(
        "--f-branch-p-threshold", metavar="float", type=float, default=0.01,
        help="P-value threshold for zeroing non-significant tree f_G values. [default=0.01]",
    )
    outputs.add_argument(
        "--write-block-table", action="store_true",
        help="Write per-block site-pattern summaries for each quartet.",
    )
    outputs.add_argument(
        "--clustering-stats", action="store_true",
        help="Add quartet-level summaries of spatial clustering across inferred blocks.",
    )

    performance = tool.add_argument_group("Performance and overwrite")
    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=1,
        help="Number of cores used during SNP filtering. [default=1]",
    )
    performance.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite existing output files.",
    )

    logging = tool.add_argument_group("Logging")
    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (TRACE, DEBUG, INFO, WARNING, ERROR) [default=INFO]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
