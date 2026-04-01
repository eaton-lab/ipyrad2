#!/usr/bin/env python

"""Population-genetic summary command line."""

import argparse
from pathlib import Path

from ..analysis.methods.popgen.models import ORDERED_STATS
from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 analysis popgen -d assembly.hdf5 -o OUT/
$ ipyrad2 analysis popgen -d assembly.hdf5 -o OUT/ --stats pi,dxy,fst,tajima_d,theta_w,fis,fit
$ ipyrad2 analysis popgen -d assembly.hdf5 -o OUT/ --stats pi,fst,fis,fit --window-size 100000 --step-size 50000
$ ipyrad2 analysis popgen -d assembly.hdf5 -o OUT/ --stats pi,fst --loci-per-window 25 --locus-step 10
$ ipyrad2 analysis popgen -d snps.hdf5 -o OUT/ --stats fst,heterozygosity,fis,fit,sfs --subsample-unlinked --seed 7
$ ipyrad2 analysis popgen -d assembly.hdf5 -o OUT/ -i POPS.txt -g MINMAP.txt
"""

SUPPORTED_STATS = ",".join(ORDERED_STATS)


def _setup_popgen_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 analysis popgen` subcommand parser."""
    tool = subparsers.add_parser(
        "popgen",
        description=header,
        help="Compute genome-wide population-genetic statistics from analysis HDF5 data.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d", "--data", metavar="Path", type=Path, required=True,
        help="Path to an analysis HDF5 file. Sequence HDF5 supports the full panel; SNP HDF5 supports the SNP-backed subset.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="popgen",
        help="Prefix name for output files. [default=popgen]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-popgen",
        help="Directory to write popgen outputs and stats. [default=output-popgen]",
    )
    core.add_argument(
        "--stats", metavar="str", type=str, default="all",
        help=(
            "Comma-separated statistics to compute. Supported names are "
            f"{SUPPORTED_STATS}. Use `all` for the full panel supported by the "
            "detected input. [default=all]"
        ),
    )

    filtering = tool.add_argument_group("Filtering and samples")
    filtering.add_argument(
        "-m", "--min-sample-coverage", metavar="int", type=int, default=1,
        help="Minimum number of samples that must have data at a site. [default=1]",
    )
    filtering.add_argument(
        "-r", "--max-sample-missing", metavar="float", type=float, default=1.0,
        help="Maximum missing-data frequency allowed in a sample. [default=1.0]",
    )
    filtering.add_argument(
        "-a", "--min-minor-allele-frequency", metavar="float", type=float, default=0.0,
        help="Minimum minor allele frequency for SNP-backed runs only. [default=0.0]",
    )
    filtering.add_argument(
        "-e", "--exclude", metavar="str", type=str, nargs="*",
        help="Exclude one or more samples by name. This takes precedence over IMAP membership and `-R`.",
    )
    filtering.add_argument(
        "-R", "--include-reference", action="store_true",
        help="Include `assembly_reference_sequence`. By default it is excluded unless IMAP already contains it.",
    )
    filtering.add_argument(
        "-i", "--imap", metavar="Path", type=Path,
        help="Sample-to-population mapping file with `sample<TAB>population` on each line.",
    )
    filtering.add_argument(
        "-g", "--minmap", metavar="Path", type=Path,
        help=(
            "Population-to-minimum-coverage mapping file with "
            "`population<TAB>min` on each line. Used with `imap` for "
            "population-aware coverage filtering; SNP-backed runs also keep "
            "the global `-m` filter."
        ),
    )

    linkage = tool.add_argument_group("SNP-backed options")
    linkage.add_argument(
        "--subsample-unlinked", action="store_true",
        help="For SNP-backed runs only, subsample one SNP per RAD locus instead of using all filtered SNPs.",
    )
    linkage.add_argument(
        "--seed", metavar="int", type=int,
        help="Random seed for SNP-backed unlinked subsampling.",
    )

    windowing = tool.add_argument_group("Windowing")
    window_mode = windowing.add_mutually_exclusive_group()
    window_mode.add_argument(
        "--window-size", metavar="int", type=int,
        help="Sequence-backed runs only. Compute additional scaffold windows of this size in bp.",
    )
    window_mode.add_argument(
        "--loci-per-window", metavar="int", type=int,
        help="Sequence-backed runs only. Compute additional anonymous RAD windows of this many consecutive loci.",
    )
    windowing.add_argument(
        "--step-size", metavar="int", type=int,
        help="Step size in bp for `--window-size`. Defaults to the same value as `--window-size`.",
    )
    windowing.add_argument(
        "--locus-step", metavar="int", type=int,
        help="Step size in loci for `--loci-per-window`. Defaults to the same value as `--loci-per-window`.",
    )

    performance = tool.add_argument_group("Performance and overwrite")
    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=1,
        help="Number of cores to use during chunked SNP filtering. [default=1]",
    )
    performance.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite existing output files with identical names.",
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
