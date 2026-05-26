#!/usr/bin/env python

"""ADMIXTURE command line."""

import argparse
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 admixture -d snps.hdf5 -o OUT/ -k 2
$ ipyrad2 admixture -d snps.hdf5 -o OUT/ --k-range 2:5
$ ipyrad2 admixture -d snps.hdf5 -o OUT/ -k 3 --impute-method none --keep-intermediates
"""


def _setup_admixture_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 admixture` subcommand parser."""
    tool = subparsers.add_parser(
        "admixture",
        description=header,
        help="Run the external ADMIXTURE program on filtered SNP HDF5 data.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d", "--data", metavar="Path", type=Path, required=True,
        help="Path to an SNP-capable HDF5 file. Convert VCF first with `ipyrad2 vcf2hdf5`.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="admixture",
        help="Prefix name for output files. [default=admixture]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-admixture",
        help="Directory to write admixture outputs and stats. [default=output-admixture]",
    )

    clustering = tool.add_argument_group("Clustering")
    k_group = clustering.add_mutually_exclusive_group(required=True)
    k_group.add_argument(
        "-k", metavar="int", type=int,
        help="Number of clusters to fit.",
    )
    k_group.add_argument(
        "--k-range", metavar="MIN:MAX", type=str,
        help="Inclusive range of K values to score by ADMIXTURE cross validation.",
    )
    clustering.add_argument(
        "--binary", metavar="Path", type=Path,
        help="Optional path to the external `admixture` binary. Defaults to PATH lookup.",
    )
    clustering.add_argument(
        "--no-subsample", action="store_true",
        help="Keep linked SNPs. By default admixture subsamples one SNP per RAD locus.",
    )
    clustering.add_argument(
        "--seed", metavar="int", type=int,
        help="Random seed for SNP subsampling and optional imputation before PLINK staging.",
    )

    filtering = tool.add_argument_group("Filtering and samples")
    filtering.add_argument(
        "-m", "--min-sample-coverage", metavar="int", type=int, default=4,
        help="Minimum number of samples that must have data at a SNP. [default=4]",
    )
    filtering.add_argument(
        "-r", "--max-sample-missing", metavar="float", type=float, default=1.0,
        help="Maximum missing-data frequency allowed in a sample. [default=1.0]",
    )
    filtering.add_argument(
        "-a", "--min-minor-allele-frequency", metavar="float", type=float, default=0.0,
        help="Minimum minor allele frequency required to retain a SNP. [default=0.0]",
    )
    filtering.add_argument(
        "--min-genotype-depth", metavar="int", type=int, default=0,
        help="Mask sample genotypes with FORMAT/DP below this threshold before site filtering. [default=0]",
    )
    filtering.add_argument(
        "--min-site-qual", metavar="float", type=float, default=0.0,
        help="Minimum VCF QUAL score required to retain a SNP site. [default=0.0]",
    )
    filtering.add_argument(
        "-I", "--impute-method", metavar="str", choices=("sample", "none"), default="sample",
        help="Imputation method used before PLINK staging. [default=sample]",
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
        help="Sample-to-population mapping file with `sample<TAB>population` or `glob<TAB>population` on each line.",
    )
    filtering.add_argument(
        "-g", "--minmap", metavar="Path", type=Path,
        help=(
            "Population-to-minimum-coverage mapping file with "
            "`population<TAB>min` on each line. Adds per-population minimum "
            "coverage checks on top of `-m` when `imap` is used."
        ),
    )

    performance = tool.add_argument_group("Performance and overwrite")
    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=1,
        help="Number of threads to pass to ADMIXTURE and chunked SNP filtering. [default=1]",
    )
    performance.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep staged PLINK files and raw ADMIXTURE outputs instead of cleaning them up.",
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
