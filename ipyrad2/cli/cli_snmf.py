#!/usr/bin/env python

"""sNMF-style clustering command line."""

import argparse
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 snmf -d snps.hdf5 -o OUT/ -k 2
$ ipyrad2 snmf -d snps.hdf5 -o OUT/ --k-range 2:5
$ ipyrad2 snmf -d snps.hdf5 -o OUT/ -k 3 --impute-method none
$ ipyrad2 snmf -d snps.hdf5 -o OUT/ --k-range 2:6 --alpha-w 1e-3 --n-init 20
"""


def _setup_snmf_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 snmf` subcommand parser."""
    tool = subparsers.add_parser(
        "snmf",
        description=header,
        help="Run sklearn-backed sNMF-style clustering on SNP HDF5 data.",
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
        "-n", "--name", metavar="str", type=str, default="snmf",
        help="Prefix name for output files. [default=snmf]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-snmf",
        help="Directory to write numerical outputs and stats. [default=output-snmf]",
    )

    clustering = tool.add_argument_group("Clustering")
    k_group = clustering.add_mutually_exclusive_group(required=True)
    k_group.add_argument(
        "-k", metavar="int", type=int,
        help="Number of clusters to fit.",
    )
    k_group.add_argument(
        "--k-range", metavar="MIN:MAX", type=str,
        help="Inclusive range of K values to score by masked-genotype cross-entropy.",
    )
    clustering.add_argument(
        "--no-subsample", action="store_true",
        help="Keep linked SNPs. By default snmf subsamples one SNP per RAD locus.",
    )
    clustering.add_argument(
        "--seed", metavar="int", type=int,
        help="Random seed for SNP subsampling, imputation, CV masking, and NMF initialization.",
    )

    regularization = tool.add_argument_group("Regularization and scoring")
    regularization.add_argument(
        "--alpha-w", metavar="float", type=float, default=1e-4,
        help="L1 regularization strength on ancestry coefficients W. [default=1e-4]",
    )
    regularization.add_argument(
        "--alpha-h", metavar="float|same", type=str, default="same",
        help="Regularization strength on genotype-loadings H, or 'same' for sklearn defaults. [default=same]",
    )
    regularization.add_argument(
        "--l1-ratio", metavar="float", type=float, default=1.0,
        help="Elastic-net mixing for NMF regularization. Use 1.0 for pure L1 sparsity. [default=1.0]",
    )
    regularization.add_argument(
        "--n-init", metavar="int", type=int, default=10,
        help="Number of NMF initializations to try per K. The best reconstruction is kept. [default=10]",
    )
    regularization.add_argument(
        "--cv-replicates", metavar="int", type=int, default=5,
        help="Number of masked-genotype cross-entropy replicates used to score each K. [default=5]",
    )
    regularization.add_argument(
        "--cv-holdout", metavar="float", type=float, default=0.1,
        help="Fraction of observed genotypes to mask per cross-entropy replicate. [default=0.1]",
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
        help="Imputation method for the numerical genotype matrix. [default=sample]",
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
