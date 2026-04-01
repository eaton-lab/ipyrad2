#!/usr/bin/env python

"""PCA-family numerical analysis command line."""

import argparse
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ --plot
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ --plot --plot-width 520 --plot-height 360
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ -M tsne --perplexity 8
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ -M umap --n-neighbors 10
$ ipyrad2 analysis pca -d snps.hdf5 -o OUT/ --no-subsample --impute-method zero
"""


def _setup_pca_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 analysis pca` subcommand parser."""
    tool = subparsers.add_parser(
        "pca",
        description=header,
        help="Run PCA, t-SNE, or UMAP on filtered SNP HDF5 data.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d", "--data", metavar="Path", type=Path, required=True,
        help="Path to an SNP-capable HDF5 file. Convert VCF first with `analysis vcf-to-hdf5`.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="pca",
        help="Prefix name for output files. [default=pca]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-pca",
        help="Directory to write numerical outputs and stats. [default=output-pca]",
    )

    method = tool.add_argument_group("Method and linkage")
    method.add_argument(
        "-M", "--method", metavar="str", choices=("pca", "tsne", "umap"), default="pca",
        help="Method to run: pca, tsne, or umap. [default=pca]",
    )
    method.add_argument(
        "--replicates", metavar="int", type=int, default=1,
        help="Number of PCA replicate runs. Only valid with `-M pca`. [default=1]",
    )
    method.add_argument(
        "--no-subsample", action="store_true",
        help="Keep linked SNPs. By default pca subsamples one SNP per RAD locus.",
    )
    method.add_argument(
        "--seed", metavar="int", type=int,
        help="Random seed for SNP subsampling, imputation, and method initialization.",
    )
    method.add_argument(
        "--perplexity", metavar="float", type=float, default=5.0,
        help="t-SNE perplexity. Used only with `-M tsne`. [default=5.0]",
    )
    method.add_argument(
        "--max-iter", metavar="int", type=int, default=1000,
        help="t-SNE maximum iterations. Used only with `-M tsne`. [default=1000]",
    )
    method.add_argument(
        "--n-neighbors", metavar="int", type=int, default=15,
        help="UMAP neighbor count. Used only with `-M umap`. [default=15]",
    )

    plotting = tool.add_argument_group("Plotting")
    plotting.add_argument(
        "--plot", action="store_true",
        help="Write an SVG plot for PCA results. Only supported with `-M pca`.",
    )
    plotting.add_argument(
        "--plot-width", metavar="int", type=int, default=400,
        help="SVG width in pixels for `--plot`. [default=400]",
    )
    plotting.add_argument(
        "--plot-height", metavar="int", type=int, default=300,
        help="SVG height in pixels for `--plot`. [default=300]",
    )
    plotting.add_argument(
        "--plot-marker-size", metavar="int", type=int, default=10,
        help="Marker size for `--plot`. [default=10]",
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
        "-I", "--impute-method", metavar="str", choices=("sample", "zero"), default="sample",
        help="Impute missing genotypes with `sample` or `zero`. [default=sample]",
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
