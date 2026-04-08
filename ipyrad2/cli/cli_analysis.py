#!/usr/bin/env python

"""Analysis command-line entrypoints."""

import argparse
import sys

from loguru import logger

from .command_log import format_logged_command
from .common import RAW_HELP_FORMATTER
from .cli_wex import _setup_wex_subparser
from .cli_lex import _setup_lex_subparser
from .cli_snpex import _setup_snpex_subparser
from .cli_vcf_to_hdf5 import _setup_vcf_to_hdf5_subparser
from .cli_pca import _setup_pca_subparser
from .cli_snmf import _setup_snmf_subparser
from .cli_dapc import _setup_dapc_subparser
from .cli_admixture import _setup_admixture_subparser
from .cli_popgen import _setup_popgen_subparser
from ..analysis.extracters.window_extracter import run_window_extracter
from ..analysis.extracters.locus_extracter import run_locus_extracter
from ..analysis.extracters.snps_extracter import run_snps_extracter
from ..analysis.converters.vcf_to_hdf5 import run_vcf_to_hdf5
from ..analysis.methods.pca import run_pca_method
from ..analysis.methods.snmf import run_snmf_method
from ..analysis.methods.dapc import run_dapc_method
from ..analysis.methods.admixture import run_admixture_method
from ..analysis.methods.popgen import run_popgen_method
from ipyrad2 import __version__ as VERSION


HEADER = f"""
-------------------------------------------------------------
ipyrad2 [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = "ipyrad2 analysis command line. Select a positional subcommand:"


EPILOG = r"""
Examples
--------
$ ipyrad2 analysis wex -d assembly.hdf5 --print-scaffold-table
$ ipyrad2 analysis wex -d assembly.hdf5 -o OUT/ -n TEST -m 10 -w Chr1
$ ipyrad2 analysis lex -d assembly.hdf5 -o OUT/ -N 100 -L 150
$ ipyrad2 analysis snpex -d assembly.hdf5 -o SNP_OUT/ -n SNPSET
$ ipyrad2 analysis vcf-to-hdf5 -d variants.vcf.gz -o SNP_HDF5/
$ ipyrad2 analysis pca -d snps.hdf5 -o PCA_OUT/
$ ipyrad2 analysis snmf -d snps.hdf5 -o SNMF_OUT/ -k 2
$ ipyrad2 analysis dapc -d snps.hdf5 -o DAPC_OUT/ --k-range 2:5
$ ipyrad2 analysis admixture -d snps.hdf5 -o ADMIX_OUT/ --k-range 2:5
$ ipyrad2 analysis popgen -d assembly.hdf5 -o POPGEN_OUT/
"""


def _setup_analysis_subparser(subparser: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 analysis` subcommand parser."""
    analysis_parser = subparser.add_parser(
        "analysis",
        #description=f"{HEADER}\n{DESCRIPTION}",
        description=header,
        help="Utilities for downstream analysis",
        formatter_class=RAW_HELP_FORMATTER)

    analysis_subparser = analysis_parser.add_subparsers(
        dest="tool",
        required=True)
    _setup_wex_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis wex: extract one alignment from selected genomic windows",
    )
    _setup_lex_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis lex: extract delimited loci from HDF5 database",
    )
    _setup_snpex_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis snpex: extract filtered SNP matrices from HDF5 database",
    )
    _setup_vcf_to_hdf5_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis vcf-to-hdf5: convert VCF to SNP-capable HDF5 database",
    )
    _setup_pca_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis pca: run PCA, t-SNE, or UMAP on SNP HDF5 data",
    )
    _setup_snmf_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis snmf: run sNMF-style clustering on SNP HDF5 data",
    )
    _setup_dapc_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis dapc: run DAPC-style clustering on SNP HDF5 data",
    )
    _setup_admixture_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis admixture: run external ADMIXTURE on SNP HDF5 data",
    )
    _setup_popgen_subparser(
        analysis_subparser,
        f"{HEADER}\nipyrad2 analysis popgen: compute genome-wide population-genetic statistics",
    )


def run_analysis_tool(args):

    # WEX: --------------------------------------------------------
    if args.tool == "wex":
        logger.info("-------------------------------------------------------")
        logger.info("----- ipyrad2 analysis wex: extract alignments from windows -----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_window_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            out_format=args.out_format,
            windows=args.windows,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            include_reference=args.include_reference,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            print_scaffold_table=args.print_scaffold_table,
            stdout=args.stdout,
            force=args.force,
        )
        sys.exit(0)

    if args.tool == "lex":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 analysis lex: extract delimited loci from HDF5 database ----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_locus_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            out_format=args.out_format,
            nloci=args.nloci,
            min_length=args.min_length,
            windows=args.windows,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            include_reference=args.include_reference,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            print_scaffold_table=args.print_scaffold_table,
            stdout=args.stdout,
            force=args.force,
        )

        sys.exit(0)

    if args.tool == "snpex":
        logger.info("------------------------------------------------------------")
        logger.info("---- ipyrad2 analysis snpex: extract filtered SNP matrices ----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_snps_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            include_reference=args.include_reference,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
            subsample=not args.no_subsample,
            random_seed=args.seed,
            write_plink=args.plink,
            write_phylip=args.phylip,
            write_nexus=args.nexus,
            write_fasta=args.fasta,
            write_treemix=args.treemix,
            write_eems=args.eems,
            impute_method=args.impute_method,
        )
        sys.exit(0)

    if args.tool == "vcf-to-hdf5":
        logger.info("-----------------------------------------------------------")
        logger.info("---- ipyrad2 analysis vcf-to-hdf5: convert VCF to HDF5 ----")
        logger.info("-----------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_vcf_to_hdf5(
            data=args.data,
            name=args.name,
            outdir=args.out,
            ld_block_size=args.ld_block_size,
            force=args.force,
        )
        sys.exit(0)

    if args.tool == "pca":
        logger.info("----------------------------------------------------------")
        logger.info("---- ipyrad2 analysis pca: numerical PCA-family methods ----")
        logger.info("----------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_pca_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            method=args.method,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            impute_method=args.impute_method,
            subsample=not args.no_subsample,
            random_seed=args.seed,
            replicates=args.replicates,
            perplexity=args.perplexity,
            max_iter=args.max_iter,
            n_neighbors=args.n_neighbors,
            plot=args.plot,
            plot_width=args.plot_width,
            plot_height=args.plot_height,
            plot_marker_size=args.plot_marker_size,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
        )
        sys.exit(0)

    if args.tool == "snmf":
        logger.info("------------------------------------------------------------")
        logger.info("---- ipyrad2 analysis snmf: numerical clustering on SNPs ----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_snmf_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            k=args.k,
            k_range=args.k_range,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            impute_method=args.impute_method,
            subsample=not args.no_subsample,
            random_seed=args.seed,
            cores=args.cores,
            force=args.force,
            alpha_w=args.alpha_w,
            alpha_h=args.alpha_h,
            l1_ratio=args.l1_ratio,
            n_init=args.n_init,
            cv_replicates=args.cv_replicates,
            cv_holdout=args.cv_holdout,
            log_level=args.log_level,
        )
        sys.exit(0)

    if args.tool == "dapc":
        logger.info("------------------------------------------------------------")
        logger.info("---- ipyrad2 analysis dapc: numerical clustering on SNPs ----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_dapc_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            k=args.k,
            k_range=args.k_range,
            n_pcs=args.n_pcs,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            impute_method=args.impute_method,
            subsample=not args.no_subsample,
            random_seed=args.seed,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
        )
        sys.exit(0)

    if args.tool == "admixture":
        logger.info("----------------------------------------------------------------")
        logger.info("---- ipyrad2 analysis admixture: external clustering on SNPs ----")
        logger.info("----------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_admixture_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            k=args.k,
            k_range=args.k_range,
            binary=None if args.binary is None else str(args.binary),
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            impute_method=args.impute_method,
            subsample=not args.no_subsample,
            random_seed=args.seed,
            keep_intermediates=args.keep_intermediates,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
        )
        sys.exit(0)

    if args.tool == "popgen":
        logger.info("------------------------------------------------------------")
        logger.info("---- ipyrad2 analysis popgen: genome-wide population genetics ----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_popgen_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            stats=args.stats,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            min_minor_allele_frequency=args.min_minor_allele_frequency,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            subsample_unlinked=args.subsample_unlinked,
            random_seed=args.seed,
            window_size=args.window_size,
            step_size=args.step_size,
            loci_per_window=args.loci_per_window,
            locus_step=args.locus_step,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
        )
        sys.exit(0)
