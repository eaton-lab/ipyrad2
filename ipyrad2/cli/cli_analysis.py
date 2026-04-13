#!/usr/bin/env python

"""Top-level export and analysis command registration and dispatch."""

from __future__ import annotations

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
from .cli_bpp import _setup_bpp_subparser
from ..analysis.extracters.window_extracter import run_window_extracter
from ..analysis.extracters.locus_extracter import run_locus_extracter
from ..analysis.extracters.snps_extracter import run_snps_extracter
from ..analysis.converters.vcf_to_hdf5 import run_vcf_to_hdf5
from ..analysis.methods.pca import run_pca_method
from ..analysis.methods.snmf import run_snmf_method
from ..analysis.methods.dapc import run_dapc_method
from ..analysis.methods.admixture import run_admixture_method
from ..analysis.methods.popgen import run_popgen_method
from ..analysis.methods.bpp import run_bpp_method
from ..utils.exceptions import IPyradError
from ipyrad2 import __version__ as VERSION


HEADER = f"""
-------------------------------------------------------------
ipyrad2 [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

ANALYSIS_TOOL_NAMES = (
    "wex",
    "lex",
    "snpex",
    "vcf2hdf5",
    "pca",
    "snmf",
    "dapc",
    "admixture",
    "popgen",
    "bpp",
)

RESERVED_TOOL_NAMES = (
    "baba",
    "treeslider",
)


def _setup_reserved_tool_subparser(
    subparsers: argparse._SubParsersAction,
    *,
    name: str,
    help_text: str,
    header: str,
) -> None:
    """Add one reserved top-level command placeholder."""
    tool = subparsers.add_parser(
        name,
        description=header,
        help=help_text,
        epilog=(
            "Reserved command placeholder.\n"
            "This command is not implemented yet."
        ),
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    logging = tool.add_argument_group("Logging")
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )


def _setup_analysis_tool_subparsers(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add top-level export and analysis subcommand parsers."""
    header = HEADER if header is None else header
    _setup_wex_subparser(
        subparsers,
        f"{header}\nipyrad2 wex: extract one alignment from selected genomic windows",
    )
    _setup_lex_subparser(
        subparsers,
        f"{header}\nipyrad2 lex: extract delimited loci from HDF5 database",
    )
    _setup_snpex_subparser(
        subparsers,
        f"{header}\nipyrad2 snpex: extract filtered SNP matrices from HDF5 database",
    )
    _setup_vcf_to_hdf5_subparser(
        subparsers,
        f"{header}\nipyrad2 vcf2hdf5: convert VCF to SNP-capable HDF5 database",
    )
    _setup_pca_subparser(
        subparsers,
        f"{header}\nipyrad2 pca: run PCA, t-SNE, or UMAP on SNP HDF5 data",
    )
    _setup_snmf_subparser(
        subparsers,
        f"{header}\nipyrad2 snmf: run sNMF-style clustering on SNP HDF5 data",
    )
    _setup_dapc_subparser(
        subparsers,
        f"{header}\nipyrad2 dapc: run DAPC-style clustering on SNP HDF5 data",
    )
    _setup_admixture_subparser(
        subparsers,
        f"{header}\nipyrad2 admixture: run external ADMIXTURE on SNP HDF5 data",
    )
    _setup_popgen_subparser(
        subparsers,
        f"{header}\nipyrad2 popgen: compute genome-wide population-genetic statistics",
    )
    _setup_bpp_subparser(
        subparsers,
        f"{header}\nipyrad2 bpp: stage one BPP analysis from sequence HDF5 data",
    )
    _setup_reserved_tool_subparser(
        subparsers,
        name="baba",
        help_text="Reserved ABBA/BABA admixture metrics command.",
        header=f"{header}\nipyrad2 baba: reserved ABBA/BABA admixture metrics command",
    )
    _setup_reserved_tool_subparser(
        subparsers,
        name="treeslider",
        help_text="Reserved per-locus or per-window gene-tree command.",
        header=f"{header}\nipyrad2 treeslider: reserved per-locus or per-window gene-tree command",
    )


def _tool_name(args) -> str:
    """Return the active top-level export or analysis command name."""
    return getattr(args, "tool", args.subcommand)


def run_analysis_tool(args, _exit: bool = True) -> None:
    """Dispatch one top-level export or analysis command."""
    tool = _tool_name(args)

    if tool == "wex":
        logger.info("-------------------------------------------------------")
        logger.info("----- ipyrad2 wex: extract alignments from windows -----")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "lex":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 lex: extract delimited loci from HDF5 database ----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_locus_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            out_format=args.out_format,
            nloci=args.max_loci,
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
        if _exit:
            sys.exit(0)
        return

    if tool == "snpex":
        logger.info("------------------------------------------------------------")
        logger.info("---- ipyrad2 snpex: extract filtered SNP matrices ----")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "vcf2hdf5":
        logger.info("--------------------------------------------------")
        logger.info("---- ipyrad2 vcf2hdf5: convert VCF to HDF5 ----")
        logger.info("--------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_vcf_to_hdf5(
            data=args.data,
            name=args.name,
            outdir=args.out,
            ld_block_size=args.ld_block_size,
            force=args.force,
        )
        if _exit:
            sys.exit(0)
        return

    if tool == "pca":
        logger.info("-------------------------------------------------")
        logger.info("---- ipyrad2 pca: numerical PCA-family methods ----")
        logger.info("-------------------------------------------------")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "snmf":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 snmf: numerical clustering on SNPs ----")
        logger.info("-------------------------------------------------------")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "dapc":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 dapc: numerical clustering on SNPs ----")
        logger.info("-------------------------------------------------------")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "admixture":
        logger.info("-----------------------------------------------------------")
        logger.info("---- ipyrad2 admixture: external clustering on SNPs ----")
        logger.info("-----------------------------------------------------------")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "popgen":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 popgen: genome-wide population genetics ----")
        logger.info("-------------------------------------------------------")
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
        if _exit:
            sys.exit(0)
        return

    if tool == "bpp":
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 bpp: single-run BPP staging and execution ----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_bpp_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            tree=args.tree,
            imap=args.imap,
            minmap=args.minmap,
            max_loci=args.max_loci,
            min_length=args.min_length,
            msc_i=args.msc_i,
            msc_m=args.msc_m,
            speciestree=args.speciestree,
            speciesdelimitation=args.speciesdelimitation,
            thetaprior=args.thetaprior,
            tauprior=args.tauprior,
            speciesmodelprior=args.speciesmodelprior,
            phiprior=args.phiprior,
            wprior=args.wprior,
            alphaprior=args.alphaprior,
            locusrate=args.locusrate,
            clock=args.clock,
            burnin=args.burnin,
            samplefreq=args.samplefreq,
            nsample=args.nsample,
            threads=args.threads,
            seed=args.seed,
            write_only=args.write_only,
            force=args.force,
            log_level=args.log_level,
        )
        if _exit:
            sys.exit(0)
        return

    if tool in RESERVED_TOOL_NAMES:
        raise IPyradError(f"`ipyrad2 {tool}` is reserved but not implemented yet.")

    raise IPyradError(f"Unknown export or analysis command: {tool}")
