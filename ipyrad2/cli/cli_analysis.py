#!/usr/bin/env python

"""Top-level export and analysis command registration and dispatch."""

from __future__ import annotations

import argparse
import importlib
import sys
from typing import Optional, Sequence

from loguru import logger

from .command_log import format_logged_command
from .common import RAW_HELP_FORMATTER
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
    "baba",
    "bpp",
    "treeslider",
)

RESERVED_TOOL_NAMES = ()


_PARSER_SPECS = {
    "wex": (".cli_wex", "_setup_wex_subparser", "ipyrad2 wex: extract one alignment from selected genomic windows"),
    "lex": (".cli_lex", "_setup_lex_subparser", "ipyrad2 lex: extract delimited loci from HDF5 database"),
    "treeslider": (".cli_treeslider", "_setup_treeslider_subparser", "ipyrad2 treeslider: extract filtered windows and infer one tree per window"),
    "snpex": (".cli_snpex", "_setup_snpex_subparser", "ipyrad2 snpex: extract filtered SNP matrices from HDF5 database"),
    "vcf2hdf5": (".cli_vcf_to_hdf5", "_setup_vcf_to_hdf5_subparser", "ipyrad2 vcf2hdf5: convert VCF to SNP-capable HDF5 database"),
    "pca": (".cli_pca", "_setup_pca_subparser", "ipyrad2 pca: run PCA, t-SNE, or UMAP on SNP HDF5 data"),
    "snmf": (".cli_snmf", "_setup_snmf_subparser", "ipyrad2 snmf: run sNMF-style clustering on SNP HDF5 data"),
    "dapc": (".cli_dapc", "_setup_dapc_subparser", "ipyrad2 dapc: run DAPC-style clustering on SNP HDF5 data"),
    "admixture": (".cli_admixture", "_setup_admixture_subparser", "ipyrad2 admixture: run external ADMIXTURE on SNP HDF5 data"),
    "popgen": (".cli_popgen", "_setup_popgen_subparser", "ipyrad2 popgen: compute genome-wide population-genetic statistics"),
    "baba": (".cli_baba", "_setup_baba_subparser", "ipyrad2 baba: compute ABBA/BABA admixture metrics from SNP HDF5 data"),
    "bpp": (".cli_bpp", "_setup_bpp_subparser", "ipyrad2 bpp: stage one BPP analysis from sequence HDF5 data"),
}


_RUNTIME_RUNNERS = {
    "wex": ("..analysis.extracters.window_extracter", "run_window_extracter"),
    "lex": ("..analysis.extracters.locus_extracter", "run_locus_extracter"),
    "treeslider": ("..analysis.methods.treeslider", "run_treeslider_method"),
    "snpex": ("..analysis.extracters.snps_extracter", "run_snps_extracter"),
    "vcf2hdf5": ("..analysis.converters.vcf_to_hdf5", "run_vcf_to_hdf5"),
    "pca": ("..analysis.methods.pca", "run_pca_method"),
    "snmf": ("..analysis.methods.snmf", "run_snmf_method"),
    "dapc": ("..analysis.methods.dapc", "run_dapc_method"),
    "admixture": ("..analysis.methods.admixture", "run_admixture_method"),
    "popgen": ("..analysis.methods.popgen.runner", "run_popgen_method"),
    "baba": ("..analysis.methods.baba.runner", "run_baba_method"),
    "bpp": ("..analysis.methods.bpp", "run_bpp_method"),
}


def _load_runner(tool: str):
    """Load one runtime analysis runner lazily on first use."""
    module_name, attr_name = _RUNTIME_RUNNERS[tool]
    module = importlib.import_module(module_name, package=__package__)
    return getattr(module, attr_name)


def _setup_tool_subparser(
    subparsers: argparse._SubParsersAction,
    tool: str,
    header: str,
) -> None:
    """Load one parser module lazily and add its subparser."""
    module_name, setup_name, description = _PARSER_SPECS[tool]
    module = importlib.import_module(module_name, package=__package__)
    setup = getattr(module, setup_name)
    setup(subparsers, f"{header}\n{description}")


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
    selected_tools: Optional[Sequence[str]] = None,
) -> None:
    """Add top-level export and analysis subcommand parsers."""
    header = HEADER if header is None else header
    selected = (
        set(ANALYSIS_TOOL_NAMES + RESERVED_TOOL_NAMES)
        if selected_tools is None
        else set(selected_tools)
    )

    if "wex" in selected:
        _setup_tool_subparser(subparsers, "wex", header)
    if "lex" in selected:
        _setup_tool_subparser(subparsers, "lex", header)
    if "snpex" in selected:
        _setup_tool_subparser(subparsers, "snpex", header)
    if "vcf2hdf5" in selected:
        _setup_tool_subparser(subparsers, "vcf2hdf5", header)
    if "pca" in selected:
        _setup_tool_subparser(subparsers, "pca", header)
    if "snmf" in selected:
        _setup_tool_subparser(subparsers, "snmf", header)
    if "dapc" in selected:
        _setup_tool_subparser(subparsers, "dapc", header)
    if "admixture" in selected:
        _setup_tool_subparser(subparsers, "admixture", header)
    if "popgen" in selected:
        _setup_tool_subparser(subparsers, "popgen", header)
    if "baba" in selected:
        _setup_tool_subparser(subparsers, "baba", header)
    if "bpp" in selected:
        _setup_tool_subparser(subparsers, "bpp", header)
    if "treeslider" in selected:
        _setup_tool_subparser(subparsers, "treeslider", header)

def _tool_name(args) -> str:
    """Return the active top-level export or analysis command name."""
    return getattr(args, "tool", args.subcommand)


def run_analysis_tool(args, _exit: bool = True) -> None:
    """Dispatch one top-level export or analysis command."""
    tool = _tool_name(args)

    if tool == "wex":
        run_window_extracter = _load_runner(tool)
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
        run_locus_extracter = _load_runner(tool)
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

    if tool == "treeslider":
        run_treeslider_method = _load_runner(tool)
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 treeslider: infer one tree per filtered window ----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_treeslider_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            window_size=args.window_size,
            slide_size=args.slide_size,
            print_scaffold_table=args.print_scaffold_table,
            scaffolds=args.scaffolds,
            min_sample_coverage=args.min_sample_coverage,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            include_reference=args.include_reference,
            min_sample_alignment_length=args.min_sample_alignment_length,
            min_alignment_length=args.min_alignment_length,
            threads=args.threads,
            workers=args.workers,
            bs_trees=args.bs_trees,
            model=args.model,
            raxml_ng_binary=args.raxml_ng_binary,
            seed=args.seed,
            force=args.force,
            redo=args.redo,
            log_level=args.log_level,
        )
        if _exit:
            sys.exit(0)
        return

    if tool == "snpex":
        run_snps_extracter = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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
        run_vcf_to_hdf5 = _load_runner(tool)
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
        run_pca_method = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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
            colors=args.plot_colors,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
        )
        if _exit:
            sys.exit(0)
        return

    if tool == "snmf":
        run_snmf_method = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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
        run_dapc_method = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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
        run_admixture_method = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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
        run_popgen_method = _load_runner(tool)
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
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
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

    if tool == "baba":
        run_baba_method = _load_runner(tool)
        logged_command = format_logged_command(sys.argv[1:])
        logger.info("-------------------------------------------------------")
        logger.info("---- ipyrad2 baba: ABBA/BABA admixture statistics ----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
        run_baba_method(
            data=args.data,
            name=args.name,
            outdir=args.out,
            tests=args.tests,
            tree=args.tree,
            imap=args.imap,
            minmap=args.minmap,
            min_sample_coverage=args.min_sample_coverage,
            min_genotype_depth=args.min_genotype_depth,
            min_site_qual=args.min_site_qual,
            exclude=args.exclude,
            include_reference=args.include_reference,
            resampling=args.resampling,
            bootstrap_replicates=args.bootstrap_replicates,
            jackknife_block_bp=args.jackknife_block_bp,
            jackknife_block_loci=args.jackknife_block_loci,
            seed=args.seed,
            f_branch=args.f_branch,
            f_branch_p_threshold=args.f_branch_p_threshold,
            write_block_table=args.write_block_table,
            clustering_stats=args.clustering_stats,
            cores=args.cores,
            force=args.force,
            log_level=args.log_level,
            logged_command=logged_command,
        )
        if _exit:
            sys.exit(0)
        return

    if tool == "bpp":
        run_bpp_method = _load_runner(tool)
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
