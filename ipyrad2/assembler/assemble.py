#!/usr/bin/env python

"""...
"""

from typing import List
import shutil
from pathlib import Path
from loguru import logger
import numpy as np
import pandas as pd
from .beds import (
    get_name_from_bam,
    get_reference_sort_order,
    get_fragment_merged_coverage_beds,
    get_coverage_bed_graphs,
    get_across_sample_loci_bed,
    get_sample_coverage_stats_in_loci_bed,
)
from .variants import (
    get_chunked_loci_beds,
    get_group_called_variants_in_vcf_chunks,
    get_concat_chunk_vcfs,
    get_filtered_vcf,
    get_vcf_with_indels_resolved,
    get_locus_and_snp_stats_in_loci_bed,
    write_vcf,
)
from .loci import (
    write_sam_faidx,
    get_reference_in_loci_beds,
    get_consensus,
    get_sample_masked_beds,
    build_locus_fasta_database,
    write_loci_and_stats_files,
)
from .write_seqs import write_seqs_hdf5
from .write_snps import write_snps_hdf5
from ..utils.parallel import run_with_pool
from ..utils.exceptions import IPyradError


def run_assembler(
    rad_bams: List[Path],
    wgs_bams: List[Path] | None,
    reference: Path,
    outdir: Path,
    name: str,
    loci_bed: Path | None,
    min_map_q: int,
    min_site_q: int,
    min_geno_q: int,
    min_base_q: int,
    min_sample_depth: int,                # sample must have depth cov or site is masked.
    min_locus_sample_coverage: int,       # locus must have data for N samples (used in locus delim)
    min_locus_trim_sample_coverage: int,  # trim r/l to region with at least N samples data (default 4)
    min_locus_length: int,
    min_locus_merge_distance: int,        # merge loci within this distance
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    populations: Path,
    masks: List[str] | None,
    exclude_reference: bool,
    cores: int,
    threads: int,
    force: bool,
    log_level: str,
    ):
    # expand paths
    loci_bed = loci_bed.expanduser().absolute() if loci_bed else None
    reference = reference.expanduser().absolute()
    outdir = outdir.expanduser().absolute()
    tmpdir = outdir / "tmpdir"

    # run this many multithreaded jobs concurrently
    workers = max(1, cores // threads)

    # check outdir for existing and raise or remove
    if (outdir / f"{name}.loci.txt").exists():
        if not force:
            raise IPyradError(f"outfiles with prefix {name} already exist in {outdir}. Use --force to overwrite.")
        else:
            # collect relevant files and rm
            logger.debug(f"removing previous ipyrad assemble files from {outdir}")
            if tmpdir.exists():
                shutil.rmtree(tmpdir)
            rfiles = [
                outdir / f"{name}.loci.txt",
                outdir / f"{name}.seqs.hdf5",
                outdir / f"{name}.snps.hdf5",
                outdir / f"{name}.stats_loci.tsv",
                outdir / f"{name}.stats_samples.tsv",
                outdir / f"{name}.stats_coverage.tsv",
            ]
            for r in rfiles:
                if r.exists():
                    r.unlink()

    # ensure directory structure
    outdir.mkdir(exist_ok=True, parents=True)
    tmpdir.mkdir(exist_ok=True)
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(exist_ok=True)
    vcf_dir = tmpdir / "vcfs"
    vcf_dir.mkdir(exist_ok=True)
    consensus_dir = tmpdir / "consensus_seqs"
    consensus_dir.mkdir(exist_ok=True)

    # raises exception if no RAD bams found. Fills {name: bam, ...}
    bam_dict = {}
    if rad_bams:
        for bam_file in rad_bams:
            sname = get_name_from_bam(bam_file)
            if sname in bam_dict:
                raise IPyradError(f"Multiple input files of sample name {sname}")
            bam_dict[sname] = bam_file.expanduser().absolute()
    if not bam_dict:
        raise IPyradError("No RAD bam files found. These are required.")
    logger.info(f"loaded {len(bam_dict)} RAD samples")

    # not required. Fills {name: bam, ...}
    wgs_dict = {}
    if wgs_bams:
        for bam_file in wgs_bams:
            sname = get_name_from_bam(bam_file)
            if sname in wgs_dict:
                raise IPyradError(f"Multiple input files of sample name {sname}")
            wgs_dict[sname] = bam_file.expanduser().absolute()
    if wgs_dict:
        logger.info(f"loaded {len(wgs_dict)} WGS samples")

    # all samples
    all_dict = wgs_dict | bam_dict
    snames = sorted(all_dict)
    all_dict = {i: all_dict[i] for i in snames}

    # ---------------------------------------------
    logger.info(f"using up to {cores} cores (up to {workers} multi-threaded jobs using {threads} threads)")
    logger.debug("fetching reference scaffold order")
    get_reference_sort_order(reference, tmpdir)

    # ------------------------------------------------------------------
    # ---- LOCUS DELIMITING --------------------------------------------
    # ------------------------------------------------------------------
    if loci_bed is not None:
        # copy input loci file to tmpdir/beds/loci.bed
        loci_bed = shutil.copy2(loci_bed, tmpdir / "beds" / "loci.bed")
        loci_bed = tmpdir / "beds" / "loci.bed"
    else:
        logger.info("delimiting sample coverage beds")
        jobs = {}
        for sname, bam_file in bam_dict.items():
            kwargs = dict(sname=sname, bam_file=bam_file, reference=reference, min_map_q=10, threads=threads, tmpdir=tmpdir)
            jobs[sname] = (get_coverage_bed_graphs, kwargs)
        run_with_pool(jobs, log_level, workers)            # multithreaded

        jobs = {}
        for sname, bam_file in bam_dict.items():
            kwargs = dict(sname=sname, tmpdir=tmpdir)
            jobs[sname] = (get_fragment_merged_coverage_beds, kwargs)
        run_with_pool(jobs, log_level, cores)              # single-threaded

        logger.info("delimiting shared coverage beds (loci)")
        args = (list(bam_dict), min_locus_sample_coverage, min_locus_merge_distance, min_locus_length, tmpdir)
        get_across_sample_loci_bed(*args)

    # if tmpdir/beds/loci.bed is empty, bail
    if (tmpdir / "beds" / "loci.bed").stat().st_size == 0:
        logger.error("No loci passed filtering")
        raise SystemExit(1)

    # --------------------------------------------------------------------
    # ---- STATS FOR DEPTH FILTERING (OPTIONAL?) -------------------------
    # --------------------------------------------------------------------
    # get the coverage stats per locus before filtering. We do this here
    # so we can apply depth-outlier filters, and report filter stats later.
    # This takes a while to run on large WGS samples, and doesn't seem to
    # often lead to many filtered loci. Maybe make it optional?
    logger.info("measuring per locus stats")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(bam_file=bam_file, loci_bed=tmpdir / "beds" / "loci.bed", min_map_q=min_map_q, ref_info=tmpdir / "REF_info.txt")
        jobs[sname] = (get_sample_coverage_stats_in_loci_bed, kwargs)
    results = run_with_pool(jobs, log_level, workers)

    # store locus stats before filtering.
    stats_before = pd.DataFrame(data={i: results[i][0] for i in snames}).T

    # store read depth outlier mask
    covs = [i[1] for i in results.values()]
    read_depth_zscores = np.vstack(covs).mean(axis=0)
    read_depth_mask = read_depth_zscores > 5.0
    # logger.warning(read_depth_mask)

    # ------------------------------------------------------------------
    # ---- VARIANT CALLING ---------------------------------------------
    # ------------------------------------------------------------------
    logger.info("calling variants in locus beds")
    # Each variant calling worker can effectively use ~2 threads
    vworkers = int(cores // 2)
    # nchunks = max(10, vworkers)
    locus_chunks = get_chunked_loci_beds(tmpdir, vworkers)
    jobs = {}
    for chunk in locus_chunks:
        kwargs = dict(
            tmpdir=tmpdir,
            reference=reference,
            bam_files=list(all_dict.values()),
            min_base_q=min_base_q,
            min_map_q=min_map_q,
            locus_chunk=chunk,
            threads=min(2, threads),
        )
        jobs[chunk] = (get_group_called_variants_in_vcf_chunks, kwargs)
    # each job does not use much RAM, so run many at 1-3 threads per job
    run_with_pool(jobs, log_level, vworkers)      # <=3 threads per job

    # write concatenated loci chunks
    get_concat_chunk_vcfs(tmpdir, threads)

    logger.info("filtering variants")
    # TODO: consider other quality filters here?
    get_filtered_vcf(tmpdir, min_sample_depth, min_geno_q, min_site_q, threads)

    logger.info("resolving indels and snps")
    get_vcf_with_indels_resolved(tmpdir, reference, threads)

    # ------------------------------------------------------------------
    # ---- CONSENSUS CALLING -------------------------------------------
    # ------------------------------------------------------------------
    logger.info("extracting reference sequence in locus beds")
    write_sam_faidx(tmpdir)
    get_reference_in_loci_beds(tmpdir, reference)

    logger.info("building coverage masks")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(sname=sname, bam_file=bam_file, min_sample_depth=min_sample_depth, tmpdir=tmpdir)
        jobs[sname] = (get_sample_masked_beds, kwargs)
    run_with_pool(jobs, log_level, cores)

    logger.info("extracting consensus sequences")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(sname=sname, reference=reference, tmpdir=tmpdir, keep_insertions=False)
        jobs[sname] = (get_consensus, kwargs)
    run_with_pool(jobs, log_level, workers)

    # ------------------------------------------------------------------
    # ---- LOCUS BUILDING ----------------------------------------------
    # ------------------------------------------------------------------
    logger.info("assembling loci")
    args = (name, snames, reference, tmpdir, exclude_reference, masks)
    build_locus_fasta_database(*args)

    # ------------------------------------------------------------------
    # ---- DATABASE WRITING --------------------------------------------
    # ------------------------------------------------------------------
    logger.info("writing outfiles (.loci, .hdf5, .bed, .stats_*)")
    jobs = {
        "loci": (write_loci_and_stats_files, dict(
            snames=snames,
            name=name,
            outdir=outdir,
            exclude_reference=exclude_reference,
            min_locus_sample_coverage=min_locus_sample_coverage,
            min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
            min_locus_length=min_locus_length,
            max_locus_hetero_frequency=max_locus_hetero_frequency,
            max_locus_variant_frequency=max_locus_variant_frequency,
            read_depth_mask=read_depth_mask,
        )),
        "seqs": (write_seqs_hdf5, dict(
            name=name,
            outdir=outdir,
            snames=snames,
            reference=reference,
            exclude_reference=exclude_reference,
            min_locus_sample_coverage=min_locus_sample_coverage,
            min_locus_trim_sample_coverage=min_locus_sample_coverage,
            min_locus_length=min_locus_length,
            max_locus_hetero_frequency=max_locus_hetero_frequency,
            max_locus_variant_frequency=max_locus_variant_frequency,
            read_depth_mask=read_depth_mask,
        ))
    }
    run_with_pool(jobs, log_level, workers)

    # get the final vcf file
    logger.info("writing variants file (.vcf.gz)")
    write_vcf(name, outdir, threads)

    # add snps dataset to the database file
    logger.info("writing snps database (.hdf5)")
    write_snps_hdf5(name, outdir, list(all_dict), reference)

    # ------------------------------------------------------------------
    # ---- STATS CALC/WRITING ------------------------------------------
    # ------------------------------------------------------------------
    # calculate and write depths per (final) loci for each sample.
    logger.info("measuring per locus stats")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(bam_file=bam_file, loci_bed=outdir / f"{name}.bed", min_map_q=min_map_q, ref_info=tmpdir / "REF_info.txt")
        jobs[sname] = (get_sample_coverage_stats_in_loci_bed, kwargs)
    results = run_with_pool(jobs, log_level, workers)
    df = pd.DataFrame(data={i: results[i][0] for i in snames}).T
    df["nloci_before_filtering"] = stats_before["nloci"]
    df["mean_depth_per_locus_with_nonzero_mapping_before_filtering"] = stats_before["mean_depth_per_locus_with_nonzero_mapping"]
    df["median_depth_per_locus_with_nonzero_mapping_before_filtering"] = stats_before["median_depth_per_locus_with_nonzero_mapping"]
    df["std_depth_per_locus_with_nonzero_mapping_before_filtering"] = stats_before["std_depth_per_locus_with_nonzero_mapping"]
    df.to_string(
        outdir / f"{name}.stats_depths.txt",
        formatters={
            "nloci": lambda x: f"{int(x)}",
            "nloci_before_filtering": lambda x: f"{int(x)}",
            "mean_depth_per_locus_total": lambda x: f"{x:.3f}",
            "median_depth_per_locus_total": lambda x: f"{int(x)}",
            "std_depth_per_locus_total": lambda x: f"{x:.3f}",
            "mean_depth_per_locus_with_nonzero_mapping": lambda x: f"{x:.3f}",
            "median_depth_per_locus_with_nonzero_mapping": lambda x: f"{int(x)}",
            "std_depth_per_locus_with_nonzero_mapping": lambda x: f"{x:.3f}",
            "mean_depth_per_locus_with_nonzero_mapping_before_filtering": lambda x: f"{x:.3f}",
            "median_depth_per_locus_with_nonzero_mapping_before_filtering": lambda x: f"{int(x)}",
            "std_depth_per_locus_with_nonzero_mapping_before_filtering": lambda x: f"{x:.3f}",
        },
    )

    # calculate and report number of variants in final loci
    # covs = [i[1] for i in results.values()]
    stats = get_locus_and_snp_stats_in_loci_bed(tmpdir, threads)
    logger.debug(f"{stats}")


if __name__ == "__main__":
    pass
