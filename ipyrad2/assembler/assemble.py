#!/usr/bin/env python

"""...
"""

from typing import List, Tuple
import sys
from pathlib import Path
from loguru import logger
from ..utils.parallel import run_with_pool
from ..utils.exceptions import IPyradError
from .beds import (
    get_name_from_bam,
    get_reference_sort_order,
    get_fragment_beds,
    get_fragment_coverage_beds,
    get_fragment_merged_coverage_beds,
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
)
from .loci import (
    write_sam_faidx,
    get_reference_in_loci_beds,
    get_consensus,
    get_sample_masked_beds,
    build_locus_fasta_database,
    write_loci_and_stats_files,
)
from .write_seqs_hdf5 import write_seqs_hdf5


def run_assembler(
    rad_bams: List[Path],
    wgs_bams: List[Path] | None,
    reference: Path,
    outdir: Path,
    name: str,
    min_gq: int,
    min_qual: int,
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
    name_parse: Tuple[str, int] | None,
    workers: int,
    threads: int,
    force: bool,
    ):
    # check reference and outdir paths
    reference = reference.expanduser().absolute()
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(exist_ok=True)

    # check outdir for existing and raise or remove
    if (outdir / f"{name}.loci.txt").exists():
        if not force:
            raise IPyradError(f"outfiles with prefix {name} already exist in {outdir}. Use --force to overwrite.")
        else:
            # collect relevant files and rm
            pass
            # for n in ["vcfs", "beds", "consensus_seqs"]:
            #     odir = outdir / n
            #     if odir.exists():
            #         shutil.rmtree(odir)

    # check bam paths and get names dicts as {name: Path, ...}
    bam_dict = {}
    if rad_bams:
        for bam_file in rad_bams:
            bam_dict[get_name_from_bam(bam_file)] = bam_file.expanduser().absolute()
    assert bam_dict, "..."
    logger.info(f"loaded {len(bam_dict)} RAD samples") # : {sorted(bam_dict)}")

    wgs_dict = {}
    if wgs_bams:
        for bam_file in wgs_bams:
            wgs_dict[get_name_from_bam(bam_file)] = bam_file.expanduser().absolute()
    if wgs_dict:
        logger.info(f"loaded {len(wgs_dict)} WGS samples") # : {sorted(wgs_dict)}")

    # all samples
    all_dict = wgs_dict | bam_dict

    # ---------------------------------------------
    logger.info(f"running up to {workers} parallel jobs each using up to {threads} threads")
    logger.debug("fetching reference scaffold order")
    get_reference_sort_order(reference, outdir)

    logger.info("delimiting sample coverage beds")
    jobs = {}
    for sname, bam_file in bam_dict.items():
        kwargs = dict(sname=sname, bam_file=bam_file, outdir=outdir)
        jobs[sname] = kwargs
    results = run_with_pool(get_fragment_beds, jobs, workers)

    jobs = {}
    for sname, bam_file in bam_dict.items():
        kwargs = dict(sname=sname, reference=reference, outdir=outdir)
        jobs[sname] = kwargs
    results = run_with_pool(get_fragment_coverage_beds, jobs, workers)

    jobs = {}
    for sname, bam_file in bam_dict.items():
        kwargs = dict(sname=sname, outdir=outdir)
        jobs[sname] = kwargs
    results = run_with_pool(get_fragment_merged_coverage_beds, jobs, workers)

    logger.info("delimiting shared coverage beds (loci)")
    get_across_sample_loci_bed(
        list(bam_dict),
        min_locus_sample_coverage,
        min_locus_merge_distance,
        min_locus_length,
        outdir,
    )

    # Maybe not necessary, we measure coverage on the filtered loci later.
    # logger.info("measuring sample locus coverage stats")
    # jobs = {}
    # for sname, bam_file in all_dict.items():
    #     jobs[sname] = dict(bam_file=bam_file, outdir=outdir)
    # results = run_with_pool(get_sample_coverage_stats_in_loci_bed, jobs, workers)
    # locus_bed_coverage_stats = results

    # ------------------------------------------------------------------
    # ---- VARIANT CALLING ---------------------------------------------
    # ------------------------------------------------------------------
    logger.info("calling variants in locus beds")
    nchunks = max(4, int(workers / threads))
    locus_chunks = get_chunked_loci_beds(outdir, nchunks)
    jobs = {}
    for chunk in locus_chunks:
        kwargs = dict(outdir=outdir, reference=reference, bam_files=list(all_dict.values()), locus_chunk=chunk, threads=max(4, threads))
        jobs[chunk] = kwargs
    results = run_with_pool(get_group_called_variants_in_vcf_chunks, jobs, workers)
    get_concat_chunk_vcfs(outdir, threads)

    logger.info("filtering variants")
    get_filtered_vcf(outdir, min_sample_depth, min_gq, min_qual, max(4, threads))

    logger.info("resolving indels and snps")
    get_vcf_with_indels_resolved(outdir, reference, max(4, threads))

    # optional: maybe wait til after locus filtering...
    stats = get_locus_and_snp_stats_in_loci_bed(outdir, max(4, threads))
    logger.warning(stats)

    # ------------------------------------------------------------------
    # ---- CONSENSUS CALLING -------------------------------------------
    # ------------------------------------------------------------------
    logger.info("extracting reference sequence in locus beds")
    write_sam_faidx(outdir)
    get_reference_in_loci_beds(outdir, reference)

    logger.info("building coverage masks")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(sname=sname, bam_file=bam_file, min_sample_depth=min_sample_depth, outdir=outdir)
        jobs[sname] = kwargs
    results = run_with_pool(get_sample_masked_beds, jobs, workers)

    logger.info("extracting consensus sequences")
    jobs = {}
    for sname, bam_file in all_dict.items():
        kwargs = dict(sname=sname, reference=reference, outdir=outdir, keep_insertions=False)
        jobs[sname] = kwargs
    results = run_with_pool(get_consensus, jobs, workers)

    # ------------------------------------------------------------------
    # ---- LOCUS BUILDING ----------------------------------------------
    # ------------------------------------------------------------------
    logger.info("assembling loci")
    build_locus_fasta_database(
        name,
        list(all_dict),
        reference,
        outdir,
        exclude_reference,
        masks,
    )

    logger.info("filtering and writing loci and stats")
    stats = write_loci_and_stats_files(
        list(all_dict),
        name,
        outdir,
        exclude_reference,
        min_locus_sample_coverage,
        min_locus_trim_sample_coverage,
        min_locus_length,
        max_locus_hetero_frequency,
        max_locus_variant_frequency,
    )

    # write seqs HDF5, snps HDF5, and final VCFs
    logger.info("writing hdf5 database files")
    write_seqs_hdf5(
        name=name,
        outdir=outdir,
        snames=list(all_dict),
        reference=reference,
        exclude_reference=exclude_reference,
        nloci=stats["nloci"],
        nsites=stats["nsites"],
    )

    # write snps HDF5


    #
    # logger.info(locus_bed_coverage_stats)



if __name__ == "__main__":
    pass
