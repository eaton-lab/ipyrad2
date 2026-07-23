#!/usr/bin/env python

"""Consensus, locus-database, and final report helpers for assemble."""

from __future__ import annotations

import json
import gzip
import csv
import re
import shutil
import sys
from pathlib import Path
from collections import Counter
import numpy as np
from loguru import logger
from .sort_utils import assemble_sort_with_args
from ..utils.seqs import comp
from ..utils.jit_funcs import snp_count_numba, max_heteros_count_numba
from ..utils.parallel import run_pipeline, run_with_pool_iter

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_BED = str(BIN / "bedtools")
HETERO_CODES = np.array(list(b"RSKYWM"), dtype=np.uint8)
REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"


def get_lowdepth_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the intermediate low-depth-only mask path for one sample."""
    return tmpdir / "beds" / f"{sname}.lowdepth.mask.bed"


def get_goodcov_bed_path(sname: str, tmpdir: Path) -> Path:
    """Return the per-sample BED of positions meeting the depth threshold."""
    return tmpdir / "beds" / f"{sname}.goodcov.bed"


def get_final_good_bed_path(sname: str, tmpdir: Path) -> Path:
    """Return the per-sample retained BED after shared paralog filtering."""
    return tmpdir / "beds" / f"{sname}.final.good.bed"


def get_paralog_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the sample-specific paralog exclusion BED path."""
    return tmpdir / "beds" / f"{sname}.paralog.mask.bed"


def get_indel_overlap_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the sample-specific overlapping-indel-cluster mask path."""
    return tmpdir / "beds" / f"{sname}.indel_overlap.mask.bed"


def get_sample_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the merged mask path consumed by consensus calling."""
    return tmpdir / "beds" / f"{sname}.mask.bed"


def get_consensus_sample_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the final-output BED path for sample-row masking in consensus outputs."""
    return tmpdir / "beds" / f"{sname}.consensus_sample.mask.bed"


def get_final_vcf_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the final merged sample-mask BED path used on the final VCF."""
    return tmpdir / "beds" / f"{sname}.final.vcf.mask.bed"


def get_retained_loci_manifest_path(name: str, tmpdir: Path) -> Path:
    """Return the manifest path describing retained/masked final loci."""
    return tmpdir / f"{name}.retained_loci.tsv"


def write_sam_faidx(tmpdir: Path) -> Path:
    """Convert the shared loci BED into samtools faidx-style region strings."""
    loci_bed = tmpdir / "beds" / "loci.bed"
    fai_path = tmpdir / "loci.faidx.txt"
    awk_prog = 'BEGIN{OFS=""}{print $1,":",$2+1,"-",$3}'
    cmd = ["awk", awk_prog, str(loci_bed)]
    run_pipeline([cmd], fai_path)
    return fai_path


def get_reference_in_loci_beds(tmpdir: Path, reference: Path) -> Path:
    """Write the locus-sliced reference FASTA reused by every consensus job."""
    loci = tmpdir / "loci.faidx.txt"
    consensus_dir = tmpdir / "consensus_seqs"
    out_fasta = consensus_dir / "assembly_reference_sequence.consensus.fa"

    # run pipeline
    cmd = [BIN_SAM, "faidx", str(reference), "-r", str(loci)]
    run_pipeline([cmd], out_fasta)
    return out_fasta


def get_consensus(
    *,
    sname: str,
    reference_fasta: Path,
    resolved_vcf: Path,
    sample_mask_bed: Path,
    out_fasta: Path,
    keep_insertions: bool,
) -> Path:
    """Write consensus sequences for one sample.

    Create FASTA for `sample_name` over the already sliced shared loci
    reference, applying variants from `vcf_gz` and masking filtered regions to N.
    """
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd1 = [
        BIN_BCF, "consensus",
        "-f", str(reference_fasta),
        "-s", f"{sname}",         # sample to apply
        "-M", "N",                # write N for missing genotypes
        "--mask", str(sample_mask_bed),  # mask zero/low-coverage intervals to N
        "--mask-with", "N",
        "--mark-del", "-",
        "--mark-ins", "lc" if keep_insertions else "+",
        "--regions-overlap", "1", # apply variants overlapping slice edges
        str(resolved_vcf),
    ]
    cmd2 = ["tr", "-d", "'+"]
    run_pipeline([cmd1, cmd2], out_fasta)

    # warn if there is no data for a sample.
    if not out_fasta.stat().st_size:
        logger.warning(f"sample {sname} has no data passed filtering and should be dropped.")
    return out_fasta


def _subtract_sorted_beds(a_bed: Path, b_bed: Path, ref_info: Path, out_bed: Path) -> Path:
    """Subtract one reference-sorted BED from another using the canonical genome order."""
    cmd1 = [
        BIN_BED,
        "subtract",
        "-a",
        str(a_bed),
        "-b",
        str(b_bed),
        "-sorted",
        "-g",
        str(ref_info),
    ]
    cmd2 = ["cut", "-f", "1-3"]
    run_pipeline([cmd1, cmd2], out_bed)
    return out_bed


def make_lowdepth_mask(
    *,
    loci_bed: Path,
    sample_bedgraph: Path,
    ref_info: Path,
    good_bed: Path,
    out_bed: Path,
    sort_tmpdir: Path,
    min_sample_depth: int,
) -> Path:
    """Build a per-bp mask of positions inside `loci_bed` where bedGraph depth < min_depth.

    Output mask contains only the A (loci) columns and is split into minimal sub-intervals
    where coverage is below threshold (including 0-coverage gaps).
    """
    # 1) Threshold bedGraph: keep depth >= min_depth, drop depth column for set ops
    cmd1 = [
        "awk",
        f'BEGIN{{OFS="\\t"}} $4>={min_sample_depth} {{print $1,$2,$3}}',
        str(sample_bedgraph),
    ]
    cmd2 = assemble_sort_with_args(["-k1,1", "-k2,2n", "-T", str(sort_tmpdir)])
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    cmd4 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], good_bed)

    return _subtract_sorted_beds(loci_bed, good_bed, ref_info, out_bed)


def make_paralog_mask(
    *,
    loci_bed: Path,
    sample_good_bed: Path,
    ref_info: Path,
    out_bed: Path,
) -> Path:
    """Write the shared loci segments excluded only for this sample by paralog filtering."""
    # Samples that had no passing per-sample BED after paralog scoring get an
    # empty mask here so the downstream merge step still has a stable filepath.
    if not sample_good_bed.exists():
        out_bed.write_text("", encoding="utf-8")
        return out_bed

    return _subtract_sorted_beds(loci_bed, sample_good_bed, ref_info, out_bed)


def merge_sample_mask_beds(
    *,
    lowdepth_bed: Path,
    paralog_bed: Path,
    indel_overlap_bed: Path,
    ref_info: Path,
    out_bed: Path,
    sort_tmpdir: Path,
) -> Path:
    """Merge all per-sample mask sources into the final consensus mask BED."""
    existing = [
        path
        for path in (lowdepth_bed, paralog_bed, indel_overlap_bed)
        if path.exists() and path.stat().st_size
    ]
    if not existing:
        out_bed.write_text("", encoding="utf-8")
        return out_bed
    if len(existing) == 1:
        shutil.copy2(existing[0], out_bed)
        return out_bed

    # Merge every active interval source into one sorted mask so consensus
    # calling sees a single BED regardless of why a site was filtered.
    cmd1 = ["cat"] + [str(path) for path in existing]
    cmd2 = assemble_sort_with_args(["-k1,1", "-k2,2n", "-T", str(sort_tmpdir)])
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    cmd4 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], out_bed)
    return out_bed


def merge_final_vcf_mask_beds(
    *,
    lowdepth_bed: Path,
    indel_overlap_bed: Path,
    consensus_sample_bed: Path,
    ref_info: Path,
    out_bed: Path,
    sort_tmpdir: Path,
) -> Path:
    """Merge the final-output-only sample masks applied to the final VCF.

    This intentionally excludes the sample-specific paralog mask because
    post-paralog calling BAMs already remove that evidence before joint calling.
    """
    existing = [
        path
        for path in (lowdepth_bed, indel_overlap_bed, consensus_sample_bed)
        if path.exists() and path.stat().st_size
    ]
    if not existing:
        out_bed.write_text("", encoding="utf-8")
        return out_bed
    if len(existing) == 1:
        shutil.copy2(existing[0], out_bed)
        return out_bed

    cmd1 = ["cat"] + [str(path) for path in existing]
    cmd2 = assemble_sort_with_args(["-k1,1", "-k2,2n", "-T", str(sort_tmpdir)])
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    cmd4 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], out_bed)
    return out_bed


def iter_fasta(fasta: Path):
    """Stream a plain-text multi-FASTA and yield `(header, sequence)` tuples."""
    # Open if a path-like was given
    fh = open(fasta, "rt", encoding="utf-8")
    header = None
    parts: list[str] = []
    try:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(parts)
                    yield header, seq.upper()
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)

        # flush last record
        if header is not None:
            seq = "".join(parts)
            yield header, seq.upper()
    finally:
        fh.close()


def iter_build_loci(fastas: list[Path]):
    """Yield one assembled locus at a time from the ordered consensus FASTAs."""
    # do not re-sort fastas here, use the input order.
    iterators = [iter_fasta(i) for i in fastas]
    names = [i.name.rsplit(".consensus.fa")[0] for i in fastas]

    while 1:
        try:
            locus = []
            for fit in iterators:
                header, seq = next(fit)
                locus.append(seq)
            yield header, names, locus
        except StopIteration:
            break


def build_locus_fasta_database(
    *,
    consensus_fastas: list[Path],
    database_fasta: Path,
    restriction_mask_bed: Path,
    masks: list[str] | None,
) -> tuple[Path, Path]:
    """Build the shared locus FASTA database and optional restriction-site mask BED."""
    fastas = list(consensus_fastas)

    # get names
    snames = [i.name.rsplit(".consensus.fa")[0] for i in fastas]

    # restriction site sequences to be masked
    re_masks = []
    if masks:
        for mask in masks:
            re_masks.append(re.compile(mask))
            re_masks.append(re.compile(comp(mask)[::-1]))

    database_fasta.parent.mkdir(parents=True, exist_ok=True)
    restriction_mask_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(database_fasta, "w") as out_fa, open(restriction_mask_bed, "w") as out_bed:

        # iterate over loci pulled from fasta files
        lit = iter_build_loci(fastas)
        for header, names, locus in lit:

            # filter cut-sites from locus
            hits = set()
            if masks:
                for seq in locus:
                    for search in re_masks:
                        for hit in search.finditer(seq):
                            hits.add((hit.start(), hit.end()))

                # store masks to bed
                for h in sorted(hits):
                    scaff, pos = header.split(":", 1)
                    start = int(pos.split("-")[0])
                    out_bed.write(f"{scaff}\t{start + h[0]}\t{start + h[1]}\n")

            # build fasta
            loc = []
            for n, seq in zip(snames, locus):
                if len(seq) > seq.count("N"):
                    # mask RE sites
                    if hits:
                        seq = list(seq)
                        for h in hits:
                            seq[h[0]:h[1]] = "N" * (h[1] - h[0])
                        seq = "".join(seq)
                    # store locus
                    loc.append(f">{header} {n}\n{seq}")

            # write locus
            out_fa.write("\n".join(loc) + "\n\n")
    return database_fasta, restriction_mask_bed


def iter_parse_loci(database_fasta: Path):
    """Yield `(header, {sample_name: sequence})` records from `database.fa`."""
    ii = iter_fasta(database_fasta)
    last_scaff_pos = None
    while 1:
        try:
            locus = {}
            for fit in ii:
                header, seq = fit
                scaff_pos, sname = header.rsplit(" ", 1)
                if scaff_pos != last_scaff_pos:
                    if last_scaff_pos:
                        yield last_scaff_pos, locus
                    locus = {}
                    last_scaff_pos = scaff_pos
                locus[sname] = seq

            # flush last record
            if locus:
                yield last_scaff_pos, locus
                break

        except StopIteration:
            break


def iter_locus_batches(database_fasta: Path, batch_size: int = 128):
    """Yield fixed-size batches of parsed loci while preserving input order."""
    batch: list[tuple[str, dict[str, str]]] = []
    batch_idx = 0
    for item in iter_parse_loci(database_fasta):
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch_idx, batch
            batch_idx += 1
            batch = []
    if batch:
        yield batch_idx, batch


def _empirical_sample_row_mask(
    tnames: list[str],
    *,
    refname: str = REFERENCE_SAMPLE_NAME,
) -> np.ndarray:
    """Return a boolean row mask that excludes the synthetic reference row."""
    return np.array([sname != refname for sname in tnames], dtype=bool)


def _trim_locus_matrix(
    seqs: np.ndarray,
    *,
    min_locus_trim_sample_coverage: int,
    empirical_row_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Trim one locus matrix to the region with sufficient sample coverage."""
    if empirical_row_mask is None:
        empirical_seqs = seqs
    else:
        empirical_seqs = seqs[empirical_row_mask]
    site_sample_covs = np.sum((empirical_seqs != 78) & (empirical_seqs != 45), axis=0)
    cov_sufficient = np.where(site_sample_covs >= min_locus_trim_sample_coverage)[0]
    try:
        trim_left = int(cov_sufficient[0])
    except IndexError:
        trim_left = 0
    try:
        trim_right = seqs.shape[1] - int(cov_sufficient[-1]) - 1
    except IndexError:
        trim_right = 0
    trimmed = seqs[:, trim_left:seqs.shape[1] - trim_right]
    trimmed_cov = site_sample_covs[trim_left:seqs.shape[1] - trim_right]
    return trimmed, trimmed_cov, trim_left, trim_right


def _sample_rows_with_data(tseqs: np.ndarray) -> np.ndarray:
    """Return a boolean row mask for samples with at least one non-missing base."""
    if tseqs.size == 0:
        return np.zeros(tseqs.shape[0], dtype=bool)
    return np.any((tseqs != 78) & (tseqs != 45), axis=1)


def _mask_high_hetero_samples(
    tnames: list[str],
    tseqs: np.ndarray,
    *,
    max_sample_hetero_frequency: float,
) -> tuple[np.ndarray, list[str], dict[str, float]]:
    """Mask samples whose observed bases exceed the per-locus heterozygosity threshold."""
    if tseqs.size == 0:
        return tseqs.copy(), [], {}

    masked = tseqs.copy()
    masked_samples: list[str] = []
    sample_props: dict[str, float] = {}
    for row_idx, sname in enumerate(tnames):
        if sname == REFERENCE_SAMPLE_NAME:
            continue
        row = masked[row_idx]
        observed = (row != 78) & (row != 45)
        explicit_hetero = np.isin(row, HETERO_CODES)
        numer = int(explicit_hetero.sum())
        denom = int(observed.sum())
        sample_props[sname] = float(numer / denom) if denom else 0.0
        if denom and sample_props[sname] > max_sample_hetero_frequency:
            masked[row_idx, :] = np.uint8(ord("N"))
            masked_samples.append(sname)
    return masked, masked_samples, sample_props


def _mask_low_observed_samples(
    tnames: list[str],
    tseqs: np.ndarray,
    *,
    min_sample_observed_fraction: float,
) -> tuple[np.ndarray, list[str], dict[str, float]]:
    """Mask samples with too little observed sequence across the trimmed locus."""
    if tseqs.size == 0:
        return tseqs.copy(), [], {}

    masked = tseqs.copy()
    masked_samples: list[str] = []
    sample_props: dict[str, float] = {}
    locus_length = int(masked.shape[1])
    for row_idx, sname in enumerate(tnames):
        if sname == REFERENCE_SAMPLE_NAME:
            continue
        row = masked[row_idx]
        observed = (row != 78) & (row != 45)
        numer = int(observed.sum())
        prop = float(numer / locus_length) if locus_length else 0.0
        sample_props[sname] = prop
        if prop < min_sample_observed_fraction:
            masked[row_idx, :] = np.uint8(ord("N"))
            masked_samples.append(sname)
    return masked, masked_samples, sample_props


def resolve_locus_for_output(
    header: str,
    locus_dict: dict[str, str],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
    min_sample_observed_fraction: float = 0.10,
    *,
    forced_masked_samples: set[str] | None = None,
):
    """Trim one locus, apply sample-level masking, and return final summary metrics."""
    scaff, pos = header.split(":")
    rstart, rend = [int(i) for i in pos.split("-")]
    tnames = list(locus_dict.keys())
    empirical_row_mask = _empirical_sample_row_mask(tnames)
    seqs = np.array([list(bytes(seq, "utf-8")) for seq in locus_dict.values()], dtype=np.uint8)

    filters = {
        "min_length": False,
        "min_samples": False,
        "max_variant_frequency": False,
        "max_shared_hetero_frequency": False,
        "max_depth_outlier": False,
    }
    stats = {
        "samples_with_data_before_final_filters": int(
            np.sum(empirical_row_mask & _sample_rows_with_data(seqs))
        ),
        "locus_cov": 0,
        "variant_sites": 0,
        "variant_phylo_informative_sites": 0,
        "nsites": 0,
        "nsites_sample_cov_greater_than_1": 0,
        "nsites_sample_cov_greater_than_2": 0,
        "nsites_sample_cov_greater_than_3": 0,
        "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 0,
        "variant_site_frequency": 0,
        "variant_site_frequency_where_sample_cov_greater_than_2": 0,
        "variant_phylo_informative_site_frequency": 0,
        "variant_phylo_informative_site_frequency_where_sample_cov_greater_than_3": 0,
        "masked_samples_by_min_sample_observed_fraction": tuple(),
        "masked_sample_count_by_min_sample_observed_fraction": 0,
        "masked_samples_by_max_sample_hetero_frequency": tuple(),
        "masked_sample_count_by_max_sample_hetero_frequency": 0,
    }

    tseqs, tsite_sample_covs, trim_left, trim_right = _trim_locus_matrix(
        seqs,
        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
        empirical_row_mask=empirical_row_mask,
    )

    if forced_masked_samples is None:
        tseqs, masked_low_obs_samples, sample_observed_props = _mask_low_observed_samples(
            tnames,
            tseqs,
            min_sample_observed_fraction=min_sample_observed_fraction,
        )
        if masked_low_obs_samples:
            tseqs, tsite_sample_covs, extra_left, extra_right = _trim_locus_matrix(
                tseqs,
                min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
                empirical_row_mask=empirical_row_mask,
            )
            trim_left += extra_left
            trim_right += extra_right

        tseqs, masked_hetero_samples, sample_hetero_props = _mask_high_hetero_samples(
            tnames,
            tseqs,
            max_sample_hetero_frequency=max_sample_hetero_frequency,
        )
        if masked_hetero_samples:
            tseqs, tsite_sample_covs, extra_left, extra_right = _trim_locus_matrix(
                tseqs,
                min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
                empirical_row_mask=empirical_row_mask,
            )
            trim_left += extra_left
            trim_right += extra_right
        masked_samples = [*masked_low_obs_samples, *masked_hetero_samples]
    else:
        tseqs = tseqs.copy()
        masked_low_obs_samples = []
        masked_hetero_samples = []
        sample_observed_props = {}
        sample_hetero_props = {}
        masked_samples = []
        for row_idx, sname in enumerate(tnames):
            if sname in forced_masked_samples:
                tseqs[row_idx, :] = np.uint8(ord("N"))
                masked_samples.append(sname)

    stats["masked_samples_by_min_sample_observed_fraction"] = tuple(masked_low_obs_samples)
    stats["masked_sample_count_by_min_sample_observed_fraction"] = len(masked_low_obs_samples)
    stats["masked_samples_by_max_sample_hetero_frequency"] = tuple(masked_hetero_samples)
    stats["masked_sample_count_by_max_sample_hetero_frequency"] = len(masked_hetero_samples)
    if sample_observed_props:
        stats["sample_observed_fractions"] = sample_observed_props
    if sample_hetero_props:
        stats["sample_hetero_frequencies"] = sample_hetero_props

    row_has_data = empirical_row_mask & _sample_rows_with_data(tseqs)
    effective_locus_cov = int(row_has_data.sum())
    stats["locus_cov"] = effective_locus_cov
    if effective_locus_cov < min_locus_sample_coverage:
        filters["min_samples"] = True

    snpsarr = snp_count_numba(tseqs)
    stats["variant_sites"] = int(np.sum(snpsarr > 0))
    stats["variant_phylo_informative_sites"] = int(np.sum(snpsarr == 2))
    stats["nsites"] = int(tseqs.shape[1])
    stats["nsites_sample_cov_greater_than_1"] = int(np.sum(tsite_sample_covs > 1))
    stats["nsites_sample_cov_greater_than_2"] = int(np.sum(tsite_sample_covs > 2))
    stats["nsites_sample_cov_greater_than_3"] = int(np.sum(tsite_sample_covs > 3))
    stats["nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage"] = int(
        np.sum(tsite_sample_covs >= min_locus_trim_sample_coverage)
    )

    if stats["nsites_sample_cov_greater_than_2"]:
        stats["variant_site_frequency_where_sample_cov_greater_than_2"] = float(
            stats["variant_sites"] / stats["nsites_sample_cov_greater_than_2"]
        )

    if min_locus_sample_coverage >= 4:
        informative_sites = stats["nsites_sample_cov_greater_than_3"]
    elif min_locus_sample_coverage == 3:
        informative_sites = stats["nsites_sample_cov_greater_than_2"]
    else:
        informative_sites = stats["nsites_sample_cov_greater_than_1"]
    if informative_sites < min_locus_length:
        filters["min_length"] = True

    if stats["variant_site_frequency_where_sample_cov_greater_than_2"] > max_locus_variant_frequency:
        filters["max_variant_frequency"] = True

    empirical_tseqs = tseqs[empirical_row_mask]
    if empirical_tseqs.size:
        max_shared_h = max_heteros_count_numba(empirical_tseqs)
        max_shared_h_prop = max_shared_h / max(1, effective_locus_cov)
        if max_shared_h_prop > max_locus_hetero_frequency:
            filters["max_shared_hetero_frequency"] = True

    header = f"{scaff}:{rstart + trim_left}-{rend - trim_right}"
    return header, tnames, tseqs, snpsarr, filters, stats


def filter_trim_locus(
    header: str,
    locus_dict: dict[str, str],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 1.0,
    min_sample_observed_fraction: float = 0.10,
):
    """Trim one locus, evaluate the final filters, and return summary metrics."""
    return resolve_locus_for_output(
        header,
        locus_dict,
        min_locus_sample_coverage,
        min_locus_trim_sample_coverage,
        min_locus_length,
        max_locus_hetero_frequency,
        max_locus_variant_frequency,
        max_sample_hetero_frequency,
        min_sample_observed_fraction,
    )


def _safe_fraction(numer: int | float, denom: int | float) -> float:
    """Return numer/denom or 0.0 when the denominator is zero."""
    if not denom:
        return 0.0
    return float(numer / denom)


def _counter_mean(counter: Counter[int]) -> float:
    """Return the mean value represented by an integer frequency counter."""
    total = sum(counter.values())
    if not total:
        return 0.0
    return float(sum(value * count for value, count in counter.items()) / total)


def _counter_median(counter: Counter[int]) -> float:
    """Return the median value represented by an integer frequency counter."""
    total = sum(counter.values())
    if not total:
        return 0.0
    left_rank = (total - 1) // 2
    right_rank = total // 2
    cumulative = 0
    left_value = None
    right_value = None
    for value in sorted(counter):
        cumulative += counter[value]
        if left_value is None and cumulative > left_rank:
            left_value = value
        if right_value is None and cumulative > right_rank:
            right_value = value
            break
    if left_value is None or right_value is None:
        raise ValueError("Counter median could not be determined from the provided frequencies.")
    return float((left_value + right_value) / 2)


def _format_count(value: int) -> str:
    """Format integer counts consistently for text reports."""
    return f"{int(value):,}"


def _format_optional_count(value: int | None) -> str:
    """Format optional integer counts consistently for text reports."""
    if value is None:
        return "N/A"
    return _format_count(value)


def _format_float(value: float, digits: int = 3) -> str:
    """Format floating-point values for text reports."""
    return f"{float(value):.{digits}f}"


def _format_fraction(value: float) -> str:
    """Format fraction-like values consistently for text reports."""
    return f"{float(value):.6f}"


def _append_key_value_section(lines: list[str], title: str, rows: list[tuple[str, str]]) -> None:
    """Append one key/value report section to a list of output lines."""
    lines.append(f"# {title}")
    if rows:
        width = max(len(key) for key, _ in rows)
        for key, value in rows:
            lines.append(f"{key.ljust(width)}  {value}")
    lines.append("")


def _append_table_section(lines: list[str], title: str, headers: list[str], rows: list[list[str]]) -> None:
    """Append one simple whitespace-aligned table section."""
    lines.append(f"# {title}")
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    lines.append("  ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)))
    for row in rows:
        lines.append("  ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)))
    lines.append("")


def _human_stats_label(key: str) -> str:
    """Return the human-readable label for one stats key."""
    labels = {
        "samples": "Samples",
        "shared_loci_before_min_sample_coverage_filter": "Shared loci before minimum sample coverage filter",
        "shared_loci_after_delimiting": "Shared loci after delimiting",
        "shared_loci_after_paralog_filtering": "Shared loci after paralog filtering",
        "final_loci_written": "Final loci written",
        "final_loci_retained_fraction_of_post_paralog": "Final loci retained fraction after paralog filtering",
        "final_loci_retained_fraction_of_delimited": "Final loci retained fraction after delimiting",
        "assembled_sites": "Assembled sites",
        "final_snp_sites_written": "Final SNP sites written",
        "variable_sites": "Variable sites",
        "phylogenetically_informative_sites": "Phylogenetically informative sites",
        "alignment_matrix_occupancy_fraction": "Alignment matrix occupancy fraction",
        "overlapping_indel_clusters_masked": "Overlapping indel clusters masked",
        "overlapping_indel_records_removed": "Overlapping indel records removed",
        "overlapping_indel_bp_masked": "Overlapping indel bases masked",
        "loci_filtered_min_length": "Loci filtered by minimum length",
        "loci_filtered_min_samples": "Loci filtered by minimum sample coverage",
        "loci_filtered_max_variant_rate": "Loci filtered by maximum variant frequency",
        "loci_filtered_max_shared_heterozygosity": "Loci filtered by maximum shared heterozygosity",
        "loci_filtered_max_depth_outlier": "Loci filtered by maximum depth outlier",
        "loci_with_samples_masked_by_min_observed_fraction": "Loci with samples masked by minimum observed fraction threshold",
        "total_masked_sample_occurrences_by_min_observed_fraction": "Sample masks triggered by minimum observed fraction threshold",
        "loci_with_samples_masked_by_max_hetero_frequency": "Loci with samples masked by sample heterozygosity threshold",
        "total_masked_sample_occurrences_by_max_hetero_frequency": "Sample masks triggered by sample heterozygosity threshold",
        "mean_locus_length": "Mean locus length",
        "median_locus_length": "Median locus length",
        "min_locus_length": "Minimum locus length",
        "max_locus_length": "Maximum locus length",
        "mean_samples_per_locus": "Mean samples per locus",
        "median_samples_per_locus": "Median samples per locus",
        "sites_with_sample_coverage_ge_2": "Sites with sample coverage >= 2",
        "sites_with_sample_coverage_ge_3": "Sites with sample coverage >= 3",
        "sites_with_sample_coverage_ge_4": "Sites with sample coverage >= 4",
        "sites_with_sample_coverage_ge_trim_min": "Sites with sample coverage >= trim minimum",
        "sample": "Sample",
        "sample_type": "Sample type",
        "read_layout": "Read layout",
        "reads_before_filtering": "Reads before filtering",
        "reads_after_filtering": "Reads after filtering",
        "loci_in_alignment": "Loci in alignment",
        "loci_in_alignment_fraction": "Loci fraction in alignment",
        "shared_loci_with_nonzero_depth": "Shared loci with nonzero depth",
        "shared_loci_with_nonzero_depth_fraction": "Shared-depth loci fraction",
        "mean_depth_shared_loci": "Mean depth in shared loci",
        "median_depth_shared_loci": "Median depth in shared loci",
        "mean_depth_nonzero_shared_loci": "Mean depth in nonzero shared loci",
        "median_depth_nonzero_shared_loci": "Median depth in nonzero shared loci",
        "masked_by_min_observed_fraction": "Masked by minimum observed fraction threshold",
        "masked_by_max_hetero_frequency": "Masked by sample heterozygosity threshold",
        "samples_with_data": "Samples with data",
        "loci": "Loci",
        "rad_loci_before_min_sample_coverage": "RAD loci before min sample coverage",
        "rad_loci_after_min_sample_coverage": "RAD loci after min sample coverage",
        "loci_after_rad_wgs_integration": "Loci after RAD/WGS integration",
        "final_loci_after_filtering": "Final loci after filtering",
        "cumulative_final_loci": "Cumulative final loci",
        "fraction_of_final_loci": "Fraction of final loci",
        "rad_samples": "RAD samples",
        "wgs_samples": "WGS samples",
        "loci_fail_paralog_rad": "Loci failed by RAD paralog QC",
        "loci_fail_paralog_wgs": "Loci failed by WGS paralog QC",
        "loci_fail_paralog_both": "Loci failed by RAD and WGS paralog QC",
        "loci_pass_paralog_rad_fail_paralog_wgs": "Loci kept by RAD but failed by WGS QC",
        "sites_supported_rad_only": "Sites supported by RAD only",
        "sites_supported_wgs_only": "Sites supported by WGS only",
        "sites_supported_both": "Sites supported by RAD and WGS",
        "sites_supported_neither": "Sites supported by neither RAD nor WGS",
        "wgs_het_genotypes_masked_by_allele_balance": "WGS heterozygous genotypes masked by allele balance",
    }
    return labels.get(key, key.replace("_", " ").capitalize())


def write_assemble_stats_report(
    *,
    name: str,
    outdir: Path,
    logged_command: str | None = None,
    snames: list[str],
    sample_types: dict[str, str],
    sample_layouts: dict[str, str],
    sample_filter_stats: dict[str, dict[str, int]],
    shared_loci_before_min_sample_coverage_filter: int | None,
    shared_loci_after_delimiting: int,
    shared_loci_after_paralog_filtering: int,
    rad_locus_occupancy_before_min_sample_coverage_filter: dict[int, int] | None,
    rad_locus_occupancy_after_min_sample_coverage_filter: dict[int, int] | None,
    loci_summary: dict[str, object],
    sample_depth_stats: dict[str, dict[str, float]],
    nsnps_written: int,
    overlap_stats: dict[str, int],
    mixed_run_summary: dict[str, int] | None = None,
) -> Path:
    """Write the final human-readable assemble summary plus JSON sidecar."""
    outpath = outdir / f"{name}.stats.txt"
    json_path = outdir / f"{name}.stats.json"

    final_loci_written = int(loci_summary["nloci_after_filtering"])
    assembled_sites = int(loci_summary["nsites_after_filtering"])
    filter_counts = dict(loci_summary["filter_counts"])
    site_totals = dict(loci_summary["site_totals"])
    sample_locus_counts = dict(loci_summary["sample_locus_counts"])
    samples_per_locus_counts = Counter(loci_summary["samples_per_locus_counts"])
    integration_samples_per_locus_counts = Counter(
        loci_summary.get("samples_per_locus_before_final_filters_counts", {})
    )
    locus_length_counts = Counter(loci_summary["locus_length_counts"])
    alignment_nonmissing_sample_bases = int(loci_summary["alignment_nonmissing_sample_bases"])

    summary_data = {
        "samples": len(snames),
        "shared_loci_before_min_sample_coverage_filter": (
            int(shared_loci_before_min_sample_coverage_filter)
            if shared_loci_before_min_sample_coverage_filter is not None
            else None
        ),
        "shared_loci_after_delimiting": shared_loci_after_delimiting,
        "shared_loci_after_paralog_filtering": shared_loci_after_paralog_filtering,
        "final_loci_written": final_loci_written,
        "final_loci_retained_fraction_of_post_paralog": _safe_fraction(
            final_loci_written, shared_loci_after_paralog_filtering
        ),
        "final_loci_retained_fraction_of_delimited": _safe_fraction(
            final_loci_written, shared_loci_after_delimiting
        ),
        "assembled_sites": assembled_sites,
        "final_snp_sites_written": nsnps_written,
        "variable_sites": int(site_totals["variant_sites"]),
        "phylogenetically_informative_sites": int(
            site_totals["variant_phylo_informative_sites"]
        ),
        "alignment_matrix_occupancy_fraction": _safe_fraction(
            alignment_nonmissing_sample_bases,
            assembled_sites * len(snames),
        ),
        "overlapping_indel_clusters_masked": int(
            overlap_stats["overlapping_indel_clusters_masked"]
        ),
        "overlapping_indel_records_removed": int(
            overlap_stats["overlapping_indel_records_removed"]
        ),
        "overlapping_indel_bp_masked": int(
            overlap_stats["overlapping_indel_bp_masked"]
        ),
    }
    filtering_data = {
        "loci_filtered_min_length": int(filter_counts["min_length"]),
        "loci_filtered_min_samples": (
            max(
                0,
                int(shared_loci_before_min_sample_coverage_filter)
                - int(shared_loci_after_delimiting),
            )
            if shared_loci_before_min_sample_coverage_filter is not None
            else None
        ),
        "loci_filtered_max_variant_rate": int(filter_counts["max_variant_frequency"]),
        "loci_filtered_max_shared_heterozygosity": int(
            filter_counts["max_shared_hetero_frequency"]
        ),
        "loci_filtered_max_depth_outlier": int(filter_counts["max_depth_outlier"]),
    }
    sample_masking_data = {
        "loci_with_samples_masked_by_min_observed_fraction": int(
            loci_summary.get("loci_with_samples_masked_by_min_observed_fraction", 0)
        ),
        "total_masked_sample_occurrences_by_min_observed_fraction": int(
            loci_summary.get("total_masked_sample_occurrences_by_min_observed_fraction", 0)
        ),
        "loci_with_samples_masked_by_max_hetero_frequency": int(
            loci_summary.get("loci_with_samples_masked_by_max_hetero_frequency", 0)
        ),
        "total_masked_sample_occurrences_by_max_hetero_frequency": int(
            loci_summary.get("total_masked_sample_occurrences_by_max_hetero_frequency", 0)
        ),
    }
    alignment_data = {
        "mean_locus_length": float(_counter_mean(locus_length_counts)),
        "median_locus_length": float(_counter_median(locus_length_counts)),
        "min_locus_length": int(min(locus_length_counts) if locus_length_counts else 0),
        "max_locus_length": int(max(locus_length_counts) if locus_length_counts else 0),
        "mean_samples_per_locus": float(_counter_mean(samples_per_locus_counts)),
        "median_samples_per_locus": float(_counter_median(samples_per_locus_counts)),
        "sites_with_sample_coverage_ge_2": int(
            site_totals["nsites_sample_cov_greater_than_1"]
        ),
        "sites_with_sample_coverage_ge_3": int(
            site_totals["nsites_sample_cov_greater_than_2"]
        ),
        "sites_with_sample_coverage_ge_4": int(
            site_totals["nsites_sample_cov_greater_than_3"]
        ),
        "sites_with_sample_coverage_ge_trim_min": int(
            site_totals[
                "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage"
            ]
        ),
    }

    sample_headers = [
        "sample",
        "sample_type",
        "read_layout",
        "reads_before_filtering",
        "reads_after_filtering",
        "loci_in_alignment",
        "loci_in_alignment_fraction",
        "shared_loci_with_nonzero_depth",
        "shared_loci_with_nonzero_depth_fraction",
        "mean_depth_shared_loci",
        "median_depth_shared_loci",
        "mean_depth_nonzero_shared_loci",
        "median_depth_nonzero_shared_loci",
        "masked_by_min_observed_fraction",
        "masked_by_max_hetero_frequency",
    ]
    sample_summary_data: list[dict[str, object]] = []
    sample_rows: list[list[str]] = []
    for sname in sorted(snames):
        depth_stats = sample_depth_stats[sname]
        loci_in_alignment = int(sample_locus_counts.get(sname, 0))
        nonzero_loci = int(depth_stats["shared_loci_with_nonzero_depth"])
        sample_record = {
            "sample": sname,
            "sample_type": sample_types.get(sname, ""),
            "read_layout": sample_layouts.get(sname, ""),
            "reads_before_filtering": int(
                sample_filter_stats.get(sname, {}).get("reads_before_filtering", 0)
            ),
            "reads_after_filtering": int(
                sample_filter_stats.get(sname, {}).get("reads_after_filtering", 0)
            ),
            "loci_in_alignment": loci_in_alignment,
            "loci_in_alignment_fraction": _safe_fraction(
                loci_in_alignment, final_loci_written
            ),
            "shared_loci_with_nonzero_depth": nonzero_loci,
            "shared_loci_with_nonzero_depth_fraction": _safe_fraction(
                nonzero_loci, final_loci_written
            ),
            "mean_depth_shared_loci": float(depth_stats["mean_depth_shared_loci"]),
            "median_depth_shared_loci": float(depth_stats["median_depth_shared_loci"]),
            "mean_depth_nonzero_shared_loci": float(
                depth_stats["mean_depth_nonzero_shared_loci"]
            ),
            "median_depth_nonzero_shared_loci": float(
                depth_stats["median_depth_nonzero_shared_loci"]
            ),
            "masked_by_min_observed_fraction": int(
                loci_summary.get("masked_by_min_observed_fraction_counts", {}).get(
                    sname, 0
                )
            ),
            "masked_by_max_hetero_frequency": int(
                loci_summary.get("masked_by_max_hetero_frequency_counts", {}).get(
                    sname, 0
                )
            ),
        }
        sample_summary_data.append(sample_record)
        sample_rows.append([
            sample_record["sample"],
            sample_record["sample_type"],
            sample_record["read_layout"],
            _format_count(sample_record["reads_before_filtering"]),
            _format_count(sample_record["reads_after_filtering"]),
            _format_count(sample_record["loci_in_alignment"]),
            _format_fraction(sample_record["loci_in_alignment_fraction"]),
            _format_count(sample_record["shared_loci_with_nonzero_depth"]),
            _format_fraction(sample_record["shared_loci_with_nonzero_depth_fraction"]),
            _format_float(sample_record["mean_depth_shared_loci"]),
            _format_float(sample_record["median_depth_shared_loci"]),
            _format_float(sample_record["mean_depth_nonzero_shared_loci"]),
            _format_float(sample_record["median_depth_nonzero_shared_loci"]),
            _format_count(sample_record["masked_by_min_observed_fraction"]),
            _format_count(sample_record["masked_by_max_hetero_frequency"]),
        ])

    occupancy_headers = [
        "samples_with_data",
        "rad_loci_before_min_sample_coverage",
        "rad_loci_after_min_sample_coverage",
        "loci_after_rad_wgs_integration",
        "final_loci_after_filtering",
        "cumulative_final_loci",
        "fraction_of_final_loci",
    ]
    pre_min_occupancy_counts = {
        int(key): int(value)
        for key, value in (rad_locus_occupancy_before_min_sample_coverage_filter or {}).items()
    }
    post_min_occupancy_counts = {
        int(key): int(value)
        for key, value in (rad_locus_occupancy_after_min_sample_coverage_filter or {}).items()
    }
    cumulative_final_loci = 0
    cumulative_final_loci_by_sample_count: dict[int, int] = {}
    for sample_count in range(len(snames) + 1):
        cumulative_final_loci += int(samples_per_locus_counts.get(sample_count, 0))
        cumulative_final_loci_by_sample_count[sample_count] = cumulative_final_loci
    locus_occupancy_data = [
        {
            "samples_with_data": int(sample_count),
            "rad_loci_before_min_sample_coverage": (
                int(pre_min_occupancy_counts.get(sample_count, 0))
                if rad_locus_occupancy_before_min_sample_coverage_filter is not None
                else None
            ),
            "rad_loci_after_min_sample_coverage": (
                int(post_min_occupancy_counts.get(sample_count, 0))
                if rad_locus_occupancy_after_min_sample_coverage_filter is not None
                else None
            ),
            "loci_after_rad_wgs_integration": int(
                integration_samples_per_locus_counts.get(sample_count, 0)
            ),
            "final_loci_after_filtering": int(
                samples_per_locus_counts.get(sample_count, 0)
            ),
            "cumulative_final_loci": int(
                cumulative_final_loci_by_sample_count.get(sample_count, 0)
            ),
            "fraction_of_final_loci": _safe_fraction(
                int(samples_per_locus_counts.get(sample_count, 0)),
                final_loci_written,
            ),
        }
        for sample_count in range(len(snames) + 1)
    ]
    occupancy_rows = [
        [
            _format_count(row["samples_with_data"]),
            _format_optional_count(row["rad_loci_before_min_sample_coverage"]),
            _format_optional_count(row["rad_loci_after_min_sample_coverage"]),
            _format_count(row["loci_after_rad_wgs_integration"]),
            _format_count(row["final_loci_after_filtering"]),
            _format_count(row["cumulative_final_loci"]),
            _format_fraction(row["fraction_of_final_loci"]),
        ]
        for row in locus_occupancy_data
    ]

    mixed_data = None
    if mixed_run_summary:
        mixed_data = {
            key: int(value)
            for key, value in mixed_run_summary.items()
        }

    stats_json: dict[str, object] = {
        "summary": summary_data,
        "locus_filtering": filtering_data,
        "sample_masking": sample_masking_data,
        "alignment_summary": alignment_data,
        "sample_summary": sample_summary_data,
        "locus_occupancy": locus_occupancy_data,
    }
    if logged_command:
        stats_json["command"] = logged_command
    if mixed_data is not None:
        stats_json["mixed_rad_wgs_diagnostics"] = mixed_data

    lines: list[str] = []
    _append_key_value_section(
        lines,
        "Assemble Summary",
        [
            (_human_stats_label("samples"), _format_count(summary_data["samples"])),
            (
                _human_stats_label("shared_loci_before_min_sample_coverage_filter"),
                _format_optional_count(
                    summary_data["shared_loci_before_min_sample_coverage_filter"]
                ),
            ),
            (
                _human_stats_label("shared_loci_after_delimiting"),
                _format_count(summary_data["shared_loci_after_delimiting"]),
            ),
            (
                _human_stats_label("shared_loci_after_paralog_filtering"),
                _format_count(summary_data["shared_loci_after_paralog_filtering"]),
            ),
            (
                _human_stats_label("final_loci_written"),
                _format_count(summary_data["final_loci_written"]),
            ),
            (
                _human_stats_label("final_loci_retained_fraction_of_post_paralog"),
                _format_fraction(
                    summary_data["final_loci_retained_fraction_of_post_paralog"]
                ),
            ),
            (
                _human_stats_label("final_loci_retained_fraction_of_delimited"),
                _format_fraction(
                    summary_data["final_loci_retained_fraction_of_delimited"]
                ),
            ),
            (
                _human_stats_label("assembled_sites"),
                _format_count(summary_data["assembled_sites"]),
            ),
            (
                _human_stats_label("final_snp_sites_written"),
                _format_count(summary_data["final_snp_sites_written"]),
            ),
            (
                _human_stats_label("variable_sites"),
                _format_count(summary_data["variable_sites"]),
            ),
            (
                _human_stats_label("phylogenetically_informative_sites"),
                _format_count(summary_data["phylogenetically_informative_sites"]),
            ),
            (
                _human_stats_label("alignment_matrix_occupancy_fraction"),
                _format_fraction(summary_data["alignment_matrix_occupancy_fraction"]),
            ),
            (
                _human_stats_label("overlapping_indel_clusters_masked"),
                _format_count(summary_data["overlapping_indel_clusters_masked"]),
            ),
            (
                _human_stats_label("overlapping_indel_records_removed"),
                _format_count(summary_data["overlapping_indel_records_removed"]),
            ),
            (
                _human_stats_label("overlapping_indel_bp_masked"),
                _format_count(summary_data["overlapping_indel_bp_masked"]),
            ),
        ],
    )
    if mixed_data is not None:
        mixed_rows = [
            (_human_stats_label(key), _format_count(value))
            for key, value in mixed_data.items()
            if key
            in {
                "rad_samples",
                "wgs_samples",
                "loci_fail_paralog_rad",
                "loci_fail_paralog_wgs",
                "loci_fail_paralog_both",
                "loci_pass_paralog_rad_fail_paralog_wgs",
                "sites_supported_rad_only",
                "sites_supported_wgs_only",
                "sites_supported_both",
                "sites_supported_neither",
                "wgs_het_genotypes_masked_by_allele_balance",
            }
        ]
        _append_key_value_section(lines, "Mixed RAD/WGS Diagnostics", mixed_rows)
    _append_key_value_section(
        lines,
        "Locus Filtering",
        [
            (
                _human_stats_label(key),
                _format_optional_count(value)
                if key == "loci_filtered_min_samples"
                else _format_count(value),
            )
            for key, value in filtering_data.items()
        ],
    )
    _append_key_value_section(
        lines,
        "Sample Masking",
        [(_human_stats_label(key), _format_count(value)) for key, value in sample_masking_data.items()],
    )
    _append_key_value_section(
        lines,
        "Alignment Summary",
        [
            (_human_stats_label("mean_locus_length"), _format_float(alignment_data["mean_locus_length"])),
            (_human_stats_label("median_locus_length"), _format_float(alignment_data["median_locus_length"])),
            (_human_stats_label("min_locus_length"), _format_count(alignment_data["min_locus_length"])),
            (_human_stats_label("max_locus_length"), _format_count(alignment_data["max_locus_length"])),
            (_human_stats_label("mean_samples_per_locus"), _format_float(alignment_data["mean_samples_per_locus"])),
            (_human_stats_label("median_samples_per_locus"), _format_float(alignment_data["median_samples_per_locus"])),
            (_human_stats_label("sites_with_sample_coverage_ge_2"), _format_count(alignment_data["sites_with_sample_coverage_ge_2"])),
            (_human_stats_label("sites_with_sample_coverage_ge_3"), _format_count(alignment_data["sites_with_sample_coverage_ge_3"])),
            (_human_stats_label("sites_with_sample_coverage_ge_4"), _format_count(alignment_data["sites_with_sample_coverage_ge_4"])),
            (_human_stats_label("sites_with_sample_coverage_ge_trim_min"), _format_count(alignment_data["sites_with_sample_coverage_ge_trim_min"])),
        ],
    )
    _append_table_section(
        lines,
        "Sample Summary",
        [_human_stats_label(header) for header in sample_headers],
        sample_rows,
    )
    _append_table_section(
        lines,
        "Locus Occupancy",
        [_human_stats_label(header) for header in occupancy_headers],
        occupancy_rows,
    )
    report_text = "\n".join(lines).rstrip() + "\n"
    if logged_command:
        report_text = f"CMD: {logged_command}\n\n{report_text}"
    outpath.write_text(report_text, encoding="utf-8")
    json_path.write_text(json.dumps(stats_json, indent=2) + "\n", encoding="utf-8")
    logger.info("wrote assemble summary report")
    logger.debug("assemble stats written to {} and {}", outpath, json_path)
    return outpath


def _build_retained_locus_outputs(
    *,
    raw_header: str,
    header: str,
    tnames: list[str],
    tseqs: np.ndarray,
    snpsarr: np.ndarray,
    stats: dict[str, object],
    padded: dict[str, str],
    refname: str,
) -> dict[str, object]:
    """Build all retained-locus outputs needed by the ordered final writer."""
    scaff, pos = header.split(":")
    pos0, pos1 = (int(i) for i in pos.split("-"))
    masked_low_observed_samples = set(
        stats.get("masked_samples_by_min_sample_observed_fraction", ())
    )
    masked_hetero_samples = set(
        stats.get("masked_samples_by_max_sample_hetero_frequency", ())
    )
    masked_samples = masked_low_observed_samples | masked_hetero_samples

    sample_mask = _empirical_sample_row_mask(tnames, refname=refname)
    sample_rows = sample_mask & _sample_rows_with_data(tseqs)
    sample_names_with_data = [
        sname
        for row_idx, sname in enumerate(tnames)
        if sname != refname and sample_rows[row_idx]
    ]
    sample_count = len(sample_names_with_data)
    nonmissing_sample_bases = 0
    if sample_count:
        sample_seqs = tseqs[sample_rows]
        nonmissing_sample_bases = int(np.sum((sample_seqs != 78) & (sample_seqs != 45)))

    # `.loci.gz` should omit samples masked by any final sample-row filter,
    # while fixed-axis outputs keep them as missing data for that locus.
    visible_rows = [
        (sname, seq)
        for sname, seq in zip(tnames, tseqs, strict=True)
        if sname not in masked_samples
    ]
    visible_empirical_row_count = sum(1 for sname, _seq in visible_rows if sname != refname)
    locus_lines = [f"{padded[sname]}{bytes(seq).decode()}" for sname, seq in visible_rows]
    snpstring_arr = snpsarr.copy()
    snpstring_arr[snpstring_arr == 0] = 32
    snpstring_arr[snpstring_arr == 1] = 45
    snpstring_arr[snpstring_arr == 2] = 42

    return {
        "raw_header": raw_header,
        "header": header,
        "tseqs": tseqs,
        "tnames": tnames,
        "locus_length": int(tseqs.shape[1]),
        "bed_row": f"{scaff}\t{pos0 - 1}\t{pos1}\t{visible_empirical_row_count}\n",
        "manifest_row": (raw_header, header, ",".join(sorted(masked_samples))),
        "masked_low_observed_samples": sorted(masked_low_observed_samples),
        "masked_hetero_samples": sorted(masked_hetero_samples),
        "masked_samples": sorted(masked_samples),
        "mask_bed_row": f"{scaff}\t{pos0 - 1}\t{pos1}\n",
        "sample_names_with_data": sample_names_with_data,
        "sample_count": sample_count,
        "nonmissing_sample_bases": nonmissing_sample_bases,
        "stats": dict(stats),
        "locus_lines": locus_lines,
        "snpstring": bytes(snpstring_arr).decode(),
        "scaff": scaff,
        "pos0": pos0,
        "pos1": pos1,
    }


def _resolve_output_batch(
    *,
    batch_idx: int,
    batch_items: list[tuple[str, dict[str, str]]],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float,
    min_sample_observed_fraction: float,
    padded: dict[str, str],
    refname: str,
    output_snames: list[str],
    scaff2idx: dict[str, int],
    map_dtype_name: str,
) -> tuple[int, dict[str, object]]:
    """Resolve one batch of loci for the combined final-output writer."""
    map_dtype = np.dtype(map_dtype_name)
    sidxs = {name: idx for idx, name in enumerate(output_snames)}
    nsamples = len(output_snames)
    retained_loci: list[dict[str, object]] = []
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    samples_per_locus_before_final_filters = Counter()
    phy_cursor = 0
    chunk_arrays: list[np.ndarray] = []
    chunk_map_rows: list[tuple[int, int, int, int, int]] = []

    for raw_header, ldict in batch_items:
        header, tnames, tseqs, snpsarr, filters, stats = resolve_locus_for_output(
            raw_header,
            ldict,
            min_locus_sample_coverage,
            min_locus_trim_sample_coverage,
            min_locus_length,
            max_locus_hetero_frequency,
            max_locus_variant_frequency,
            max_sample_hetero_frequency,
            min_sample_observed_fraction,
        )
        for key in total_filters:
            total_filters[key] += int(filters[key])
        samples_per_locus_before_final_filters[
            int(stats["samples_with_data_before_final_filters"])
        ] += 1
        if sum(filters.values()):
            continue

        locus_output = _build_retained_locus_outputs(
            raw_header=raw_header,
            header=header,
            tnames=tnames,
            tseqs=tseqs,
            snpsarr=snpsarr,
            stats=stats,
            padded=padded,
            refname=refname,
        )
        locus_arr = np.full((nsamples, locus_output["locus_length"]), np.uint8(ord("N")), dtype=np.uint8)
        for row_idx, sname in enumerate(tnames):
            locus_arr[sidxs[sname], :] = tseqs[row_idx]
        chunk_arrays.append(locus_arr)
        chunk_map_rows.append(
            (
                int(scaff2idx[locus_output["scaff"]]),
                phy_cursor,
                phy_cursor + int(locus_output["locus_length"]),
                int(locus_output["pos0"]),
                int(locus_output["pos1"]),
            )
        )
        phy_cursor += int(locus_output["locus_length"])
        retained_loci.append(locus_output)

    if chunk_arrays:
        chunkarr = np.concatenate(chunk_arrays, axis=1)
    else:
        chunkarr = np.empty((nsamples, 0), dtype=np.uint8)
    if chunk_map_rows:
        chunkmap = np.array(chunk_map_rows, dtype=map_dtype)
    else:
        chunkmap = np.empty((0, 5), dtype=map_dtype)

    return {
        "nloci_before_filtering": len(batch_items),
        "filter_counts": total_filters,
        "samples_per_locus_before_final_filters_counts": dict(
            samples_per_locus_before_final_filters
        ),
        "retained_loci": retained_loci,
        "chunkarr": chunkarr,
        "chunkmap": chunkmap,
    }


def write_final_outputs(
    *,
    snames: list[str],
    name: str,
    outdir: Path,
    reference: Path,
    database_fasta: Path,
    retained_loci_manifest: Path,
    consensus_sample_mask_beds: dict[str, Path],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
    min_sample_observed_fraction: float = 0.10,
    cores: int = 1,
    log_level: str = "INFO",
    batch_size: int = 128,
) -> dict[str, object]:
    """Resolve final loci once, then write `.loci`, BED, stats, and sequence HDF5."""
    from .write_seqs import (
        append_seqs_hdf5_chunk,
        close_seqs_hdf5_writer,
        finalize_seqs_hdf5_writer,
        open_seqs_hdf5_writer,
    )

    loci_file = outdir / f"{name}.loci.gz"
    final_loci_bed = outdir / f"{name}.bed"
    for path in consensus_sample_mask_beds.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    real_snames = list(snames)
    refname = REFERENCE_SAMPLE_NAME
    output_snames = [refname] + sorted(real_snames)
    max_len = max(len(i) for i in output_snames) + 2
    padded = {n: n + (" " * (max_len - len(n))) for n in output_snames}

    if final_loci_bed.exists():
        final_loci_bed.unlink()
    writer = open_seqs_hdf5_writer(
        name=name,
        outdir=outdir,
        snames=real_snames,
        reference=reference,
        loci_bed=final_loci_bed,
    )
    scaff2idx = {str(scaff): idx for idx, scaff in enumerate(writer.io5.attrs["scaffold_names"])}

    samples_per_locus = Counter()
    samples_per_locus_before_final_filters = Counter()
    locus_length_counts = Counter()
    per_sample_locus_counts = {i: 0 for i in real_snames}
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    masked_by_min_observed_fraction_counts = {i: 0 for i in real_snames}
    masked_by_max_hetero_frequency_counts = {i: 0 for i in real_snames}
    loci_with_samples_masked_by_min_observed_fraction = 0
    loci_with_samples_masked_by_max_hetero_frequency = 0
    total_masked_sample_occurrences_by_min_observed_fraction = 0
    total_masked_sample_occurrences_by_max_hetero_frequency = 0
    total_stats = {
        "variant_sites": 0,
        "variant_phylo_informative_sites": 0,
        "nsites": 0,
        "nsites_sample_cov_greater_than_1": 0,
        "nsites_sample_cov_greater_than_2": 0,
        "nsites_sample_cov_greater_than_3": 0,
        "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 0,
    }
    alignment_nonmissing_sample_bases = 0
    nloci_before_filtering = 0
    flidx = 0
    retained_scaffold_names_seen: set[str] = set()

    nbatches = sum(1 for _ in iter_locus_batches(database_fasta, batch_size=batch_size))
    output_workers = min(max(1, cores - 1), max(1, nbatches))

    def _iter_jobs():
        for batch_idx, batch_items in iter_locus_batches(database_fasta, batch_size=batch_size):
            yield batch_idx, (
                _resolve_output_batch,
                dict(
                    batch_idx=batch_idx,
                    batch_items=batch_items,
                    min_locus_sample_coverage=min_locus_sample_coverage,
                    min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
                    min_locus_length=min_locus_length,
                    max_locus_hetero_frequency=max_locus_hetero_frequency,
                    max_locus_variant_frequency=max_locus_variant_frequency,
                    max_sample_hetero_frequency=max_sample_hetero_frequency,
                    min_sample_observed_fraction=min_sample_observed_fraction,
                    padded=padded,
                    refname=refname,
                    output_snames=output_snames,
                    scaff2idx=scaff2idx,
                    map_dtype_name=writer.phymap.dtype.name,
                ),
            )

    pending_results: dict[int, dict[str, object]] = {}
    next_batch_idx = 0
    with gzip.open(loci_file, "wt", encoding="utf-8", compresslevel=1) as out, final_loci_bed.open(
        "w",
        encoding="utf-8",
    ) as out_bed, retained_loci_manifest.open("w", encoding="utf-8", newline="") as manifest_handle:
        manifest_writer = csv.writer(manifest_handle, delimiter="\t")
        manifest_writer.writerow(["raw_header", "final_header", "masked_samples"])
        mask_handles = {
            sname: consensus_sample_mask_beds[sname].open("w", encoding="utf-8")
            for sname in real_snames
        }

        def _flush_pending() -> None:
            nonlocal next_batch_idx, nloci_before_filtering, flidx
            nonlocal loci_with_samples_masked_by_min_observed_fraction
            nonlocal loci_with_samples_masked_by_max_hetero_frequency
            nonlocal total_masked_sample_occurrences_by_min_observed_fraction
            nonlocal total_masked_sample_occurrences_by_max_hetero_frequency
            nonlocal alignment_nonmissing_sample_bases
            while next_batch_idx in pending_results:
                batch_result = pending_results.pop(next_batch_idx)
                next_batch_idx += 1
                nloci_before_filtering += int(batch_result["nloci_before_filtering"])
                for key in total_filters:
                    total_filters[key] += int(batch_result["filter_counts"][key])
                samples_per_locus_before_final_filters.update(
                    batch_result[
                        "samples_per_locus_before_final_filters_counts"
                    ]
                )
                retained_loci = batch_result["retained_loci"]
                for locus_output in retained_loci:
                    out_bed.write(str(locus_output["bed_row"]))
                    retained_scaffold_names_seen.add(str(locus_output["scaff"]))
                    manifest_writer.writerow(list(locus_output["manifest_row"]))
                    masked_low_obs_samples = list(locus_output["masked_low_observed_samples"])
                    if masked_low_obs_samples:
                        loci_with_samples_masked_by_min_observed_fraction += 1
                        total_masked_sample_occurrences_by_min_observed_fraction += len(
                            masked_low_obs_samples
                        )
                        for sname in masked_low_obs_samples:
                            masked_by_min_observed_fraction_counts[sname] += 1
                    masked_hetero_samples = list(locus_output["masked_hetero_samples"])
                    masked_samples = list(locus_output["masked_samples"])
                    if masked_hetero_samples:
                        loci_with_samples_masked_by_max_hetero_frequency += 1
                        total_masked_sample_occurrences_by_max_hetero_frequency += len(
                            masked_hetero_samples
                        )
                        for sname in masked_hetero_samples:
                            masked_by_max_hetero_frequency_counts[sname] += 1
                    if masked_samples:
                        for sname in masked_samples:
                            mask_handles[sname].write(str(locus_output["mask_bed_row"]))

                    for sname in locus_output["sample_names_with_data"]:
                        per_sample_locus_counts[sname] += 1
                    samples_per_locus[int(locus_output["sample_count"])] += 1
                    locus_length_counts[int(locus_output["locus_length"])] += 1
                    alignment_nonmissing_sample_bases += int(locus_output["nonmissing_sample_bases"])
                    for stat in total_stats:
                        total_stats[stat] += int(locus_output["stats"][stat])

                    out.write(
                        "\n".join(locus_output["locus_lines"])
                        + f"\n//{' ' * (max_len - 2)}{locus_output['snpstring']}|{flidx}:{locus_output['header']}\n"
                    )
                    flidx += 1

                chunkarr = batch_result["chunkarr"]
                chunkmap = batch_result["chunkmap"]
                if chunkmap.size:
                    chunkmap = chunkmap.copy()
                    chunkmap[:, 1] += writer.nsites
                    chunkmap[:, 2] += writer.nsites
                append_seqs_hdf5_chunk(writer, chunkarr, chunkmap)

        try:
            if nbatches == 0:
                pass
            elif output_workers == 1 or nbatches == 1:
                for batch_idx, batch_items in iter_locus_batches(database_fasta, batch_size=batch_size):
                    result = _resolve_output_batch(
                        batch_idx=batch_idx,
                        batch_items=batch_items,
                        min_locus_sample_coverage=min_locus_sample_coverage,
                        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
                        min_locus_length=min_locus_length,
                        max_locus_hetero_frequency=max_locus_hetero_frequency,
                        max_locus_variant_frequency=max_locus_variant_frequency,
                        max_sample_hetero_frequency=max_sample_hetero_frequency,
                        min_sample_observed_fraction=min_sample_observed_fraction,
                        padded=padded,
                        refname=refname,
                        output_snames=output_snames,
                        scaff2idx=scaff2idx,
                        map_dtype_name=writer.phymap.dtype.name,
                    )
                    pending_results[batch_idx] = result
                    _flush_pending()
            else:
                for batch_idx, result in run_with_pool_iter(
                    _iter_jobs(),
                    log_level=log_level,
                    max_workers=output_workers,
                    msg="Resolving and writing final loci",
                    njobs=nbatches,
                ):
                    pending_results[batch_idx] = result
                    _flush_pending()
        except Exception:
            for handle in mask_handles.values():
                handle.close()
            close_seqs_hdf5_writer(writer)
            raise
        else:
            for handle in mask_handles.values():
                handle.close()

    all_scaffold_names = [str(value) for value in writer.io5.attrs["scaffold_names"]]
    all_scaffold_lengths = [int(value) for value in writer.io5.attrs["scaffold_lengths"]]
    retained_scaffold_names = [
        name
        for name in all_scaffold_names
        if name in retained_scaffold_names_seen
    ]
    length_lookup = dict(zip(all_scaffold_names, all_scaffold_lengths, strict=True))
    retained_scaffold_lengths = [length_lookup[name] for name in retained_scaffold_names]
    try:
        finalize_seqs_hdf5_writer(
            writer,
            expected_nsites=int(total_stats["nsites"]),
            expected_nloci=flidx,
            retained_scaffold_names=retained_scaffold_names,
            retained_scaffold_lengths=retained_scaffold_lengths,
        )
    finally:
        close_seqs_hdf5_writer(writer)

    logger.debug("wrote final loci, BED, and sequence HDF5 outputs to {}", outdir)
    return {
        "nloci_before_filtering": nloci_before_filtering,
        "nloci_after_filtering": flidx,
        "nsites_after_filtering": int(total_stats["nsites"]),
        "filter_counts": dict(total_filters),
        "site_totals": dict(total_stats),
        "sample_locus_counts": dict(per_sample_locus_counts),
        "samples_per_locus_before_final_filters_counts": dict(
            samples_per_locus_before_final_filters
        ),
        "masked_by_min_observed_fraction_counts": dict(masked_by_min_observed_fraction_counts),
        "loci_with_samples_masked_by_min_observed_fraction": loci_with_samples_masked_by_min_observed_fraction,
        "total_masked_sample_occurrences_by_min_observed_fraction": total_masked_sample_occurrences_by_min_observed_fraction,
        "masked_by_max_hetero_frequency_counts": dict(masked_by_max_hetero_frequency_counts),
        "loci_with_samples_masked_by_max_hetero_frequency": loci_with_samples_masked_by_max_hetero_frequency,
        "total_masked_sample_occurrences_by_max_hetero_frequency": total_masked_sample_occurrences_by_max_hetero_frequency,
        "samples_per_locus_counts": dict(samples_per_locus),
        "locus_length_counts": dict(locus_length_counts),
        "alignment_nonmissing_sample_bases": alignment_nonmissing_sample_bases,
    }


def write_loci_and_stats_files(
    snames: list[str],
    name: str,
    outdir: Path,
    tmpdir: Path,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
    min_sample_observed_fraction: float = 0.10,
)-> dict[str, object]:
    """Write the final `.loci.gz` and `.bed` outputs and collect report counters."""
    # database file is in the tmpdir inside outdir
    database = tmpdir / f"{name}.database.fa"
    loci_file = outdir / f"{name}.loci.gz"
    loci_compresslevel = 1

    # Add the reference as a synthetic sample in the written database/stats
    # without mutating the caller-owned sample-name list.
    real_snames = list(snames)
    refname = REFERENCE_SAMPLE_NAME
    write_snames = real_snames + [refname]

    # get name padding for loci file
    max_len = max(len(i) for i in write_snames) + 2
    padded = {n: n + (" " * (max_len - len(n))) for n in write_snames}

    # Collect final-locus summary counters in the same streaming pass that
    # writes the filtered loci and BED outputs, so the final stats report does
    # not need to reparse the full database later.
    samples_per_locus = Counter()
    samples_per_locus_before_final_filters = Counter()
    locus_length_counts = Counter()
    per_sample_locus_counts = {i: 0 for i in real_snames}
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    masked_by_min_observed_fraction_counts = {i: 0 for i in real_snames}
    masked_by_max_hetero_frequency_counts = {i: 0 for i in real_snames}
    loci_with_samples_masked_by_min_observed_fraction = 0
    loci_with_samples_masked_by_max_hetero_frequency = 0
    total_masked_sample_occurrences_by_min_observed_fraction = 0
    total_masked_sample_occurrences_by_max_hetero_frequency = 0
    total_stats = {
        "variant_sites": 0,
        "variant_phylo_informative_sites": 0,
        "nsites": 0,
        "nsites_sample_cov_greater_than_1": 0,
        "nsites_sample_cov_greater_than_2": 0,
        "nsites_sample_cov_greater_than_3": 0,
        "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 0,
    }
    alignment_nonmissing_sample_bases = 0

    # Build the final loci and BED outputs in one streaming pass over the
    # consensus database so large assemblies do not accumulate all loci text or
    # all retained BED rows in memory before writing them.
    loci = []
    lidx = 0    # counter of all loci
    flidx = 0   # counter of loci that passed filters
    manifest_path = get_retained_loci_manifest_path(name, tmpdir)
    (tmpdir / "beds").mkdir(parents=True, exist_ok=True)
    with gzip.open(loci_file, "wt", encoding="utf-8", compresslevel=loci_compresslevel) as out, open(
        outdir / f"{name}.bed",
        "w",
        encoding="utf-8",
    ) as out_bed, manifest_path.open("w", encoding="utf-8", newline="") as manifest_handle:
        manifest_writer = csv.writer(manifest_handle, delimiter="\t")
        manifest_writer.writerow(["raw_header", "final_header", "masked_samples"])
        mask_handles = {
            sname: get_consensus_sample_mask_path(sname, tmpdir).open("w", encoding="utf-8")
            for sname in real_snames
        }
        for oheader, ldict in iter_parse_loci(database):

            # apply trim and filters to locus
            args = (
                oheader,
                ldict,
                min_locus_sample_coverage,
                min_locus_trim_sample_coverage,
                min_locus_length,
                max_locus_hetero_frequency,
                max_locus_variant_frequency,
                max_sample_hetero_frequency,
                min_sample_observed_fraction,
            )
            result = filter_trim_locus(*args)
            header, tnames, tseqs, snpsarr, filters, stats = result

            # update total dicts
            for key in total_filters:
                total_filters[key] += int(result[4][key])
            samples_per_locus_before_final_filters[
                int(stats["samples_with_data_before_final_filters"])
            ] += 1

            # store for writing if locus passed filters
            if not sum(filters.values()):
                # Stream each retained locus BED row immediately so the BED file
                # stays synchronized with the filtered loci text output.
                scaff, pos = header.split(":")
                pos0, pos1 = (int(i) for i in pos.split("-"))
                masked_low_obs_samples = list(
                    stats.get("masked_samples_by_min_sample_observed_fraction", ())
                )
                masked_hetero_samples = list(
                    stats.get("masked_samples_by_max_sample_hetero_frequency", ())
                )
                masked_samples = [
                    *masked_low_obs_samples,
                    *masked_hetero_samples,
                ]
                visible_sample_rows = [
                    (sname, seq)
                    for sname, seq in zip(tnames, tseqs, strict=True)
                    if sname not in masked_samples
                ]
                visible_empirical_row_count = sum(
                    1 for sname, _seq in visible_sample_rows if sname != refname
                )
                out_bed.write(
                    f"{scaff}\t{pos0 - 1}\t{pos1}\t{visible_empirical_row_count}\n"
                )
                manifest_writer.writerow([oheader, header, ",".join(masked_samples)])
                if masked_low_obs_samples:
                    loci_with_samples_masked_by_min_observed_fraction += 1
                    total_masked_sample_occurrences_by_min_observed_fraction += len(
                        masked_low_obs_samples
                    )
                    for sname in masked_low_obs_samples:
                        masked_by_min_observed_fraction_counts[sname] += 1
                if masked_hetero_samples:
                    loci_with_samples_masked_by_max_hetero_frequency += 1
                    total_masked_sample_occurrences_by_max_hetero_frequency += len(masked_hetero_samples)
                    for sname in masked_hetero_samples:
                        masked_by_max_hetero_frequency_counts[sname] += 1
                if masked_samples:
                    for sname in masked_samples:
                        mask_handles[sname].write(f"{scaff}\t{pos0 - 1}\t{pos1}\n")

                # Count only empirical samples in the report summaries so the
                # synthetic reference row does not inflate occupancy metrics.
                sample_mask = _empirical_sample_row_mask(tnames, refname=refname)
                sample_rows = sample_mask & _sample_rows_with_data(tseqs)
                sample_count = int(np.sum(sample_rows))
                for row_idx, sname in enumerate(tnames):
                    if sname != refname and sample_rows[row_idx]:
                        per_sample_locus_counts[sname] += 1
                samples_per_locus[sample_count] += 1
                locus_length_counts[int(tseqs.shape[1])] += 1
                if sample_count:
                    sample_seqs = tseqs[sample_rows]
                    alignment_nonmissing_sample_bases += int(
                        np.sum((sample_seqs != 78) & (sample_seqs != 45))
                    )
                for stat in total_stats:
                    total_stats[stat] += stats[stat]

                # build locus with snpstring
                locus = []
                for sname, seq in visible_sample_rows:
                    locus.append(f"{padded[sname]}{bytes(seq).decode()}")
                snpsarr[snpsarr == 0] = 32
                snpsarr[snpsarr == 1] = 45
                snpsarr[snpsarr == 2] = 42
                snpstring = bytes(snpsarr).decode()
                locus.append(f"//{' ' * (max_len - 2)}{snpstring}|{flidx}:{header}\n")

                # store
                loci.append("\n".join(locus))
                flidx += 1
            lidx += 1

            # Flush the accumulated loci records in coarse batches so gzip still
            # sees large writes, but we do not hold the whole loci file text in
            # memory while filtering the database.
            if not flidx % 5000:
                if loci:
                    out.write("".join(loci))
                    loci = []

        # write last chunk
        if loci:
            out.write("".join(loci))
        for handle in mask_handles.values():
            handle.close()

    logger.debug("wrote final loci and BED outputs to {}", outdir)
    return {
        "nloci_before_filtering": lidx,
        "nloci_after_filtering": flidx,
        "nsites_after_filtering": int(total_stats["nsites"]),
        "filter_counts": dict(total_filters),
        "site_totals": dict(total_stats),
        "sample_locus_counts": dict(per_sample_locus_counts),
        "samples_per_locus_before_final_filters_counts": dict(
            samples_per_locus_before_final_filters
        ),
        "masked_by_min_observed_fraction_counts": dict(masked_by_min_observed_fraction_counts),
        "loci_with_samples_masked_by_min_observed_fraction": loci_with_samples_masked_by_min_observed_fraction,
        "total_masked_sample_occurrences_by_min_observed_fraction": total_masked_sample_occurrences_by_min_observed_fraction,
        "masked_by_max_hetero_frequency_counts": dict(masked_by_max_hetero_frequency_counts),
        "loci_with_samples_masked_by_max_hetero_frequency": loci_with_samples_masked_by_max_hetero_frequency,
        "total_masked_sample_occurrences_by_max_hetero_frequency": total_masked_sample_occurrences_by_max_hetero_frequency,
        "samples_per_locus_counts": dict(samples_per_locus),
        "locus_length_counts": dict(locus_length_counts),
        "alignment_nonmissing_sample_bases": alignment_nonmissing_sample_bases,
    }
