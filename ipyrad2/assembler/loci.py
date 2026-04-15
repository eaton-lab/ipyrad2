#!/usr/bin/env python

"""Consensus, locus-database, and final report helpers for assemble."""

from __future__ import annotations

import gzip
import csv
import re
import shutil
import sys
from pathlib import Path
from collections import Counter
import numpy as np
from loguru import logger
from ..utils.seqs import comp
from ..utils.jit_funcs import snp_count_numba, max_heteros_count_numba
from ..utils.parallel import run_pipeline, run_with_pool_iter

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_BED = str(BIN / "bedtools")
HETERO_CODES = np.array(list(b"RSKYWM"), dtype=np.uint8)


def get_lowdepth_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the intermediate low-depth-only mask path for one sample."""
    return tmpdir / "beds" / f"{sname}.lowdepth.mask.bed"


def get_paralog_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the sample-specific paralog exclusion BED path."""
    return tmpdir / "beds" / f"{sname}.paralog.mask.bed"


def get_indel_overlap_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the sample-specific overlapping-indel-cluster mask path."""
    return tmpdir / "beds" / f"{sname}.indel_overlap.mask.bed"


def get_sample_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the merged mask path consumed by consensus calling."""
    return tmpdir / "beds" / f"{sname}.mask.bed"


def get_consensus_hetero_mask_path(sname: str, tmpdir: Path) -> Path:
    """Return the final-output BED path for consensus heterozygosity masking."""
    return tmpdir / "beds" / f"{sname}.consensus_hetero.mask.bed"


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
    sname: str,
    reference_fasta: Path,
    tmpdir: Path,
    keep_insertions: bool,
) -> Path:
    """Write consensus sequences for one sample.

    Create FASTA for `sample_name` over the already sliced shared loci
    reference, applying variants from `vcf_gz` and masking filtered regions to N.
    """
    # The shared locus-sliced reference FASTA is written once per assemble run
    # and then reused here for every sample to avoid re-running `samtools faidx`
    # inside each consensus worker.
    vcf_gz = tmpdir / "vcfs" / "variants.resolved.vcf.gz"
    consensus_dir = tmpdir / "consensus_seqs"
    consensus_dir.mkdir(parents=True, exist_ok=True)

    # sample files
    mask_bed = get_sample_mask_path(sname, tmpdir)
    out_fasta = consensus_dir / f"{sname}.consensus.fa"

    cmd1 = [
        BIN_BCF, "consensus",
        "-f", str(reference_fasta),
        "-s", f"{sname}",         # sample to apply
        "-M", "N",                # write N for missing genotypes
        "--mask", str(mask_bed),  # mask zero/low-coverage intervals to N
        "--mask-with", "N",
        "--mark-del", "-",
        "--mark-ins", "lc" if keep_insertions else "+",
        "--regions-overlap", "1", # apply variants overlapping slice edges
        str(vcf_gz),
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


def make_lowdepth_mask(sname: str, min_sample_depth: int, tmpdir: Path):
    """Build a per-bp mask of positions inside `loci_bed` where bedGraph depth < min_depth.

    Output mask contains only the A (loci) columns and is split into minimal sub-intervals
    where coverage is below threshold (including 0-coverage gaps).
    """
    bed_dir = tmpdir / "beds"
    loci_bed = bed_dir / "loci.bed"
    ref_info = tmpdir / "REF_info.txt"
    sample_bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
    good_bed = bed_dir / f"{sname}.goodcov.bed"
    out_bed = get_lowdepth_mask_path(sname, tmpdir)

    # 1) Threshold bedGraph: keep depth >= min_depth, drop depth column for set ops
    cmd1 = [
        "awk",
        f'BEGIN{{OFS="\\t"}} $4>={min_sample_depth} {{print $1,$2,$3}}',
        str(sample_bedgraph),
    ]
    cmd2 = ["sort", "-k1,1", "-k2,2n", "-T", str(tmpdir)]
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    cmd4 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], good_bed)

    return _subtract_sorted_beds(loci_bed, good_bed, ref_info, out_bed)


def make_paralog_mask(sname: str, tmpdir: Path) -> Path:
    """Write the shared loci segments excluded only for this sample by paralog filtering."""
    bed_dir = tmpdir / "beds"
    loci_bed = bed_dir / "loci.bed"
    ref_info = tmpdir / "REF_info.txt"
    sample_good_bed = bed_dir / f"{sname}.final.good.bed"
    out_bed = get_paralog_mask_path(sname, tmpdir)

    # Samples that had no passing per-sample BED after paralog scoring get an
    # empty mask here so the downstream merge step still has a stable filepath.
    if not sample_good_bed.exists():
        out_bed.write_text("", encoding="utf-8")
        return out_bed

    return _subtract_sorted_beds(loci_bed, sample_good_bed, ref_info, out_bed)


def merge_sample_mask_beds(sname: str, tmpdir: Path) -> Path:
    """Merge all per-sample mask sources into the final consensus mask BED."""
    ref_info = tmpdir / "REF_info.txt"
    lowdepth_bed = get_lowdepth_mask_path(sname, tmpdir)
    paralog_bed = get_paralog_mask_path(sname, tmpdir)
    indel_overlap_bed = get_indel_overlap_mask_path(sname, tmpdir)
    out_bed = get_sample_mask_path(sname, tmpdir)

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
    cmd2 = ["sort", "-k1,1", "-k2,2n", "-T", str(tmpdir)]
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    cmd4 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], out_bed)
    return out_bed


def merge_final_vcf_mask_beds(sname: str, tmpdir: Path) -> Path:
    """Merge consensus-time and final-output sample masks for final VCF masking."""
    ref_info = tmpdir / "REF_info.txt"
    out_bed = get_final_vcf_mask_path(sname, tmpdir)
    existing = [
        path
        for path in (get_sample_mask_path(sname, tmpdir), get_consensus_hetero_mask_path(sname, tmpdir))
        if path.exists() and path.stat().st_size
    ]
    if not existing:
        out_bed.write_text("", encoding="utf-8")
        return out_bed
    if len(existing) == 1:
        shutil.copy2(existing[0], out_bed)
        return out_bed

    cmd1 = ["cat"] + [str(path) for path in existing]
    cmd2 = ["sort", "-k1,1", "-k2,2n", "-T", str(tmpdir)]
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
    name: str,
    snames: list[str],
    reference: Path,
    tmpdir: Path,
    masks: list[str] | None,
) -> tuple[Path, Path]:
    """Build the shared locus FASTA database and optional restriction-site mask BED."""
    # get sorted consensus fastas with reference on top
    consensus_dir = tmpdir / "consensus_seqs"
    fastas = [consensus_dir / f"{i}.consensus.fa" for i in sorted(snames)]

    # insert reference as first sample unless explicitly excluded
    reference_fa = consensus_dir / "assembly_reference_sequence.consensus.fa"
    fastas = [reference_fa] + fastas

    # get names
    snames = [i.name.rsplit(".consensus.fa")[0] for i in fastas]

    # file paths
    database = tmpdir / f"{name}.database.fa"
    bed_mask = tmpdir / f"{name}.re_mask.bed"

    # restriction site sequences to be masked
    re_masks = []
    if masks:
        for mask in masks:
            re_masks.append(re.compile(mask))
            re_masks.append(re.compile(comp(mask)[::-1]))

    with open(database, "w") as out_fa, open(bed_mask, "w") as out_bed:

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
    return database, bed_mask


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


def _trim_locus_matrix(
    seqs: np.ndarray,
    *,
    min_locus_trim_sample_coverage: int,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Trim one locus matrix to the region with sufficient sample coverage."""
    site_sample_covs = np.sum((seqs != 78) & (seqs != 45), axis=0)
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
        if sname == "assembly_reference_sequence":
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


def resolve_locus_for_output(
    header: str,
    locus_dict: dict[str, str],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
    *,
    forced_masked_samples: set[str] | None = None,
):
    """Trim one locus, apply sample-level masking, and return final summary metrics."""
    scaff, pos = header.split(":")
    rstart, rend = [int(i) for i in pos.split("-")]
    tnames = list(locus_dict.keys())
    seqs = np.array([list(bytes(seq, "utf-8")) for seq in locus_dict.values()], dtype=np.uint8)

    filters = {
        "min_length": False,
        "min_samples": False,
        "max_variant_frequency": False,
        "max_shared_hetero_frequency": False,
        "max_depth_outlier": False,
    }
    stats = {
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
        "masked_samples_by_max_sample_hetero_frequency": tuple(),
        "masked_sample_count_by_max_sample_hetero_frequency": 0,
    }

    tseqs, tsite_sample_covs, trim_left, trim_right = _trim_locus_matrix(
        seqs,
        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
    )

    if forced_masked_samples is None:
        tseqs, masked_samples, sample_props = _mask_high_hetero_samples(
            tnames,
            tseqs,
            max_sample_hetero_frequency=max_sample_hetero_frequency,
        )
    else:
        tseqs = tseqs.copy()
        masked_samples = []
        sample_props = {}
        for row_idx, sname in enumerate(tnames):
            if sname in forced_masked_samples:
                tseqs[row_idx, :] = np.uint8(ord("N"))
                masked_samples.append(sname)

    if masked_samples:
        tseqs, tsite_sample_covs, extra_left, extra_right = _trim_locus_matrix(
            tseqs,
            min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
        )
        trim_left += extra_left
        trim_right += extra_right

    stats["masked_samples_by_max_sample_hetero_frequency"] = tuple(masked_samples)
    stats["masked_sample_count_by_max_sample_hetero_frequency"] = len(masked_samples)
    if sample_props:
        stats["sample_hetero_frequencies"] = sample_props

    row_has_data = _sample_rows_with_data(tseqs)
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

    if tseqs.size:
        max_shared_h = max_heteros_count_numba(tseqs)
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


def write_assemble_stats_report(
    *,
    name: str,
    outdir: Path,
    snames: list[str],
    shared_loci_after_delimiting: int,
    shared_loci_after_paralog_filtering: int,
    loci_summary: dict[str, object],
    sample_depth_stats: dict[str, dict[str, float]],
    nsnps_written: int,
    overlap_stats: dict[str, int],
    mixed_run_summary: dict[str, int] | None = None,
) -> Path:
    """Write the final human-readable assemble summary report."""
    outpath = outdir / f"{name}.stats.txt"

    final_loci_written = int(loci_summary["nloci_after_filtering"])
    assembled_sites = int(loci_summary["nsites_after_filtering"])
    filter_counts = dict(loci_summary["filter_counts"])
    site_totals = dict(loci_summary["site_totals"])
    sample_locus_counts = dict(loci_summary["sample_locus_counts"])
    samples_per_locus_counts = Counter(loci_summary["samples_per_locus_counts"])
    locus_length_counts = Counter(loci_summary["locus_length_counts"])
    alignment_nonmissing_sample_bases = int(loci_summary["alignment_nonmissing_sample_bases"])

    summary_rows = [
        ("samples", _format_count(len(snames))),
        ("shared_loci_after_delimiting", _format_count(shared_loci_after_delimiting)),
        ("shared_loci_after_paralog_filtering", _format_count(shared_loci_after_paralog_filtering)),
        ("final_loci_written", _format_count(final_loci_written)),
        (
            "final_loci_retained_fraction_of_post_paralog",
            _format_fraction(_safe_fraction(final_loci_written, shared_loci_after_paralog_filtering)),
        ),
        (
            "final_loci_retained_fraction_of_delimited",
            _format_fraction(_safe_fraction(final_loci_written, shared_loci_after_delimiting)),
        ),
        ("assembled_sites", _format_count(assembled_sites)),
        ("final_snp_sites_written", _format_count(nsnps_written)),
        ("variable_sites", _format_count(int(site_totals["variant_sites"]))),
        (
            "phylogenetically_informative_sites",
            _format_count(int(site_totals["variant_phylo_informative_sites"])),
        ),
        (
            "alignment_matrix_occupancy_fraction",
            _format_fraction(_safe_fraction(alignment_nonmissing_sample_bases, assembled_sites * len(snames))),
        ),
        (
            "overlapping_indel_clusters_masked",
            _format_count(overlap_stats["overlapping_indel_clusters_masked"]),
        ),
        (
            "overlapping_indel_records_removed",
            _format_count(overlap_stats["overlapping_indel_records_removed"]),
        ),
        (
            "overlapping_indel_bp_masked",
            _format_count(overlap_stats["overlapping_indel_bp_masked"]),
        ),
    ]

    filtering_rows = [
        ("loci_filtered_min_length", _format_count(int(filter_counts["min_length"]))),
        ("loci_filtered_min_samples", _format_count(int(filter_counts["min_samples"]))),
        ("loci_filtered_max_variant_rate", _format_count(int(filter_counts["max_variant_frequency"]))),
        (
            "loci_filtered_max_shared_heterozygosity",
            _format_count(int(filter_counts["max_shared_hetero_frequency"])),
        ),
        ("loci_filtered_max_depth_outlier", _format_count(int(filter_counts["max_depth_outlier"]))),
    ]
    sample_masking_rows = [
        (
            "loci_with_samples_masked_by_max_hetero_frequency",
            _format_count(int(loci_summary.get("loci_with_samples_masked_by_max_hetero_frequency", 0))),
        ),
        (
            "total_masked_sample_occurrences_by_max_hetero_frequency",
            _format_count(int(loci_summary.get("total_masked_sample_occurrences_by_max_hetero_frequency", 0))),
        ),
    ]

    alignment_rows = [
        ("mean_locus_length", _format_float(_counter_mean(locus_length_counts))),
        ("median_locus_length", _format_float(_counter_median(locus_length_counts))),
        ("min_locus_length", _format_count(min(locus_length_counts) if locus_length_counts else 0)),
        ("max_locus_length", _format_count(max(locus_length_counts) if locus_length_counts else 0)),
        ("mean_samples_per_locus", _format_float(_counter_mean(samples_per_locus_counts))),
        ("median_samples_per_locus", _format_float(_counter_median(samples_per_locus_counts))),
        (
            "sites_with_sample_coverage_ge_2",
            _format_count(int(site_totals["nsites_sample_cov_greater_than_1"])),
        ),
        (
            "sites_with_sample_coverage_ge_3",
            _format_count(int(site_totals["nsites_sample_cov_greater_than_2"])),
        ),
        (
            "sites_with_sample_coverage_ge_4",
            _format_count(int(site_totals["nsites_sample_cov_greater_than_3"])),
        ),
        (
            "sites_with_sample_coverage_ge_trim_min",
            _format_count(int(site_totals["nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage"])),
        ),
    ]

    sample_headers = [
        "sample",
        "loci_in_alignment",
        "loci_in_alignment_fraction",
        "shared_loci_with_nonzero_depth",
        "shared_loci_with_nonzero_depth_fraction",
        "mean_depth_shared_loci",
        "median_depth_shared_loci",
        "mean_depth_nonzero_shared_loci",
        "median_depth_nonzero_shared_loci",
        "masked_by_max_hetero_frequency",
    ]
    sample_rows: list[list[str]] = []
    for sname in sorted(snames):
        depth_stats = sample_depth_stats[sname]
        loci_in_alignment = int(sample_locus_counts.get(sname, 0))
        nonzero_loci = int(depth_stats["shared_loci_with_nonzero_depth"])
        sample_rows.append([
            sname,
            _format_count(loci_in_alignment),
            _format_fraction(_safe_fraction(loci_in_alignment, final_loci_written)),
            _format_count(nonzero_loci),
            _format_fraction(_safe_fraction(nonzero_loci, final_loci_written)),
            _format_float(float(depth_stats["mean_depth_shared_loci"])),
            _format_float(float(depth_stats["median_depth_shared_loci"])),
            _format_float(float(depth_stats["mean_depth_nonzero_shared_loci"])),
            _format_float(float(depth_stats["median_depth_nonzero_shared_loci"])),
            _format_count(int(loci_summary.get("masked_by_max_hetero_frequency_counts", {}).get(sname, 0))),
        ])

    occupancy_headers = ["samples_with_data", "loci", "fraction_of_final_loci"]
    occupancy_rows = [
        [
            _format_count(sample_count),
            _format_count(int(samples_per_locus_counts.get(sample_count, 0))),
            _format_fraction(_safe_fraction(int(samples_per_locus_counts.get(sample_count, 0)), final_loci_written)),
        ]
        for sample_count in range(len(snames) + 1)
    ]

    lines: list[str] = []
    _append_key_value_section(lines, "Assemble Summary", summary_rows)
    if mixed_run_summary:
        mixed_rows = [
            ("rad_samples", _format_count(int(mixed_run_summary["rad_samples"]))),
            ("wgs_samples", _format_count(int(mixed_run_summary["wgs_samples"]))),
            ("loci_fail_paralog_rad", _format_count(int(mixed_run_summary["loci_fail_paralog_rad"]))),
            ("loci_fail_paralog_wgs", _format_count(int(mixed_run_summary["loci_fail_paralog_wgs"]))),
            ("loci_fail_paralog_both", _format_count(int(mixed_run_summary["loci_fail_paralog_both"]))),
            (
                "loci_pass_paralog_rad_fail_paralog_wgs",
                _format_count(int(mixed_run_summary["loci_pass_paralog_rad_fail_paralog_wgs"])),
            ),
            ("sites_supported_rad_only", _format_count(int(mixed_run_summary["sites_supported_rad_only"]))),
            ("sites_supported_wgs_only", _format_count(int(mixed_run_summary["sites_supported_wgs_only"]))),
            ("sites_supported_both", _format_count(int(mixed_run_summary["sites_supported_both"]))),
            (
                "wgs_het_genotypes_masked_by_allele_balance",
                _format_count(int(mixed_run_summary["wgs_het_genotypes_masked_by_allele_balance"])),
            ),
        ]
        _append_key_value_section(lines, "Mixed RAD/WGS Summary", mixed_rows)
    _append_key_value_section(lines, "Locus Filtering", filtering_rows)
    _append_key_value_section(lines, "Sample Masking", sample_masking_rows)
    _append_key_value_section(lines, "Alignment Summary", alignment_rows)
    _append_table_section(lines, "Sample Summary", sample_headers, sample_rows)
    _append_table_section(lines, "Locus Occupancy", occupancy_headers, occupancy_rows)
    outpath.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    logger.info("assemble stats written to {}", outpath)
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
    masked_samples = set(stats.get("masked_samples_by_max_sample_hetero_frequency", ()))

    sample_mask = np.array([sname != refname for sname in tnames], dtype=bool)
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

    # `.loci.gz` should omit samples masked by the max-sample-heterozygosity rule,
    # while fixed-axis outputs keep them as missing data for that locus.
    visible_rows = [
        (sname, seq)
        for sname, seq in zip(tnames, tseqs, strict=True)
        if sname not in masked_samples
    ]
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
        "bed_row": f"{scaff}\t{pos0 - 1}\t{pos1}\t{len(visible_rows)}\n",
        "manifest_row": (raw_header, header, ",".join(sorted(masked_samples))),
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
        )
        for key in total_filters:
            total_filters[key] += int(filters[key])
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
        "retained_loci": retained_loci,
        "chunkarr": chunkarr,
        "chunkmap": chunkmap,
    }


def write_final_outputs(
    *,
    snames: list[str],
    name: str,
    outdir: Path,
    tmpdir: Path,
    reference: Path,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
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

    database = tmpdir / f"{name}.database.fa"
    loci_file = outdir / f"{name}.loci.gz"
    manifest_path = get_retained_loci_manifest_path(name, tmpdir)
    final_loci_bed = outdir / f"{name}.bed"
    (tmpdir / "beds").mkdir(parents=True, exist_ok=True)

    real_snames = list(snames)
    refname = "assembly_reference_sequence"
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
    locus_length_counts = Counter()
    per_sample_locus_counts = {i: 0 for i in real_snames}
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    masked_by_max_hetero_frequency_counts = {i: 0 for i in real_snames}
    loci_with_samples_masked_by_max_hetero_frequency = 0
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

    nbatches = sum(1 for _ in iter_locus_batches(database, batch_size=batch_size))
    output_workers = min(max(1, cores - 1), max(1, nbatches))

    def _iter_jobs():
        for batch_idx, batch_items in iter_locus_batches(database, batch_size=batch_size):
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
    ) as out_bed, manifest_path.open("w", encoding="utf-8", newline="") as manifest_handle:
        manifest_writer = csv.writer(manifest_handle, delimiter="\t")
        manifest_writer.writerow(["raw_header", "final_header", "masked_samples"])
        mask_handles = {
            sname: get_consensus_hetero_mask_path(sname, tmpdir).open("w", encoding="utf-8")
            for sname in real_snames
        }

        def _flush_pending() -> None:
            nonlocal next_batch_idx, nloci_before_filtering, flidx
            nonlocal loci_with_samples_masked_by_max_hetero_frequency
            nonlocal total_masked_sample_occurrences_by_max_hetero_frequency
            nonlocal alignment_nonmissing_sample_bases
            while next_batch_idx in pending_results:
                batch_result = pending_results.pop(next_batch_idx)
                next_batch_idx += 1
                nloci_before_filtering += int(batch_result["nloci_before_filtering"])
                for key in total_filters:
                    total_filters[key] += int(batch_result["filter_counts"][key])
                retained_loci = batch_result["retained_loci"]
                for locus_output in retained_loci:
                    out_bed.write(str(locus_output["bed_row"]))
                    retained_scaffold_names_seen.add(str(locus_output["scaff"]))
                    manifest_writer.writerow(list(locus_output["manifest_row"]))
                    masked_samples = list(locus_output["masked_samples"])
                    if masked_samples:
                        loci_with_samples_masked_by_max_hetero_frequency += 1
                        total_masked_sample_occurrences_by_max_hetero_frequency += len(masked_samples)
                        for sname in masked_samples:
                            masked_by_max_hetero_frequency_counts[sname] += 1
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
                for batch_idx, batch_items in iter_locus_batches(database, batch_size=batch_size):
                    result = _resolve_output_batch(
                        batch_idx=batch_idx,
                        batch_items=batch_items,
                        min_locus_sample_coverage=min_locus_sample_coverage,
                        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
                        min_locus_length=min_locus_length,
                        max_locus_hetero_frequency=max_locus_hetero_frequency,
                        max_locus_variant_frequency=max_locus_variant_frequency,
                        max_sample_hetero_frequency=max_sample_hetero_frequency,
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
)-> dict[str, object]:
    """Write the final `.loci.gz` and `.bed` outputs and collect report counters."""
    # database file is in the tmpdir inside outdir
    database = tmpdir / f"{name}.database.fa"
    loci_file = outdir / f"{name}.loci.gz"
    loci_compresslevel = 1

    # Add the reference as a synthetic sample in the written database/stats
    # without mutating the caller-owned sample-name list.
    real_snames = list(snames)
    refname = "assembly_reference_sequence"
    write_snames = real_snames + [refname]

    # get name padding for loci file
    max_len = max(len(i) for i in write_snames) + 2
    padded = {n: n + (" " * (max_len - len(n))) for n in write_snames}

    # Collect final-locus summary counters in the same streaming pass that
    # writes the filtered loci and BED outputs, so the final stats report does
    # not need to reparse the full database later.
    samples_per_locus = Counter()
    locus_length_counts = Counter()
    per_sample_locus_counts = {i: 0 for i in real_snames}
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    masked_by_max_hetero_frequency_counts = {i: 0 for i in real_snames}
    loci_with_samples_masked_by_max_hetero_frequency = 0
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
            sname: get_consensus_hetero_mask_path(sname, tmpdir).open("w", encoding="utf-8")
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
            )
            result = filter_trim_locus(*args)
            header, tnames, tseqs, snpsarr, filters, stats = result

            # update total dicts
            for key in total_filters:
                total_filters[key] += int(result[4][key])

            # store for writing if locus passed filters
            if not sum(filters.values()):
                # Stream each retained locus BED row immediately so the BED file
                # stays synchronized with the filtered loci text output.
                scaff, pos = header.split(":")
                pos0, pos1 = (int(i) for i in pos.split("-"))
                masked_samples = list(stats.get("masked_samples_by_max_sample_hetero_frequency", ()))
                visible_sample_rows = [
                    (sname, seq)
                    for sname, seq in zip(tnames, tseqs, strict=True)
                    if sname not in masked_samples
                ]
                out_bed.write(f"{scaff}\t{pos0 - 1}\t{pos1}\t{len(visible_sample_rows)}\n")
                manifest_writer.writerow([oheader, header, ",".join(masked_samples)])
                if masked_samples:
                    loci_with_samples_masked_by_max_hetero_frequency += 1
                    total_masked_sample_occurrences_by_max_hetero_frequency += len(masked_samples)
                    for sname in masked_samples:
                        masked_by_max_hetero_frequency_counts[sname] += 1
                        mask_handles[sname].write(f"{scaff}\t{pos0 - 1}\t{pos1}\n")

                # Count only empirical samples in the report summaries so the
                # synthetic reference row does not inflate occupancy metrics.
                sample_mask = np.array([sname != refname for sname in tnames], dtype=bool)
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
        "masked_by_max_hetero_frequency_counts": dict(masked_by_max_hetero_frequency_counts),
        "loci_with_samples_masked_by_max_hetero_frequency": loci_with_samples_masked_by_max_hetero_frequency,
        "total_masked_sample_occurrences_by_max_hetero_frequency": total_masked_sample_occurrences_by_max_hetero_frequency,
        "samples_per_locus_counts": dict(samples_per_locus),
        "locus_length_counts": dict(locus_length_counts),
        "alignment_nonmissing_sample_bases": alignment_nonmissing_sample_bases,
    }
