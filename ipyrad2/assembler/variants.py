#!/usr/bin/env python

"""Variant calling and VCF shaping for the active `ipyrad2 assemble` path.

This module owns the current bcftools-based workflow for:
- chunked joint variant calling across delimited loci
- project-level genotype and site filtering
- SNP/indel resolution into one canonical resolved VCF
- masking ambiguous overlapping-indel clusters
- writing the final SNP-only output VCF
"""

from collections import defaultdict
import gzip
import os
import re
import sys
from pathlib import Path
from loguru import logger
from ..utils.parallel import run_pipeline, run_with_pool, stream_pipeline_lines
from ..utils.exceptions import IPyradError
from .hdf5_utils import write_retained_fai
from .loci import get_indel_overlap_mask_path
from .sort_utils import assemble_sort_with_args

BIN = Path(sys.prefix) / "bin"
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")

# ==========================================================================


def _require_existing_file(path: Path, description: str) -> Path:
    """Return an existing file path or raise a clear active-workflow error."""
    path = Path(path)
    if not path.exists():
        raise IPyradError(f"{description} not found: {path}")
    return path


def _require_nonempty_file(path: Path, description: str) -> Path:
    """Return a file path only when it exists and is not zero bytes."""
    path = _require_existing_file(path, description)
    if path.stat().st_size == 0:
        raise IPyradError(f"{description} is empty: {path}")
    return path


def _read_vcf_sample_names(vcf_gz: Path) -> list[str]:
    """Read sample names from the #CHROM header line of one gzipped VCF."""
    vcf_gz = _require_nonempty_file(vcf_gz, "VCF file")
    try:
        with gzip.open(vcf_gz, "rt", encoding="utf-8") as handle:
            for raw_line in handle:
                if raw_line.startswith("#CHROM"):
                    return raw_line.rstrip("\n").split("\t")[9:]
    except OSError as exc:
        raise IPyradError(f"Could not read VCF header from {vcf_gz}: {exc}") from exc
    raise IPyradError(f"VCF header is missing the #CHROM line: {vcf_gz}")


def _get_sorted_chunk_vcfs(vcf_dir: Path, prefix: str = "chunk") -> list[Path]:
    """Return chunk VCFs in numeric order or raise a clear error if none exist."""
    chunk_vcfs: list[tuple[int, Path]] = []
    for path in vcf_dir.glob(f"{prefix}-*.vcf.gz"):
        match = re.fullmatch(rf"{re.escape(prefix)}-(\d+)\.vcf\.gz", path.name)
        if match:
            chunk_vcfs.append((int(match.group(1)), path))
    if not chunk_vcfs:
        raise IPyradError(
            f"No chunk VCFs found in {vcf_dir}. Expected files like {prefix}-0.vcf.gz."
        )
    return [path for _, path in sorted(chunk_vcfs)]


def get_indel_overlap_clusters_bed_path(tmpdir: Path) -> Path:
    """Return the shared BED path for overlapping-indel clusters."""
    return tmpdir / "beds" / "indel.overlap_clusters.bed"


def get_variant_resolution_stats_path(tmpdir: Path) -> Path:
    """Return the tmp stats path for overlap-cluster resolution details."""
    return tmpdir / "vcfs" / "variants.resolution.stats.tsv"


def get_variant_postfilter_stats_path(tmpdir: Path) -> Path:
    """Return the tmp stats path for post-filter mixed-run variant summaries."""
    return tmpdir / "vcfs" / "variants.postfilter.stats.tsv"


def _write_variant_resolution_stats(tmpdir: Path, **stats: int) -> None:
    """Write a tiny key/value summary for the variant-resolution stage."""
    path = get_variant_resolution_stats_path(tmpdir)
    lines = [f"{key}\t{int(value)}" for key, value in stats.items()]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_variant_postfilter_stats(tmpdir: Path, **stats: int) -> None:
    """Write a tiny key/value summary for mixed-run post-filter variant steps."""
    path = get_variant_postfilter_stats_path(tmpdir)
    lines = [f"{key}\t{int(value)}" for key, value in stats.items()]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def load_variant_resolution_stats(tmpdir: Path) -> dict[str, int]:
    """Load overlap-cluster summary stats written during variant resolution."""
    path = get_variant_resolution_stats_path(tmpdir)
    stats = {
        "overlapping_indel_clusters_masked": 0,
        "overlapping_indel_records_removed": 0,
        "overlapping_indel_bp_masked": 0,
        "indel_records_inspected": 0,
    }
    if not path.exists():
        return stats
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        key, value = raw_line.split("\t", 1)
        stats[key] = int(value)
    return stats


def load_variant_postfilter_stats(tmpdir: Path) -> dict[str, int]:
    """Load mixed-run post-filter variant summary stats when present."""
    path = get_variant_postfilter_stats_path(tmpdir)
    stats = {
        "wgs_het_genotypes_masked_by_allele_balance": 0,
        "wgs_het_genotypes_examined_for_allele_balance": 0,
    }
    if not path.exists():
        return stats
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        key, value = raw_line.split("\t", 1)
        if key in stats:
            stats[key] = int(value)
    return stats


def _merge_sorted_intervals(
    intervals: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Merge sorted BED intervals using bedtools-like book-ended semantics."""
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        if not merged:
            merged.append((chrom, start, end))
            continue
        last_chrom, last_start, last_end = merged[-1]
        if chrom == last_chrom and start <= last_end:
            merged[-1] = (last_chrom, last_start, max(last_end, end))
        else:
            merged.append((chrom, start, end))
    return merged


def _parse_gt_alleles(sample_field: str) -> list[int]:
    """Return called allele indexes from one sample FORMAT field."""
    gt_field = sample_field.split(":", 1)[0]
    if gt_field in {".", "./.", ".|."}:
        return []
    alleles: list[int] = []
    for token in gt_field.replace("|", "/").split("/"):
        if token in {"", "."}:
            continue
        try:
            alleles.append(int(token))
        except ValueError:
            continue
    return alleles


def _is_heterozygous_gt(gt_field: str) -> bool:
    """Return True when a diploid genotype contains two distinct called alleles."""
    alleles = _parse_gt_alleles(gt_field)
    return len(alleles) == 2 and alleles[0] != alleles[1]


def _parse_ad_counts(ad_field: str) -> list[int] | None:
    """Return integer allele depths from one AD field or None when unavailable."""
    if ad_field in {"", "."}:
        return None
    counts: list[int] = []
    for token in ad_field.split(","):
        if token in {"", "."}:
            return None
        try:
            counts.append(int(token))
        except ValueError:
            return None
    return counts or None


def apply_wgs_het_allele_balance_mask(
    vcf_gz: Path,
    sample_names: list[str],
    *,
    low: float,
    high: float,
) -> dict[str, int]:
    """Mask WGS heterozygous genotypes when called-allele balance is too extreme."""
    stats = {
        "wgs_het_genotypes_masked_by_allele_balance": 0,
        "wgs_het_genotypes_examined_for_allele_balance": 0,
    }
    if not sample_names:
        return stats

    vcf_gz = _require_nonempty_file(vcf_gz, "Filtered project VCF")
    vcf_samples = _read_vcf_sample_names(vcf_gz)
    sample_to_index = {name: idx for idx, name in enumerate(vcf_samples)}
    unknown = sorted(set(sample_names).difference(sample_to_index))
    if unknown:
        raise IPyradError(
            f"Cannot apply WGS allele-balance masks because these samples are not present in {vcf_gz}: "
            f"{', '.join(unknown)}"
        )

    target_indexes = [sample_to_index[name] for name in sample_names]
    tmp_plain = vcf_gz.with_suffix(".allele_balance.tmp.vcf")
    tmp_gz = vcf_gz.with_suffix(".allele_balance.tmp.vcf.gz")
    tmp_index = tmp_gz.with_suffix(tmp_gz.suffix + ".csi")

    try:
        with gzip.open(vcf_gz, "rt", encoding="utf-8") as in_handle, tmp_plain.open(
            "w",
            encoding="utf-8",
        ) as out_handle:
            for raw_line in in_handle:
                if raw_line.startswith("#"):
                    out_handle.write(raw_line)
                    continue

                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    out_handle.write(raw_line)
                    continue

                format_keys = fields[8].split(":")
                try:
                    gt_idx = format_keys.index("GT")
                    ad_idx = format_keys.index("AD")
                except ValueError:
                    out_handle.write(raw_line)
                    continue

                for sample_idx in target_indexes:
                    field_idx = 9 + sample_idx
                    sample_field = fields[field_idx]
                    parts = sample_field.split(":")
                    if gt_idx >= len(parts) or ad_idx >= len(parts):
                        continue
                    gt_field = parts[gt_idx]
                    if not _is_heterozygous_gt(gt_field):
                        continue

                    alleles = _parse_gt_alleles(gt_field)
                    ad_counts = _parse_ad_counts(parts[ad_idx])
                    if ad_counts is None or max(alleles) >= len(ad_counts):
                        continue

                    called_depth = sum(ad_counts[allele] for allele in alleles)
                    if called_depth <= 0:
                        continue

                    stats["wgs_het_genotypes_examined_for_allele_balance"] += 1
                    if any(
                        (ad_counts[allele] / called_depth) < low
                        or (ad_counts[allele] / called_depth) > high
                        for allele in alleles
                    ):
                        parts[gt_idx] = "./."
                        fields[field_idx] = ":".join(parts)
                        stats["wgs_het_genotypes_masked_by_allele_balance"] += 1

                out_handle.write("\t".join(fields) + "\n")

        if stats["wgs_het_genotypes_masked_by_allele_balance"] == 0:
            return stats

        cmd1 = [
            BIN_BCF,
            "+fill-tags",
            str(tmp_plain),
            "--",
            "-t",
            "AC,AN,AF,MAF,F_MISSING",
        ]
        cmd2 = [
            BIN_BCF,
            "view",
            "-Oz",
            "-o",
            str(tmp_gz),
            "-",
        ]
        run_pipeline([cmd1, cmd2])
        run_pipeline([[BIN_BCF, "index", "-f", "-c", str(tmp_gz)]])
        os.replace(tmp_gz, vcf_gz)
        if tmp_index.exists():
            os.replace(tmp_index, vcf_gz.with_suffix(vcf_gz.suffix + ".csi"))
        return stats
    finally:
        tmp_plain.unlink(missing_ok=True)
        tmp_gz.unlink(missing_ok=True)
        tmp_index.unlink(missing_ok=True)


def summarize_variant_support_by_sample_type(
    vcf_gz: Path,
    rad_samples: list[str],
    wgs_samples: list[str],
) -> dict[str, int]:
    """Count final SNP records supported by RAD only, WGS only, both, or neither."""
    stats = {
        "sites_supported_rad_only": 0,
        "sites_supported_wgs_only": 0,
        "sites_supported_both": 0,
        "sites_supported_neither": 0,
    }
    if not rad_samples or not wgs_samples:
        return stats

    vcf_gz = _require_nonempty_file(vcf_gz, "Final SNP VCF")
    vcf_samples = _read_vcf_sample_names(vcf_gz)
    sample_to_index = {name: idx for idx, name in enumerate(vcf_samples)}
    unknown = sorted(set(rad_samples + wgs_samples).difference(sample_to_index))
    if unknown:
        raise IPyradError(
            f"Cannot summarize mixed RAD/WGS site support because these samples are not present in {vcf_gz}: "
            f"{', '.join(unknown)}"
        )

    rad_indexes = [sample_to_index[name] for name in rad_samples]
    wgs_indexes = [sample_to_index[name] for name in wgs_samples]

    with gzip.open(vcf_gz, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            sample_fields = fields[9:]
            rad_support = any(
                any(allele > 0 for allele in _parse_gt_alleles(sample_fields[idx]))
                for idx in rad_indexes
            )
            wgs_support = any(
                any(allele > 0 for allele in _parse_gt_alleles(sample_fields[idx]))
                for idx in wgs_indexes
            )
            if rad_support and wgs_support:
                stats["sites_supported_both"] += 1
            elif rad_support:
                stats["sites_supported_rad_only"] += 1
            elif wgs_support:
                stats["sites_supported_wgs_only"] += 1
            else:
                stats["sites_supported_neither"] += 1
    return stats


def _get_indel_alt_spans(
    pos1: int, ref: str, alts: list[str]
) -> dict[int, tuple[int, int]]:
    """Return exact reference spans for indel ALT alleles in one VCF record."""
    pos0 = pos1 - 1
    spans: dict[int, tuple[int, int]] = {}
    for alt_idx, alt in enumerate(alts, start=1):
        if alt == "*":
            continue
        if len(ref) > len(alt):
            spans[alt_idx] = (pos0, pos0 + len(ref))
        elif len(alt) > len(ref):
            spans[alt_idx] = (pos0, pos0 + 1)
    return spans


def _write_overlapping_indel_cluster_masks(tmpdir: Path) -> dict[str, Path]:
    """Prune overlapping-indel clusters and write affected-sample mask BEDs."""
    vcf_dir = tmpdir / "vcfs"
    resolved_vcf = _require_nonempty_file(
        vcf_dir / "variants.resolved.vcf.gz", "Resolved project VCF"
    )
    resolved_index = resolved_vcf.with_suffix(resolved_vcf.suffix + ".csi")
    pre_overlap_vcf = vcf_dir / "variants.resolved.pre_overlap_clusters.vcf.gz"
    pre_overlap_index = pre_overlap_vcf.with_suffix(pre_overlap_vcf.suffix + ".csi")
    overlap_bed = get_indel_overlap_clusters_bed_path(tmpdir)

    sample_names: list[str] = []
    clusters: list[tuple[str, int, int, list[dict[str, object]]]] = []
    current_cluster: list[dict[str, object]] = []
    cluster_chrom = ""
    cluster_end = -1
    indel_record_total = 0

    def flush_cluster() -> None:
        if len(current_cluster) > 1:
            clusters.append(
                (
                    cluster_chrom,
                    int(current_cluster[0]["start0"]),
                    cluster_end,
                    list(current_cluster),
                )
            )

    # Scan the resolved VCF once, keeping only indel records and grouping any
    # connected set of overlapping indel spans into one ambiguity cluster.
    with gzip.open(resolved_vcf, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith("##"):
                continue
            if raw_line.startswith("#CHROM"):
                sample_names = raw_line.rstrip("\n").split("\t")[9:]
                continue

            fields = raw_line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt_field = fields[:5]
            alt_spans = _get_indel_alt_spans(int(pos), ref, alt_field.split(","))
            if not alt_spans:
                continue

            indel_record_total += 1
            start0 = min(start for start, _ in alt_spans.values())
            end0 = max(end for _, end in alt_spans.values())
            record = {
                "chrom": chrom,
                "start0": start0,
                "end0": end0,
                "sample_fields": fields[9:],
                "indel_alt_indexes": set(alt_spans),
            }

            if not current_cluster:
                current_cluster = [record]
                cluster_chrom = chrom
                cluster_end = end0
                continue

            if chrom == cluster_chrom and start0 < cluster_end:
                current_cluster.append(record)
                cluster_end = max(cluster_end, end0)
                continue

            flush_cluster()
            current_cluster = [record]
            cluster_chrom = chrom
            cluster_end = end0
    flush_cluster()

    # Turn each overlap cluster into one merged BED interval, then decide which
    # samples should be masked there based on whether they actually carry an ALT
    # allele at any indel record inside that cluster.
    per_sample_intervals: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    overlap_bed_lines: list[str] = []
    overlap_record_total = 0
    overlap_bp_total = 0

    for chrom, start0, end0, records in clusters:
        overlap_bed_lines.append(f"{chrom}\t{start0}\t{end0}")
        overlap_record_total += len(records)
        overlap_bp_total += end0 - start0

        # Mark the full ambiguous cluster span only in samples that carry any
        # ALT call inside that cluster. Samples that are reference or missing
        # across the cluster keep the reference-coincident sequence.
        for sample_idx, sname in enumerate(sample_names):
            if any(
                any(
                    allele in record["indel_alt_indexes"]
                    for allele in _parse_gt_alleles(record["sample_fields"][sample_idx])
                )
                for record in records
            ):
                per_sample_intervals[sname].append((chrom, start0, end0))

    overlap_bed.write_text(
        "\n".join(overlap_bed_lines) + ("\n" if overlap_bed_lines else ""),
        encoding="utf-8",
    )

    # Write one stable per-sample overlap-mask BED, leaving unaffected samples
    # with empty files so downstream mask merging can use a uniform filepath
    # convention.
    sample_mask_paths: dict[str, Path] = {}
    for sname in sample_names:
        mask_path = get_indel_overlap_mask_path(sname, tmpdir)
        merged = _merge_sorted_intervals(per_sample_intervals[sname])
        mask_path.write_text(
            "\n".join(f"{chrom}\t{start}\t{end}" for chrom, start, end in merged)
            + ("\n" if merged else ""),
            encoding="utf-8",
        )
        sample_mask_paths[sname] = mask_path

    if not clusters:
        _write_variant_resolution_stats(
            tmpdir,
            overlapping_indel_clusters_masked=0,
            overlapping_indel_records_removed=0,
            overlapping_indel_bp_masked=0,
            indel_records_inspected=indel_record_total,
        )
        logger.info("no overlapping indel clusters found in resolved VCF")
        return sample_mask_paths

    os.replace(resolved_vcf, pre_overlap_vcf)
    if resolved_index.exists():
        os.replace(resolved_index, pre_overlap_index)
    tmp_vcf = vcf_dir / "variants.resolved.tmp.vcf.gz"

    # Drop the ambiguous overlap-cluster records from the canonical resolved VCF
    # before consensus so bcftools never has to arbitrate between them.
    cmd = [
        BIN_BCF,
        "view",
        "-T",
        f"^{str(overlap_bed)}",
        "-Oz",
        "-o",
        str(tmp_vcf),
        str(pre_overlap_vcf),
    ]
    run_pipeline([cmd])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(tmp_vcf)]])
    os.replace(tmp_vcf, resolved_vcf)
    tmp_index = tmp_vcf.with_suffix(tmp_vcf.suffix + ".csi")
    if tmp_index.exists():
        os.replace(tmp_index, resolved_index)

    logger.info(
        "masked {} overlapping-indel clusters ({} records, {} bp)",
        len(clusters),
        overlap_record_total,
        overlap_bp_total,
    )
    logger.debug(
        "overlapping indel scan inspected {} indel records across {} samples",
        indel_record_total,
        len(sample_names),
    )
    _write_variant_resolution_stats(
        tmpdir,
        overlapping_indel_clusters_masked=len(clusters),
        overlapping_indel_records_removed=overlap_record_total,
        overlapping_indel_bp_masked=overlap_bp_total,
        indel_records_inspected=indel_record_total,
    )
    return sample_mask_paths


def get_chunked_loci_beds(
    tmpdir: Path,
    nchunks: int,
    source_bed: Path | None = None,
    prefix: str = "chunk",
) -> list[Path]:
    """Split the selected loci BED into approximately even chunk BEDs."""
    loci_bed = _require_nonempty_file(
        Path(source_bed) if source_bed is not None else tmpdir / "beds" / "loci.bed",
        "Callable loci BED" if source_bed is not None else "Canonical loci BED",
    )
    lines = [line for line in loci_bed.read_text().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"No loci found in {loci_bed}.")
    nchunks = max(1, min(int(nchunks), len(lines)))
    q, r = divmod(len(lines), nchunks)

    paths = []
    i = 0
    for k in range(nchunks):
        chunk_bed = tmpdir / "beds" / f"{prefix}-{i}.bed"
        size = q + (1 if k < r else 0)
        chunk = lines[i : i + size]
        with open(chunk_bed, "w", encoding="utf-8") as out:
            out.write("\n".join(chunk))
        paths.append(chunk_bed)
        i += size
    return paths


def _count_nonempty_bed_rows(path: Path) -> int:
    """Return the number of non-empty BED rows in one file."""
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def _normalize_effective_sample_masks(
    vcf_gz: Path,
    sample_masks: dict[str, Path],
) -> tuple[list[str], dict[str, int], dict[str, Path]]:
    """Return validated non-empty sample masks keyed by VCF sample order."""
    if not sample_masks:
        return [], {}, {}

    effective_masks = {
        sname: Path(mask_bed)
        for sname, mask_bed in sample_masks.items()
        if Path(mask_bed).exists() and Path(mask_bed).stat().st_size > 0
    }
    if not effective_masks:
        return [], {}, {}

    vcf_samples = _read_vcf_sample_names(vcf_gz)
    sample_to_index = {name: idx for idx, name in enumerate(vcf_samples)}
    unknown_samples = sorted(set(effective_masks) - set(vcf_samples))
    if unknown_samples:
        raise IPyradError(
            f"Cannot apply sample masks because these samples are not present in {vcf_gz}: "
            f"{', '.join(unknown_samples)}"
        )
    return vcf_samples, sample_to_index, effective_masks


def _masked_gt_value(gt_field: str) -> str:
    """Return the missing-genotype token matching one sample GT delimiter."""
    return ".|." if "|" in gt_field else "./."


def _load_mask_entries_from_paths(
    effective_masks: dict[str, Path],
    sample_to_index: dict[str, int],
) -> dict[str, list[dict[str, object]]]:
    """Load all sample masks into per-chrom interval cursors."""
    chrom_to_samples: dict[str, dict[int, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for sname, mask_bed in sorted(effective_masks.items()):
        sample_idx = sample_to_index[sname]
        with mask_bed.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                if not raw_line.strip():
                    continue
                chrom, start, end, *_rest = raw_line.rstrip("\n").split("\t")
                chrom_to_samples[chrom][sample_idx].append((int(start), int(end)))

    return {
        chrom: [
            {"field_idx": 9 + sample_idx, "intervals": intervals, "cursor": 0}
            for sample_idx, intervals in sorted(sample_map.items())
        ]
        for chrom, sample_map in chrom_to_samples.items()
    }


def _load_mask_entries_from_shard(
    mask_shard: Path | None,
) -> dict[str, list[dict[str, object]]]:
    """Load one chunk-local sample-mask shard into per-chrom interval cursors."""
    if mask_shard is None or not mask_shard.exists() or mask_shard.stat().st_size == 0:
        return {}

    chrom_to_samples: dict[str, dict[int, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with mask_shard.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            sample_idx, chrom, start, end = raw_line.rstrip("\n").split("\t", 3)
            chrom_to_samples[chrom][int(sample_idx)].append((int(start), int(end)))

    return {
        chrom: [
            {"field_idx": 9 + sample_idx, "intervals": intervals, "cursor": 0}
            for sample_idx, intervals in sorted(sample_map.items())
        ]
        for chrom, sample_map in chrom_to_samples.items()
    }


def _write_masked_vcf_plain(
    cmds: list[list[str]],
    mask_entries: dict[str, list[dict[str, object]]],
    out_plain: Path,
) -> tuple[bool, int]:
    """Stream VCF text through one masking pass and write a plain VCF file."""
    header_lines: list[str] = []
    out_handle = None
    saw_records = False
    masked_gt_total = 0

    try:
        for raw_line in stream_pipeline_lines(cmds):
            if raw_line.startswith("#"):
                if not saw_records:
                    header_lines.append(raw_line)
                continue

            fields = raw_line.split("\t")
            if len(fields) > 9:
                chrom = fields[0]
                chrom_entries = mask_entries.get(chrom, ())
                if chrom_entries:
                    try:
                        pos0 = int(fields[1]) - 1
                    except ValueError:
                        pos0 = -1
                    if pos0 >= 0:
                        format_keys = fields[8].split(":")
                        try:
                            gt_idx = format_keys.index("GT")
                        except ValueError:
                            gt_idx = -1
                        if gt_idx >= 0:
                            for entry in chrom_entries:
                                intervals = entry["intervals"]
                                cursor = int(entry["cursor"])
                                while cursor < len(intervals) and intervals[cursor][1] <= pos0:
                                    cursor += 1
                                entry["cursor"] = cursor
                                if cursor >= len(intervals):
                                    continue
                                start, end = intervals[cursor]
                                if not (start <= pos0 < end):
                                    continue
                                field_idx = int(entry["field_idx"])
                                if field_idx >= len(fields):
                                    continue
                                sample_parts = fields[field_idx].split(":")
                                if gt_idx >= len(sample_parts):
                                    continue
                                new_gt = _masked_gt_value(sample_parts[gt_idx])
                                if sample_parts[gt_idx] != new_gt:
                                    sample_parts[gt_idx] = new_gt
                                    fields[field_idx] = ":".join(sample_parts)
                                    masked_gt_total += 1

            if out_handle is None:
                out_handle = out_plain.open("w", encoding="utf-8")
                for header_line in header_lines:
                    out_handle.write(header_line + "\n")
            saw_records = True
            out_handle.write("\t".join(fields) + "\n")
    finally:
        if out_handle is not None:
            out_handle.close()

    return saw_records, masked_gt_total


def _write_final_vcf_mask_manifest(
    tmpdir: Path,
    effective_masks: dict[str, Path],
    sample_to_index: dict[str, int],
) -> Path:
    """Write one sorted BED4 combining every sample-specific final VCF mask."""
    vcf_dir = tmpdir / "vcfs"
    unsorted_path = vcf_dir / "final-vcf-masks.unsorted.bed"
    out_path = vcf_dir / "final-vcf-masks.bed"
    ref_info = tmpdir / "REF_info.txt"

    with unsorted_path.open("w", encoding="utf-8") as out:
        for sname, mask_bed in sorted(effective_masks.items()):
            sample_idx = sample_to_index[sname]
            with mask_bed.open("r", encoding="utf-8") as handle:
                for raw_line in handle:
                    if not raw_line.strip():
                        continue
                    chrom, start, end, *_rest = raw_line.rstrip("\n").split("\t")
                    out.write(f"{chrom}\t{start}\t{end}\t{sample_idx}\n")

    try:
        cmd1 = assemble_sort_with_args(
            ["-k1,1", "-k2,2n", "-T", str(tmpdir), str(unsorted_path)]
        )
        cmd2 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
        run_pipeline([cmd1, cmd2], out_path)
    finally:
        unsorted_path.unlink(missing_ok=True)
    return out_path


def _write_final_vcf_chunk_manifest(tmpdir: Path, chunk_beds: list[Path]) -> Path:
    """Write one BED4 describing every final-VCF chunk interval."""
    path = tmpdir / "vcfs" / "final-vcf-chunks.bed"
    with path.open("w", encoding="utf-8") as out:
        for chunk_idx, chunk_bed in enumerate(chunk_beds):
            with chunk_bed.open("r", encoding="utf-8") as handle:
                for raw_line in handle:
                    if not raw_line.strip():
                        continue
                    chrom, start, end, *_rest = raw_line.rstrip("\n").split("\t")
                    out.write(f"{chrom}\t{start}\t{end}\t{chunk_idx}\n")
    return path


def _write_final_vcf_mask_shards(
    tmpdir: Path,
    effective_masks: dict[str, Path],
    sample_to_index: dict[str, int],
    chunk_beds: list[Path],
) -> dict[int, Path]:
    """Split merged sample-mask BEDs into chunk-local mask shard files."""
    if not effective_masks or not chunk_beds:
        return {}

    ref_info = tmpdir / "REF_info.txt"
    masks_path = _write_final_vcf_mask_manifest(tmpdir, effective_masks, sample_to_index)
    chunks_path = _write_final_vcf_chunk_manifest(tmpdir, chunk_beds)
    shard_paths: dict[int, Path] = {}
    handles: dict[int, object] = {}

    cmd = [
        BIN_BED,
        "intersect",
        "-a",
        str(masks_path),
        "-b",
        str(chunks_path),
        "-wa",
        "-wb",
        "-sorted",
        "-g",
        str(ref_info),
    ]

    try:
        for raw_line in stream_pipeline_lines([cmd]):
            if not raw_line:
                continue
            parts = raw_line.split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            sample_idx = parts[3]
            chunk_start = int(parts[5])
            chunk_end = int(parts[6])
            chunk_idx = int(parts[7])
            overlap_start = max(start, chunk_start)
            overlap_end = min(end, chunk_end)
            if overlap_start >= overlap_end:
                continue
            shard_path = shard_paths.setdefault(
                chunk_idx,
                tmpdir / "vcfs" / f"final-vcf-mask-shard-{chunk_idx}.tsv",
            )
            handle = handles.get(chunk_idx)
            if handle is None:
                handle = shard_path.open("w", encoding="utf-8")
                handles[chunk_idx] = handle
            handle.write(f"{sample_idx}\t{chrom}\t{overlap_start}\t{overlap_end}\n")
    finally:
        for handle in handles.values():
            handle.close()
        masks_path.unlink(missing_ok=True)
        chunks_path.unlink(missing_ok=True)

    return shard_paths


def _mask_final_vcf_chunk(
    vcf_gz: Path,
    chunk_bed: Path,
    out_vcf_path: Path,
    mask_shard: Path | None,
) -> Path | None:
    """Write one plain masked final-VCF chunk or return None when it has no records."""
    cmds = [
        [
            BIN_BCF,
            "view",
            "-T",
            str(chunk_bed),
            "-f",
            "PASS",
            "-V",
            "indels",
            str(vcf_gz),
        ]
    ]
    mask_entries = _load_mask_entries_from_shard(mask_shard)

    saw_records, _masked_gt_total = _write_masked_vcf_plain(
        cmds,
        mask_entries,
        out_vcf_path,
    )
    if not saw_records:
        out_vcf_path.unlink(missing_ok=True)
        return None
    return out_vcf_path


def _concatenate_plain_vcf_chunks(chunk_vcfs: list[Path], out_plain: Path) -> Path:
    """Concatenate plain chunk VCFs, writing headers once and rows in order."""
    header_written = False
    with out_plain.open("w", encoding="utf-8") as out:
        for chunk_vcf in chunk_vcfs:
            with chunk_vcf.open("r", encoding="utf-8") as handle:
                for raw_line in handle:
                    if raw_line.startswith("#"):
                        if not header_written:
                            out.write(raw_line)
                        continue
                    out.write(raw_line)
            header_written = True
    return out_plain


def apply_sample_region_masks_to_resolved_vcf(
    tmpdir: Path,
    sample_masks: dict[str, Path],
    vcf_gz: Path | None = None,
) -> Path:
    """Set per-sample genotypes to missing inside sample-specific BED masks.

    By default this edits the canonical resolved VCF in tmpdir, but callers can
    also point it at a smaller downstream VCF when they want the same per-sample
    genotype-missing semantics without repeatedly rewriting the larger project
    VCF earlier in the assemble workflow.
    """
    current_vcf = (
        Path(vcf_gz)
        if vcf_gz is not None
        else tmpdir / "vcfs" / "variants.resolved.vcf.gz"
    )
    current_index = current_vcf.with_suffix(current_vcf.suffix + ".csi")

    if not sample_masks:
        return current_vcf

    _require_nonempty_file(current_vcf, "Target VCF for sample masking")
    _vcf_samples, sample_to_index, effective_masks = _normalize_effective_sample_masks(
        current_vcf,
        sample_masks,
    )
    if not effective_masks:
        return current_vcf

    mask_entries = _load_mask_entries_from_paths(effective_masks, sample_to_index)
    tmp_plain = current_vcf.with_suffix(".sample_masks.tmp.vcf")
    tmp_gz = current_vcf.with_suffix(".sample_masks.tmp.vcf.gz")
    tmp_index = tmp_gz.with_suffix(tmp_gz.suffix + ".csi")

    try:
        saw_records, masked_gt_total = _write_masked_vcf_plain(
            [[BIN_BCF, "view", str(current_vcf)]],
            mask_entries,
            tmp_plain,
        )
        if not saw_records or masked_gt_total == 0:
            return current_vcf

        cmd = [
            BIN_BCF,
            "view",
            "-Oz",
            "-o",
            str(tmp_gz),
            str(tmp_plain),
        ]
        run_pipeline([cmd])
        run_pipeline([[BIN_BCF, "index", "-f", "-c", str(tmp_gz)]])
        os.replace(tmp_gz, current_vcf)
        if tmp_index.exists():
            os.replace(tmp_index, current_index)
        return current_vcf
    finally:
        tmp_plain.unlink(missing_ok=True)
        tmp_gz.unlink(missing_ok=True)
        tmp_index.unlink(missing_ok=True)
    return current_vcf


def get_group_called_variants_in_vcf_chunks(
    tmpdir: Path,
    reference: Path,
    bam_files: list[Path],
    locus_chunk: Path,
    min_map_q: int,
    min_base_q: int,
    threads: int,
    group_samples_file: Path | None = None,
):
    """Joint-call one loci BED chunk across all filtered analysis BAMs."""
    _require_existing_file(reference, "Reference FASTA")
    _require_nonempty_file(locus_chunk, "Variant-calling locus chunk BED")
    if not bam_files:
        raise IPyradError("No BAM files were provided for joint variant calling.")
    if group_samples_file is not None:
        group_samples_file = _require_nonempty_file(
            group_samples_file,
            "Grouped-calling populations file",
        )

    out_vcf_gz = tmpdir / "vcfs" / locus_chunk.with_suffix(".vcf.gz").name

    threads_mpileup = max(1, threads)

    # Compute genotype likelihoods only inside this chunk BED, then call and
    # keep just SNP and indel records in one per-chunk compressed VCF.
    cmd1 = [
        BIN_BCF,
        "mpileup",
        "-f",
        str(reference),
        "-q",
        str(min_map_q),
        "-Q",
        str(min_base_q),
        "-d",
        str(5000),
        "-a",
        "FMT/DP,FMT/AD",
        "-R",
        str(locus_chunk),
        "--threads",
        str(threads_mpileup),
        "-Ou",
    ] + [str(i) for i in bam_files]

    cmd2 = [
        BIN_BCF,
        "call",
        "-m",
        "-a",
        "GQ",
        "-G",
        str(group_samples_file) if group_samples_file is not None else "-",
        "--ploidy",
        "2",
        "-Ou",
        "--threads",
        "1",
    ]
    cmd3 = [
        BIN_BCF,
        "view",
        "-v",
        "snps,indels",
        "--threads",
        "1",
        "-Oz",
        "-o",
        str(out_vcf_gz),
    ]
    run_pipeline([cmd1, cmd2, cmd3])
    return out_vcf_gz


def get_concat_chunk_vcfs(tmpdir: Path, threads: int):
    """Concatenate numerically ordered chunk VCFs into the raw project VCF."""
    vcf_dir = tmpdir / "vcfs"
    out_vcf_gz = vcf_dir / "loci.raw.vcf.gz"
    sorted_vcfs = _get_sorted_chunk_vcfs(vcf_dir)
    for chunk_vcf in sorted_vcfs:
        _require_existing_file(chunk_vcf, "Chunk VCF")

    cmd = [
        BIN_BCF,
        "concat",
        "--threads",
        str(threads),
        "-Oz",
        "-o",
        str(out_vcf_gz),
        "-W",
    ] + [str(i) for i in sorted_vcfs]
    run_pipeline([cmd])
    return out_vcf_gz


def get_filtered_vcf(
    tmpdir: Path, min_read_depth: int, min_geno_q: int, min_site_q: int, threads: int
) -> Path:
    """Mask low-confidence genotypes and annotate the filtered project VCF."""
    in_vcf_gz = _require_nonempty_file(
        tmpdir / "vcfs" / "loci.raw.vcf.gz", "Raw project VCF"
    )
    out_vcf_gz = tmpdir / "vcfs" / "loci.filtered.vcf.gz"
    out_vcf_tmp = out_vcf_gz.with_suffix(out_vcf_gz.suffix + ".tmp")

    dp_min: int = min_read_depth
    gq_min: int = min_geno_q
    qual_min: int = min_site_q
    threads = max(1, int(threads / 2))

    expr_gt_mask = f"FMT/DP<{dp_min} | FMT/GQ<{gq_min}"
    cmd1 = [
        BIN_BCF,
        "+setGT",
        str(in_vcf_gz),
        "--",
        "-t",
        "q",
        "-n",
        ".",
        "-i",
        expr_gt_mask,
    ]

    expr_site_mask = f"QUAL<{qual_min}"
    cmd2 = [
        BIN_BCF,
        "filter",
        "-S",
        ".",
        "-s",
        "lowQual",
        "-e",
        expr_site_mask,
        "--threads",
        str(threads),
        "-Ou",
        "-",
    ]

    cmd3 = [
        BIN_BCF,
        "+fill-tags",
        "-",
        "--",
        "-t",
        "AC,AN,AF,MAF,F_MISSING",
    ]

    remove_tags = "FORMAT/PL,FORMAT/GQ,INFO/RPBZ,INFO/SCBZ,INFO/MQBZ,INFO/BQBZ,INFO/MQSBZ,INFO/DP4,INFO/VDB,INFO/MQ0F,INFO/SGB"
    cmd4 = [
        BIN_BCF,
        "annotate",
        "-x",
        remove_tags,
        "-Oz",
        "-o",
        str(out_vcf_tmp),
        "--threads",
        str(threads),
        "-",
    ]

    run_pipeline([cmd1, cmd2, cmd3, cmd4])

    os.replace(out_vcf_tmp, out_vcf_gz)
    cmd = [BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)]
    run_pipeline([cmd])
    return out_vcf_gz


def get_vcf_with_indels_resolved(tmpdir: Path, reference: Path, threads: int) -> Path:
    """Build the canonical resolved VCF used by consensus and final outputs."""
    in_vcf_gz = _require_nonempty_file(
        tmpdir / "vcfs" / "loci.filtered.vcf.gz", "Filtered project VCF"
    )
    _require_existing_file(reference, "Reference FASTA")
    out_vcf_gz = tmpdir / "vcfs" / "variants.resolved.vcf.gz"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    indel_beds = bed_dir / "indel.regions.bed"

    # Normalize into one-record-per-ALT form so SNP/indel separation and
    # span-based masking are based on explicit, stable records.
    cmd1 = [
        BIN_BCF,
        "norm",
        "-f",
        str(reference),
        "-m",
        "-both",
        "--threads",
        str(threads),
        "-W",
        "-Oz",
        "-o",
        str(vcf_dir / "norm.vcf.gz"),
        str(in_vcf_gz),
    ]
    run_pipeline([cmd1])

    # Split the normalized records so simple SNP cleanup can happen without
    # losing the original indel records.
    cmd1 = [
        BIN_BCF,
        "view",
        "-v",
        "snps",
        "-Oz",
        "-o",
        str(vcf_dir / "snps.vcf.gz"),
        "--threads",
        str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    run_pipeline([cmd1])

    cmd1 = [
        BIN_BCF,
        "view",
        "-v",
        "indels",
        "-Oz",
        "-o",
        str(vcf_dir / "indels.vcf.gz"),
        "--threads",
        str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    run_pipeline([cmd1])

    # Convert every indel into the exact affected reference span. Deletions
    # cover the deleted REF bases; insertions cover only their 1 bp anchor.
    awk_prog = (
        r'BEGIN{OFS="\t"}'
        r'{chrom=$1; pos0=$2; ref=$4; n=split($5,alts,",");'
        r" for(i=1;i<=n;i++){alt=alts[i];"
        r"  if(length(ref)>length(alt)){print chrom, pos0, pos0+length(ref);} "
        r"  else if(length(alt)>length(ref)){print chrom, pos0, pos0+1;} "
        r" }}"
    )
    cmd1 = [
        BIN_BCF,
        "query",
        "-f",
        r"%CHROM\t%POS0\t%POS\t%REF\t%ALT\n",
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = ["awk", awk_prog]
    cmd3 = assemble_sort_with_args(["-k1,1", "-k2,2n", "-T", str(vcf_dir)])
    cmd4 = [BIN_BED, "merge", "-i", "-"]
    run_pipeline([cmd1, cmd2, cmd3, cmd4], indel_beds)

    if indel_beds.stat().st_size == 0:
        logger.info("no indels found after filtering; resolved VCF remains SNP-only")
        cmd1 = [
            BIN_BCF,
            "norm",
            "-m",
            "+snps",
            "--threads",
            str(threads),
            str(vcf_dir / "snps.vcf.gz"),
        ]
        cmd2 = [
            BIN_BCF,
            "sort",
            "-Oz",
            "-o",
            str(out_vcf_gz),
            "-W",
        ]
        run_pipeline([cmd1, cmd2])
        run_pipeline([[BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)]])
        _write_overlapping_indel_cluster_masks(tmpdir)
        return out_vcf_gz

    # Drop only SNPs whose reference spans overlap indel-affected intervals.
    cmd1 = [
        BIN_BCF,
        "view",
        "-T",
        f"^{str(indel_beds)}",
        "-Oz",
        "-o",
        str(vcf_dir / "snps.clean.vcf.gz"),
        "--threads",
        str(threads),
        "-W",
        str(vcf_dir / "snps.vcf.gz"),
    ]
    run_pipeline([cmd1])

    # Recombine cleaned SNPs with indels, then collapse same-position biallelic
    # records back to multiallelic form for the canonical resolved VCF.
    cmd1 = [
        BIN_BCF,
        "concat",
        "-a",
        "-Oz",
        "-o",
        str(vcf_dir / "combined.vcf.gz"),
        "--threads",
        str(threads),
        str(vcf_dir / "snps.clean.vcf.gz"),
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF,
        "sort",
        "-Oz",
        "-o",
        str(vcf_dir / "combined.sorted.vcf.gz"),
        "-T",
        str(vcf_dir),
        "-W",
        str(vcf_dir / "combined.vcf.gz"),
    ]
    run_pipeline([cmd1])
    run_pipeline([cmd2])

    cmd1 = [
        BIN_BCF,
        "norm",
        "-m",
        "+both",
        "--threads",
        str(threads),
        str(vcf_dir / "combined.sorted.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF,
        "sort",
        "-Oz",
        "-o",
        str(out_vcf_gz),
        "-W",
    ]
    run_pipeline([cmd1, cmd2])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)]])

    # Remove only genuinely ambiguous overlap clusters. Simple indels stay in
    # the resolved VCF, while overlapping clusters are pruned and propagated
    # into per-sample BED masks for downstream consensus masking.
    _write_overlapping_indel_cluster_masks(tmpdir)
    return out_vcf_gz


def compact_resolved_vcf_to_final_loci_contigs(
    tmpdir: Path, reference: Path, loci_bed: Path
) -> Path:
    """Trim resolved-VCF contig headers down to the final retained locus scaffolds."""
    loci_bed = _require_nonempty_file(loci_bed, "Final loci BED")
    resolved_vcf = _require_nonempty_file(
        tmpdir / "vcfs" / "variants.resolved.vcf.gz", "Resolved project VCF"
    )
    current_index = resolved_vcf.with_suffix(resolved_vcf.suffix + ".csi")
    vcf_dir = tmpdir / "vcfs"
    subset_fai = write_retained_fai(
        reference, loci_bed, vcf_dir / "variants.resolved.contigs.fai"
    )
    tmp_vcf = vcf_dir / "variants.resolved.reheadered.vcf.gz"
    tmp_index = tmp_vcf.with_suffix(tmp_vcf.suffix + ".csi")

    cmd = [
        BIN_BCF,
        "reheader",
        "-f",
        str(subset_fai),
        "-o",
        str(tmp_vcf),
        str(resolved_vcf),
    ]
    run_pipeline([cmd])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(tmp_vcf)]])
    os.replace(tmp_vcf, resolved_vcf)
    if tmp_index.exists():
        os.replace(tmp_index, current_index)
    return resolved_vcf


def _write_unmasked_final_vcf(
    *,
    loci_bed: Path,
    in_vcf_gz: Path,
    out_vcf_gz: Path,
    threads: int,
) -> Path:
    """Write the final SNP-only VCF without any extra sample-region masking."""
    cmd = [
        BIN_BCF,
        "view",
        "-T",
        str(loci_bed),
        "-Oz",
        "-o",
        str(out_vcf_gz),
        "--threads",
        str(threads),
        "-f",
        "PASS",
        "-V",
        "indels",
        str(in_vcf_gz),
    ]
    run_pipeline([cmd])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)]])
    return out_vcf_gz


def write_vcf(
    name: str,
    outdir: Path,
    tmpdir: Path,
    threads: int,
    *,
    sample_masks: dict[str, Path] | None = None,
    cores: int | None = None,
    log_level: str = "INFO",
) -> Path:
    """Write the final SNP-only project VCF trimmed to the final loci BED."""
    loci_bed = _require_nonempty_file(outdir / f"{name}.bed", "Final loci BED")
    out_vcf_gz = outdir / f"{name}.vcf.gz"
    in_vcf_gz = _require_nonempty_file(
        tmpdir / "vcfs" / "variants.resolved.vcf.gz", "Resolved project VCF"
    )

    _vcf_samples, sample_to_index, effective_masks = _normalize_effective_sample_masks(
        in_vcf_gz,
        sample_masks or {},
    )
    if not effective_masks:
        return _write_unmasked_final_vcf(
            loci_bed=loci_bed,
            in_vcf_gz=in_vcf_gz,
            out_vcf_gz=out_vcf_gz,
            threads=threads,
        )

    nloci = _count_nonempty_bed_rows(loci_bed)
    chunk_workers = max(1, min(int(cores or 1), nloci))
    nchunks = max(1, min(nloci, chunk_workers * 4))
    chunk_beds = get_chunked_loci_beds(
        tmpdir,
        nchunks=nchunks,
        source_bed=loci_bed,
        prefix="final-vcf-chunk",
    )
    shard_paths = _write_final_vcf_mask_shards(
        tmpdir,
        effective_masks,
        sample_to_index,
        chunk_beds,
    )
    vcf_dir = tmpdir / "vcfs"
    jobs = {
        chunk_idx: (
            _mask_final_vcf_chunk,
            dict(
                vcf_gz=in_vcf_gz,
                chunk_bed=chunk_bed,
                out_vcf_path=vcf_dir / f"final-vcf-chunk-{chunk_idx}.vcf",
                mask_shard=shard_paths.get(chunk_idx),
            ),
        )
        for chunk_idx, chunk_bed in enumerate(chunk_beds)
    }

    logger.info("masking final VCF chunks")
    if len(jobs) == 1 or chunk_workers == 1:
        chunk_results = {
            chunk_idx: func(**kwargs)
            for chunk_idx, (func, kwargs) in jobs.items()
        }
    else:
        chunk_results = run_with_pool(
            jobs,
            log_level,
            max_workers=chunk_workers,
            msg="Masking final VCF chunks",
        )

    chunk_vcfs = [
        chunk_results[chunk_idx]
        for chunk_idx in sorted(chunk_results)
        if chunk_results[chunk_idx] is not None
    ]
    if not chunk_vcfs:
        return _write_unmasked_final_vcf(
            loci_bed=loci_bed,
            in_vcf_gz=in_vcf_gz,
            out_vcf_gz=out_vcf_gz,
            threads=threads,
        )

    logger.info("concatenating masked final VCF chunks")
    tmp_plain = vcf_dir / "final-vcf.masked.concat.vcf"
    try:
        _concatenate_plain_vcf_chunks(chunk_vcfs, tmp_plain)
        cmd = [
            BIN_BCF,
            "view",
            "-Oz",
            "-o",
            str(out_vcf_gz),
            "--threads",
            str(max(1, threads)),
            str(tmp_plain),
        ]
        run_pipeline([cmd])
        run_pipeline([[BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)]])
    finally:
        tmp_plain.unlink(missing_ok=True)

    for chunk_vcf in chunk_vcfs:
        Path(chunk_vcf).unlink(missing_ok=True)
    for chunk_bed in chunk_beds:
        Path(chunk_bed).unlink(missing_ok=True)
    for shard_path in shard_paths.values():
        Path(shard_path).unlink(missing_ok=True)
    return out_vcf_gz
