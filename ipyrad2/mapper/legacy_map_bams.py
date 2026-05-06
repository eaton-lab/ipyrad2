#!/usr/bin/env python

"""Migrate legacy ipyrad2 map BAMs and stats to current naming conventions."""

from __future__ import annotations

from dataclasses import dataclass
import json
import sys
import tempfile
from pathlib import Path
from typing import Iterable

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.names import normalize_workflow_sample_name
from ..utils.parallel import run_pipeline
from .map_stats import render_map_stats_payload_report


BIN = Path(sys.prefix) / "bin"
BIN_SAMTOOLS = str(BIN / "samtools")
LEGACY_BAM_SUFFIX = ".trimmed.filtered.bam"
CURRENT_BAM_SUFFIX = ".trimmed.sorted.bam"
MAP_STATS_JSON_GLOB = "ipyrad_map_stats_*.json"
MAP_STATS_TXT_GLOB = "ipyrad_map_stats_*.txt"


@dataclass(frozen=True)
class LegacyBamMigration:
    """Planned migration for one legacy map BAM."""

    source_bam: Path
    output_bam: Path
    output_index: Path
    legacy_sample: str
    canonical_sample: str


@dataclass(frozen=True)
class MapStatsMigration:
    """Planned rewrite for one mapper stats sidecar set."""

    source_json: Path
    output_json: Path
    output_txt: Path


@dataclass(frozen=True)
class LegacyMapMigrationPlan:
    """Full migration plan for one input directory."""

    bam_migrations: list[LegacyBamMigration]
    skipped_bams: list[Path]
    stats_migrations: list[MapStatsMigration]
    orphan_stats_txt: list[Path]


def _resolve_dir(path: str | Path) -> Path:
    """Return an absolute directory path after user expansion."""
    return Path(path).expanduser().resolve()


def is_legacy_map_bam(path: Path) -> bool:
    """Return True when a BAM name matches the legacy/current map patterns."""
    return path.name.endswith((LEGACY_BAM_SUFFIX, CURRENT_BAM_SUFFIX))


def migrated_bam_name(name: str) -> str:
    """Return the current BAM basename for one legacy/current BAM filename."""
    if name.endswith(LEGACY_BAM_SUFFIX):
        return name[:-len(LEGACY_BAM_SUFFIX)] + CURRENT_BAM_SUFFIX
    if name.endswith(CURRENT_BAM_SUFFIX):
        return name
    raise IPyradError(
        "Expected a legacy/current ipyrad2 map BAM ending in "
        f"'{LEGACY_BAM_SUFFIX}' or '{CURRENT_BAM_SUFFIX}', not '{name}'."
    )


def _parse_samtools_samples_output(text: str, bam_path: Path) -> str:
    """Return the single effective sample name reported by `samtools samples`."""
    seen: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split()
        if not fields:
            continue
        sample_name = fields[0].strip()
        if not sample_name or sample_name == ".":
            continue
        if sample_name not in seen:
            seen.append(sample_name)

    if not seen:
        raise IPyradError(
            f"Could not determine a sample name from BAM header: {bam_path}"
        )
    if len(seen) != 1:
        joined = ", ".join(seen)
        raise IPyradError(
            "Legacy BAM migration expects one effective sample per BAM, but "
            f"`samtools samples` reported multiple names for {bam_path}: {joined}"
        )
    return seen[0]


def read_bam_sample_name(bam_path: Path) -> str:
    """Read the single effective sample name recorded in one BAM header."""
    cmd = [BIN_SAMTOOLS, "samples", str(bam_path)]
    _, out, _ = run_pipeline([cmd])
    text = out.decode() if isinstance(out, bytes) else str(out)
    return _parse_samtools_samples_output(text, bam_path)


def discover_legacy_map_bams(indir: Path) -> tuple[list[Path], list[Path]]:
    """Return recognized map BAMs and unrelated BAMs in one directory."""
    bam_paths = sorted(
        path for path in indir.iterdir()
        if path.is_file() and path.suffix == ".bam"
    )
    selected = [path for path in bam_paths if is_legacy_map_bam(path)]
    skipped = [path for path in bam_paths if path not in selected]
    return selected, skipped


def _plan_bam_migrations(
    indir: Path,
    outdir: Path,
    bam_paths: Iterable[Path],
) -> list[LegacyBamMigration]:
    """Plan BAM rewrites and validate canonical-name collisions."""
    migrations: list[LegacyBamMigration] = []
    canonical_sources: dict[str, list[Path]] = {}

    for bam_path in bam_paths:
        legacy_sample = read_bam_sample_name(bam_path)
        canonical_sample = normalize_workflow_sample_name(legacy_sample)
        output_bam = outdir / migrated_bam_name(bam_path.name)
        migrations.append(
            LegacyBamMigration(
                source_bam=bam_path,
                output_bam=output_bam,
                output_index=output_bam.with_suffix(output_bam.suffix + ".csi"),
                legacy_sample=legacy_sample,
                canonical_sample=canonical_sample,
            )
        )
        canonical_sources.setdefault(canonical_sample, []).append(bam_path)

    collisions = {
        canonical_sample: sorted(str(path) for path in bam_paths_for_sample)
        for canonical_sample, bam_paths_for_sample in canonical_sources.items()
        if len(bam_paths_for_sample) > 1
    }
    if collisions:
        detail = "; ".join(
            f"{canonical_sample} <- {', '.join(paths)}"
            for canonical_sample, paths in sorted(collisions.items())
        )
        raise IPyradError(
            "Legacy BAM sample names collide after stripping one terminal "
            f"'.trimmed' suffix: {detail}"
        )

    output_sources: dict[Path, list[Path]] = {}
    for migration in migrations:
        output_sources.setdefault(migration.output_bam, []).append(migration.source_bam)
    duplicate_outputs = {
        output_bam: sorted(str(path) for path in source_paths)
        for output_bam, source_paths in output_sources.items()
        if len(source_paths) > 1
    }
    if duplicate_outputs:
        detail = "; ".join(
            f"{output_bam} <- {', '.join(paths)}"
            for output_bam, paths in sorted(duplicate_outputs.items(), key=lambda item: str(item[0]))
        )
        raise IPyradError(
            "Legacy BAM migration would overwrite multiple inputs onto the same "
            f"output BAM path: {detail}"
        )

    return migrations


def _discover_stats_migrations(indir: Path, outdir: Path) -> tuple[list[MapStatsMigration], list[Path]]:
    """Return matching map-stats JSON/TXT rewrites and orphan TXT files."""
    json_paths = sorted(path for path in indir.glob(MAP_STATS_JSON_GLOB) if path.is_file())
    stats_migrations = [
        MapStatsMigration(
            source_json=path,
            output_json=outdir / path.name,
            output_txt=outdir / (path.stem + ".txt"),
        )
        for path in json_paths
    ]
    json_names = {path.name for path in json_paths}
    orphan_txt = [
        path
        for path in sorted(indir.glob(MAP_STATS_TXT_GLOB))
        if path.is_file() and path.with_suffix(".json").name not in json_names
    ]
    return stats_migrations, orphan_txt


def plan_legacy_map_migration(indir: str | Path, outdir: str | Path) -> LegacyMapMigrationPlan:
    """Plan a legacy map migration without writing any files."""
    resolved_indir = _resolve_dir(indir)
    resolved_outdir = _resolve_dir(outdir)
    if resolved_indir == resolved_outdir:
        raise IPyradError("--indir and --outdir must be different directories.")
    if not resolved_indir.is_dir():
        raise IPyradError(f"--indir does not exist or is not a directory: {resolved_indir}")

    bam_paths, skipped_bams = discover_legacy_map_bams(resolved_indir)
    bam_migrations = _plan_bam_migrations(resolved_indir, resolved_outdir, bam_paths)
    stats_migrations, orphan_stats_txt = _discover_stats_migrations(resolved_indir, resolved_outdir)

    if not bam_migrations and not stats_migrations:
        raise IPyradError(
            "No legacy/current ipyrad2 map BAMs or map stats JSON files were found in "
            f"{resolved_indir}"
        )

    return LegacyMapMigrationPlan(
        bam_migrations=bam_migrations,
        skipped_bams=skipped_bams,
        stats_migrations=stats_migrations,
        orphan_stats_txt=orphan_stats_txt,
    )


def _ensure_targets_writable(paths: Iterable[Path], force: bool) -> None:
    """Raise unless the requested output paths can be written."""
    existing = sorted(path for path in paths if path.exists())
    if existing and not force:
        shown = ", ".join(str(path) for path in existing[:5])
        if len(existing) > 5:
            shown += f", ... ({len(existing) - 5} more)"
        raise IPyradError(
            "Refusing to overwrite existing migration outputs without --force: "
            f"{shown}"
        )


def _build_canonical_rg_line(sample_name: str) -> str:
    """Return the single canonical RG header line for one migrated BAM."""
    return f"@RG\tID:{sample_name}\tSM:{sample_name}"


def _rewrite_header_with_single_rg(header_text: str, sample_name: str) -> str:
    """Replace all @RG lines in one SAM header with a single canonical RG line."""
    lines = [line for line in header_text.splitlines() if line]
    kept_lines = [line for line in lines if not line.startswith("@RG\t")]
    rg_line = _build_canonical_rg_line(sample_name)

    insert_at = 0
    for idx, line in enumerate(kept_lines):
        if line.startswith(("@HD\t", "@SQ\t")):
            insert_at = idx + 1
    kept_lines.insert(insert_at, rg_line)
    return "\n".join(kept_lines) + "\n"


def _rewrite_bam_rg_tags(migration: LegacyBamMigration, threads: int) -> None:
    """Rewrite one BAM RG header/read tags to the canonical sample name."""
    with tempfile.TemporaryDirectory(prefix="legacy-map-bam-") as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        rg_tmp = tmpdir / "rg_updated.bam"
        header_tmp = tmpdir / "header.sam"

        add_rg_cmd = [
            BIN_SAMTOOLS,
            "addreplacerg",
            "-m",
            "overwrite_all",
            "-r",
            _build_canonical_rg_line(migration.canonical_sample),
            "-o",
            str(rg_tmp),
            "-@",
            str(max(0, threads - 1)),
            str(migration.source_bam),
        ]
        run_pipeline([add_rg_cmd])

        header_cmd = [BIN_SAMTOOLS, "view", "-H", str(rg_tmp)]
        _, out, _ = run_pipeline([header_cmd])
        header_text = out.decode() if isinstance(out, bytes) else str(out)
        header_tmp.write_text(
            _rewrite_header_with_single_rg(header_text, migration.canonical_sample),
            encoding="utf-8",
        )

        reheader_cmd = [
            BIN_SAMTOOLS,
            "reheader",
            "-P",
            str(header_tmp),
            str(rg_tmp),
        ]
        run_pipeline([reheader_cmd], outfile=migration.output_bam)

    index_cmd = [
        BIN_SAMTOOLS,
        "index",
        "-c",
        "-@",
        str(max(0, threads - 1)),
        str(migration.output_bam),
    ]
    run_pipeline([index_cmd])


def _normalize_stats_rows(rows: list[dict]) -> list[dict]:
    """Return migrated stats rows with canonical sample labels."""
    normalized_rows: list[dict] = []
    for row in rows:
        new_row = dict(row)
        sample_name = new_row.get("sample")
        if isinstance(sample_name, str):
            new_row["sample"] = normalize_workflow_sample_name(sample_name)
        normalized_rows.append(new_row)
    return normalized_rows


def rewrite_map_stats_payload(payload: dict[str, object]) -> dict[str, object]:
    """Return a migrated map-stats payload with canonical sample labels."""
    migrated = dict(payload)
    if "applied_mapping_summary" in migrated:
        migrated["applied_mapping_summary"] = _normalize_stats_rows(
            list(migrated["applied_mapping_summary"])
        )

    preview = dict(migrated.get("assemble_read_filter_preview", {}))
    if "filter_effects" in preview:
        preview["filter_effects"] = _normalize_stats_rows(list(preview["filter_effects"]))
    if "metric_summaries" in preview:
        preview["metric_summaries"] = _normalize_stats_rows(list(preview["metric_summaries"]))
    migrated["assemble_read_filter_preview"] = preview
    return migrated


def _rewrite_stats_sidecar(migration: MapStatsMigration) -> None:
    """Rewrite one mapper stats JSON sidecar and regenerate its TXT report."""
    payload = json.loads(migration.source_json.read_text(encoding="utf-8"))
    rewritten_payload = rewrite_map_stats_payload(payload)

    migration.output_json.parent.mkdir(parents=True, exist_ok=True)
    migration.output_json.write_text(
        json.dumps(rewritten_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    migration.output_txt.write_text(
        render_map_stats_payload_report(rewritten_payload),
        encoding="utf-8",
    )


def _log_migration_plan(plan: LegacyMapMigrationPlan, *, dry_run: bool) -> None:
    """Log the planned BAM and sidecar migrations."""
    logger.info(
        "{} {} BAM migration(s), {} stats rewrite(s), {} unrelated BAM(s) skipped",
        "planned" if dry_run else "running",
        len(plan.bam_migrations),
        len(plan.stats_migrations),
        len(plan.skipped_bams),
    )
    for skipped in plan.skipped_bams:
        logger.info("skipping unrelated BAM: {}", skipped.name)
    for orphan_txt in plan.orphan_stats_txt:
        logger.warning(
            "skipping orphan map stats text file without matching JSON: {}",
            orphan_txt.name,
        )
    for migration in plan.bam_migrations:
        logger.info(
            "{} -> {} | header sample {} -> {}",
            migration.source_bam.name,
            migration.output_bam.name,
            migration.legacy_sample,
            migration.canonical_sample,
        )
    for migration in plan.stats_migrations:
        logger.info(
            "stats {} -> {} and {}",
            migration.source_json.name,
            migration.output_json.name,
            migration.output_txt.name,
        )


def migrate_legacy_map_outputs(
    *,
    indir: str | Path,
    outdir: str | Path,
    threads: int = 1,
    dry_run: bool = False,
    force: bool = False,
) -> LegacyMapMigrationPlan:
    """Migrate legacy map BAMs and sidecars into a new output directory."""
    if threads < 1:
        raise IPyradError("--threads must be at least 1.")

    plan = plan_legacy_map_migration(indir, outdir)
    target_paths = [
        path
        for migration in plan.bam_migrations
        for path in (migration.output_bam, migration.output_index)
    ] + [
        path
        for migration in plan.stats_migrations
        for path in (migration.output_json, migration.output_txt)
    ]
    _ensure_targets_writable(target_paths, force=force)
    _log_migration_plan(plan, dry_run=dry_run)

    if dry_run:
        return plan

    resolved_outdir = _resolve_dir(outdir)
    resolved_outdir.mkdir(parents=True, exist_ok=True)

    for migration in plan.bam_migrations:
        logger.info("rewriting BAM RG tags and index for {}", migration.source_bam.name)
        _rewrite_bam_rg_tags(migration, threads=threads)

    for migration in plan.stats_migrations:
        logger.info("rewriting map stats sidecars for {}", migration.source_json.name)
        _rewrite_stats_sidecar(migration)

    logger.info(
        "completed legacy map migration into {} with {} BAM(s) and {} stats file set(s)",
        resolved_outdir,
        len(plan.bam_migrations),
        len(plan.stats_migrations),
    )
    return plan
