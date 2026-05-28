#!/usr/bin/env python

"""Checkpointed per-window tree inference for sequence-backed HDF5 inputs."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from fnmatch import fnmatchcase
from pathlib import Path
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd
from loguru import logger

from ..extracters.window_extracter import (
    MISSING_BASE,
    WindowExtracter,
    count_snps,
    filter_block_by_minmap,
)
from ..extracters.sequence_common import load_sequence_chunk_from_phy
from ...utils.exceptions import IPyradError


FINAL_TREE_STATUSES = {"polytomy_written", "tree_completed"}
TERMINAL_STATUSES = FINAL_TREE_STATUSES | {
    "skipped_no_data",
    "skipped_short_alignment",
    "skipped_few_samples",
    "tree_failed",
}
MANIFEST_COLUMNS = [
    "window_id",
    "window_name",
    "window_mode",
    "scaffold",
    "start",
    "end",
    "first_locus",
    "last_locus",
    "nloci",
    "sites_total",
    "status",
    "status_detail",
    "nsamples_before_filtering",
    "nsites_before_filtering",
    "nvariants_before_filtering",
    "nsites_after_site_filter",
    "nsamples_after_sample_length_filter",
    "nsites_after_sample_length_filter",
    "nsamples_after_filtering",
    "nsites_after_filtering",
    "nvariants_after_filtering",
    "samples_dropped_by_min_sample_alignment_length",
    "retained_sample_names",
    "alignment_path",
    "tree_newick",
    "tree_source",
    "tree_error",
]


@dataclass(frozen=True)
class TreeSliderWindowSpec:
    """One planned sequence window."""

    window_id: int
    window_name: str
    window_mode: str
    scaffold: str
    start: int
    end: int
    first_locus: int | None
    last_locus: int | None
    nloci: int
    sites_total: int
    spans: tuple[tuple[int, int], ...]


def _merge_spans(spans: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Merge contiguous or overlapping spans."""
    if not spans:
        return tuple()

    ordered = sorted(spans)
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return tuple(merged)


def _require_sequence_hdf5(data: Path | str) -> None:
    """Raise one clear error if the input does not look sequence-backed."""
    with h5py.File(data, "r") as io5:
        if "phy" not in io5 or "phymap" not in io5:
            raise IPyradError(
                "treeslider requires a sequence-backed HDF5 with `phy` and `phymap` datasets."
            )
        if io5.attrs["version"] < 2.0:
            raise IPyradError("hdf5 database version must be >= 2.0")


def _normalize_positive_int(value: int | None, label: str) -> int | None:
    """Validate one optional positive integer."""
    if value is None:
        return None
    if int(value) < 1:
        raise IPyradError(f"{label} must be at least 1.")
    return int(value)


def _normalize_nonnegative_int(value: int | None, label: str) -> int | None:
    """Validate one optional non-negative integer."""
    if value is None:
        return None
    if int(value) < 0:
        raise IPyradError(f"{label} must be >= 0.")
    return int(value)


def _resolve_scaffold_subset(
    scaffold_table: pd.DataFrame,
    patterns: list[str] | None,
) -> list[str]:
    """Return selected scaffold names in table order."""
    all_names = scaffold_table["scaffold_name"].tolist()
    if not patterns:
        return all_names

    selected: list[str] = []
    seen: set[str] = set()
    for pattern in patterns:
        matched = [name for name in all_names if fnmatchcase(name, pattern)]
        if not matched:
            raise IPyradError(
                f"No scaffold names match '{pattern}'. Use --print-scaffold-table to inspect available scaffolds."
            )
        for name in matched:
            if name not in seen:
                seen.add(name)
                selected.append(name)
    return selected


def _load_phymap(data: Path | str) -> np.ndarray:
    """Return the phymap matrix as one dense numpy array."""
    with h5py.File(data, "r") as io5:
        return np.asarray(io5["phymap"], dtype=np.int64)


def _plan_locus_windows(
    phymap: np.ndarray,
    scaffold_names: list[str],
    selected_scaffolds: list[str],
) -> list[TreeSliderWindowSpec]:
    """Plan one tree window per phymap locus."""
    allowed = set(selected_scaffolds)
    specs: list[TreeSliderWindowSpec] = []
    for row_idx, row in enumerate(phymap):
        scaffold = scaffold_names[int(row[0])]
        if scaffold not in allowed:
            continue
        start = int(row[3])
        end = int(row[4])
        specs.append(
            TreeSliderWindowSpec(
                window_id=len(specs) + 1,
                window_name=f"{scaffold}:{start}-{end}",
                window_mode="locus",
                scaffold=scaffold,
                start=start,
                end=end,
                first_locus=row_idx + 1,
                last_locus=row_idx + 1,
                nloci=1,
                sites_total=int(row[2] - row[1]),
                spans=((int(row[1]), int(row[2])),),
            )
        )
    return specs


def _plan_genomic_windows(
    phymap: np.ndarray,
    scaffold_names: list[str],
    scaffold_lengths: list[int],
    selected_scaffolds: list[str],
    *,
    window_size: int,
    slide_size: int,
) -> list[TreeSliderWindowSpec]:
    """Plan scaffold-coordinate sliding windows."""
    selected = set(selected_scaffolds)
    specs: list[TreeSliderWindowSpec] = []
    row_numbers = np.arange(phymap.shape[0], dtype=np.int64)

    for scaffold_idx, scaffold in enumerate(scaffold_names):
        if scaffold not in selected:
            continue
        scaffold_rows = phymap[phymap[:, 0] == scaffold_idx]
        scaffold_row_numbers = row_numbers[phymap[:, 0] == scaffold_idx]
        scaffold_length = int(scaffold_lengths[scaffold_idx])
        locus_start0 = (
            scaffold_rows[:, 3].astype(np.int64, copy=False) - 1
            if scaffold_rows.size
            else np.zeros(0, dtype=np.int64)
        )
        locus_end0 = (
            scaffold_rows[:, 4].astype(np.int64, copy=False)
            if scaffold_rows.size
            else np.zeros(0, dtype=np.int64)
        )

        for start0 in range(0, scaffold_length, slide_size):
            end0 = min(start0 + window_size, scaffold_length)
            spans: list[tuple[int, int]] = []
            first_locus = None
            last_locus = None
            nloci = 0

            if scaffold_rows.size:
                overlap_mask = (locus_end0 > start0) & (locus_start0 < end0)
                selected_rows = scaffold_rows[overlap_mask]
                selected_row_numbers = scaffold_row_numbers[overlap_mask]
                nloci = int(selected_rows.shape[0])
                if nloci:
                    first_locus = int(selected_row_numbers[0]) + 1
                    last_locus = int(selected_row_numbers[-1]) + 1
                for row in selected_rows:
                    locus_start = int(row[3]) - 1
                    locus_end = int(row[4])
                    overlap_start = max(locus_start, start0)
                    overlap_end = min(locus_end, end0)
                    if overlap_start >= overlap_end:
                        continue
                    phy_start = int(row[1]) + (overlap_start - locus_start)
                    phy_end = phy_start + (overlap_end - overlap_start)
                    spans.append((phy_start, phy_end))

            merged_spans = _merge_spans(spans)
            specs.append(
                TreeSliderWindowSpec(
                    window_id=len(specs) + 1,
                    window_name=f"{scaffold}:{start0 + 1}-{end0}",
                    window_mode="genomic",
                    scaffold=scaffold,
                    start=start0 + 1,
                    end=end0,
                    first_locus=first_locus,
                    last_locus=last_locus,
                    nloci=nloci,
                    sites_total=sum(stop - start for start, stop in merged_spans),
                    spans=merged_spans,
                )
            )
    return specs


def _spec_record(spec: TreeSliderWindowSpec) -> dict[str, object]:
    """Return one manifest-ready planning record."""
    return {
        "window_id": spec.window_id,
        "window_name": spec.window_name,
        "window_mode": spec.window_mode,
        "scaffold": spec.scaffold,
        "start": spec.start,
        "end": spec.end,
        "first_locus": spec.first_locus if spec.first_locus is not None else "",
        "last_locus": spec.last_locus if spec.last_locus is not None else "",
        "nloci": spec.nloci,
        "sites_total": spec.sites_total,
        "status": "planned",
        "status_detail": "",
        "nsamples_before_filtering": "",
        "nsites_before_filtering": "",
        "nvariants_before_filtering": "",
        "nsites_after_site_filter": "",
        "nsamples_after_sample_length_filter": "",
        "nsites_after_sample_length_filter": "",
        "nsamples_after_filtering": "",
        "nsites_after_filtering": "",
        "nvariants_after_filtering": "",
        "samples_dropped_by_min_sample_alignment_length": "",
        "retained_sample_names": "",
        "alignment_path": "",
        "tree_newick": "",
        "tree_source": "",
        "tree_error": "",
    }


def _write_manifest(manifest: pd.DataFrame, path: Path) -> None:
    """Persist the manifest in stable column order."""
    output = manifest.reset_index(drop=True).copy()
    for column in MANIFEST_COLUMNS:
        if column not in output.columns:
            output[column] = ""
    output = output[MANIFEST_COLUMNS]
    output.to_csv(path, sep="\t", index=False)


def _load_manifest(path: Path) -> pd.DataFrame:
    """Load one existing manifest."""
    manifest = pd.read_csv(path, sep="\t", keep_default_na=False)
    for column in MANIFEST_COLUMNS:
        if column not in manifest.columns:
            manifest[column] = ""
    manifest = manifest[MANIFEST_COLUMNS]
    manifest["window_id"] = manifest["window_id"].astype(int)
    manifest = manifest.set_index("window_id", drop=False)
    return manifest


def _initialize_manifest(specs: list[TreeSliderWindowSpec]) -> pd.DataFrame:
    """Return one new manifest dataframe."""
    manifest = pd.DataFrame([_spec_record(spec) for spec in specs], columns=MANIFEST_COLUMNS)
    manifest["window_id"] = manifest["window_id"].astype(int)
    return manifest.set_index("window_id", drop=False)


def _validate_manifest_against_specs(
    manifest: pd.DataFrame,
    specs: list[TreeSliderWindowSpec],
) -> None:
    """Ensure an existing manifest matches the current requested plan."""
    planned = pd.DataFrame([_spec_record(spec) for spec in specs])[[
        "window_id",
        "window_name",
        "window_mode",
        "scaffold",
        "start",
        "end",
        "first_locus",
        "last_locus",
        "nloci",
        "sites_total",
    ]]
    current = (
        manifest.reset_index(drop=True)[planned.columns]
        .replace({np.nan: ""})
    )
    if len(planned) != len(current) or not planned.equals(current):
        raise IPyradError(
            "Existing treeslider manifest does not match the current requested windows. "
            "Use --force to restart from scratch."
        )


def _clear_tree_fields(manifest: pd.DataFrame, window_id: int) -> None:
    """Reset the tree-related columns for one window."""
    for key in ("tree_newick", "tree_source", "tree_error"):
        manifest.at[window_id, key] = ""


def _quote_newick_label(name: str) -> str:
    """Return one Newick-safe taxon label."""
    if name and all(ch.isalnum() or ch in "._-" for ch in name):
        return name
    escaped = name.replace("'", "''")
    return f"'{escaped}'"


def _build_polytomy(names: list[str]) -> str:
    """Return one multifurcating tree over the retained taxa."""
    labels = ",".join(_quote_newick_label(name) for name in names)
    return f"({labels});"


def _write_fasta(path: Path, names: list[str], seqarr: np.ndarray) -> None:
    """Write one FASTA alignment."""
    with open(path, "w", encoding="utf-8") as out:
        for idx, name in enumerate(names):
            out.write(f">{name}\n{seqarr[idx].tobytes().decode('utf-8')}\n")


def _rebuild_imap_row_indices(
    names: list[str],
    imap: dict[str, list[str]],
) -> dict[str, np.ndarray]:
    """Project one IMAP onto the retained sample order."""
    name_to_idx = {name: idx for idx, name in enumerate(names)}
    return {
        pop: np.array([name_to_idx[name] for name in pop_names if name in name_to_idx], dtype=np.int64)
        for pop, pop_names in imap.items()
    }


def _extract_window_alignment(
    *,
    data: Path,
    extracter: WindowExtracter,
    spec: TreeSliderWindowSpec,
    min_sample_alignment_length: int,
    min_alignment_length: int,
    alignment_path: Path,
) -> dict[str, object]:
    """Filter one planned window and optionally write its staged FASTA."""
    with h5py.File(data, "r") as io5:
        block = load_sequence_chunk_from_phy(io5["phy"], extracter.sidxs, spec.spans)

    nsites_before = int(block.shape[1])
    nvariants_before = count_snps(block)
    site_filtered = extracter._filter_block_sites(block)
    nsites_after_site = int(site_filtered.shape[1])
    if nsites_after_site == 0:
        return {
            "status": "skipped_no_data",
            "status_detail": "No sites remained after min-sample-coverage filtering.",
            "nsamples_before_filtering": len(extracter.snames),
            "nsites_before_filtering": nsites_before,
            "nvariants_before_filtering": nvariants_before,
            "nsites_after_site_filter": 0,
            "nsamples_after_sample_length_filter": 0,
            "nsites_after_sample_length_filter": 0,
            "nsamples_after_filtering": 0,
            "nsites_after_filtering": 0,
            "nvariants_after_filtering": 0,
            "samples_dropped_by_min_sample_alignment_length": "",
            "retained_sample_names": "",
            "alignment_path": "",
        }

    observed_counts = np.sum(site_filtered != MISSING_BASE, axis=1).astype(np.int64, copy=False)
    keep_mask = observed_counts >= min_sample_alignment_length
    kept_names = [extracter.snames[idx] for idx in np.flatnonzero(keep_mask)]
    dropped_names = [extracter.snames[idx] for idx in np.flatnonzero(~keep_mask)]
    sample_filtered = site_filtered[keep_mask, :]
    nsites_after_sample = int(sample_filtered.shape[1]) if sample_filtered.size else 0

    final_alignment = sample_filtered
    if sample_filtered.size and dropped_names:
        final_alignment = filter_block_by_minmap(
            sample_filtered,
            _rebuild_imap_row_indices(kept_names, extracter.imap),
            extracter.minmap,
        )

    nsites_after_final = int(final_alignment.shape[1]) if final_alignment.size else 0
    nvariants_after_final = count_snps(final_alignment) if final_alignment.size else 0
    filtered_names = kept_names

    result = {
        "nsamples_before_filtering": len(extracter.snames),
        "nsites_before_filtering": nsites_before,
        "nvariants_before_filtering": nvariants_before,
        "nsites_after_site_filter": nsites_after_site,
        "nsamples_after_sample_length_filter": len(kept_names),
        "nsites_after_sample_length_filter": nsites_after_sample,
        "nsamples_after_filtering": len(filtered_names),
        "nsites_after_filtering": nsites_after_final,
        "nvariants_after_filtering": nvariants_after_final,
        "samples_dropped_by_min_sample_alignment_length": ",".join(dropped_names),
        "retained_sample_names": ",".join(filtered_names),
    }

    if len(filtered_names) < 3:
        result.update(
            status="skipped_few_samples",
            status_detail="Fewer than 3 samples remained after filtering.",
            alignment_path="",
        )
        return result

    if nsites_after_final < min_alignment_length:
        result.update(
            status="skipped_short_alignment",
            status_detail="Alignment length after filtering is below --min-alignment-length.",
            alignment_path="",
        )
        return result

    alignment_path.parent.mkdir(parents=True, exist_ok=True)
    _write_fasta(alignment_path, filtered_names, final_alignment)
    result.update(
        status="accepted_pending_tree",
        status_detail="Alignment staged for tree inference.",
        alignment_path=str(alignment_path),
        filtered_names=filtered_names,
        filtered_alignment=final_alignment,
    )
    return result


def _resolve_binary(binary: Path | str | None) -> str:
    """Resolve the raxml-ng executable path."""
    if binary is not None:
        path = Path(binary).expanduser()
        if not path.exists():
            raise IPyradError(f"Could not find the requested raxml-ng binary: {path}")
        return str(path)

    prefix_binary = Path(sys.prefix) / "bin" / "raxml-ng"
    if prefix_binary.exists():
        return str(prefix_binary)

    resolved = shutil.which("raxml-ng")
    if not resolved:
        raise IPyradError(
            "Could not find the `raxml-ng` binary. Install it in the active environment, "
            "ensure it is on PATH, or pass an explicit path with `--raxml-ng-binary`."
        )
    return resolved


def _resolve_parallelism(
    *,
    threads: int | str,
    workers: int | str,
    pending_jobs: int,
) -> tuple[int, int]:
    """Resolve `auto` thread and worker settings."""
    available = max(1, os.cpu_count() or 1)
    if pending_jobs < 1:
        return 1, 1

    if threads == "auto" and workers == "auto":
        resolved_threads = 1
        resolved_workers = min(available, pending_jobs)
    elif threads == "auto":
        resolved_workers = int(workers)
        resolved_threads = max(1, available // resolved_workers)
    elif workers == "auto":
        resolved_threads = int(threads)
        resolved_workers = max(1, min(pending_jobs, available // resolved_threads))
    else:
        resolved_threads = int(threads)
        resolved_workers = int(workers)

    return max(1, resolved_threads), max(1, min(resolved_workers, pending_jobs))


def _run_tree_job(
    *,
    binary: str,
    spec: TreeSliderWindowSpec,
    alignment_path: Path,
    stage_dir: Path,
    threads: int,
    model: str,
    bs_trees: int,
    seed: int | None,
    redo: bool,
) -> dict[str, object]:
    """Run one raxml-ng tree job and return one terminal result record."""
    prefix = f"window_{spec.window_id:06d}"
    workdir = stage_dir / "raxml" / prefix
    workdir.mkdir(parents=True, exist_ok=True)

    cmd = [binary]
    if bs_trees > 0:
        cmd.extend(["--all", "--bs-trees", str(bs_trees)])
    else:
        cmd.append("--search")
    cmd.extend(
        [
            "--msa", str(alignment_path),
            "--model", model,
            "--prefix", prefix,
            "--threads", str(threads),
            "--log", "INFO",
        ]
    )
    if seed is not None:
        cmd.extend(["--seed", str(seed + spec.window_id)])
    if redo:
        cmd.append("--redo")

    try:
        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            check=False,
        )
    except Exception as exc:  # pragma: no cover - defensive
        return {
            "window_id": spec.window_id,
            "status": "tree_failed",
            "tree_newick": "",
            "tree_source": "",
            "tree_error": str(exc),
            "workdir": workdir,
            "alignment_path": alignment_path,
        }

    combined = "\n".join(part for part in (result.stdout, result.stderr) if part).strip()
    best_tree = workdir / f"{prefix}.raxml.bestTree"
    support_tree = workdir / f"{prefix}.raxml.support"
    final_tree_path = support_tree if bs_trees > 0 and support_tree.exists() else best_tree
    if result.returncode != 0:
        return {
            "window_id": spec.window_id,
            "status": "tree_failed",
            "tree_newick": "",
            "tree_source": "",
            "tree_error": combined or f"raxml-ng exited with code {result.returncode}.",
            "workdir": workdir,
            "alignment_path": alignment_path,
        }
    if not final_tree_path.exists():
        return {
            "window_id": spec.window_id,
            "status": "tree_failed",
            "tree_newick": "",
            "tree_source": "",
            "tree_error": "raxml-ng completed but did not produce the expected final tree output.",
            "workdir": workdir,
            "alignment_path": alignment_path,
        }

    tree_text = final_tree_path.read_text(encoding="utf-8").strip()
    return {
        "window_id": spec.window_id,
        "status": "tree_completed",
        "tree_newick": tree_text,
        "tree_source": "raxml-ng",
        "tree_error": "",
        "workdir": workdir,
        "alignment_path": alignment_path,
    }


def _cleanup_window_stage(alignment_path: Path | str, workdir: Path | str) -> None:
    """Remove one window's staged alignment and raxml work directory."""
    if alignment_path:
        Path(alignment_path).unlink(missing_ok=True)
    if workdir:
        shutil.rmtree(Path(workdir), ignore_errors=True)


def _prune_empty_dirs(root: Path) -> None:
    """Remove empty stage directories from the bottom up."""
    if not root.exists():
        return
    for path in sorted(root.glob("**/*"), reverse=True):
        if path.is_dir():
            try:
                path.rmdir()
            except OSError:
                continue
    try:
        root.rmdir()
    except OSError:
        pass


def _write_trees_nexus(manifest: pd.DataFrame, path: Path) -> None:
    """Write one multitree Nexus file from terminal tree rows."""
    rows = (
        manifest[manifest["status"].isin(FINAL_TREE_STATUSES)]
        .reset_index(drop=True)
        .sort_values("window_id")
    )
    lines = ["#NEXUS", "Begin trees;"]
    for row in rows.itertuples(index=False):
        lines.append(f"  Tree window_{int(row.window_id):06d} = [&U] {row.tree_newick}")
    lines.append("End;")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_stats_text(
    *,
    path: Path,
    data: Path,
    name: str,
    outdir: Path,
    mode: str,
    selected_scaffolds: list[str],
    binary: str,
    threads: int | str,
    workers: int | str,
    resolved_threads: int,
    resolved_workers: int,
    bs_trees: int,
    model: str,
    min_sample_coverage: int,
    min_sample_alignment_length: int,
    min_alignment_length: int,
    redo: bool,
    force: bool,
    manifest: pd.DataFrame,
) -> None:
    """Write one human-readable summary report."""
    counts = manifest["status"].value_counts().to_dict()
    lines = [
        "Summary",
        "-------",
        f"tool: treeslider",
        f"infile: {data}",
        f"name: {name}",
        f"outdir: {outdir}",
        f"window_mode: {mode}",
        f"scaffolds_selected: {', '.join(selected_scaffolds)}",
        f"raxml_ng_binary: {binary}",
        f"threads_requested: {threads}",
        f"workers_requested: {workers}",
        f"threads_resolved: {resolved_threads}",
        f"workers_resolved: {resolved_workers}",
        f"bs_trees: {bs_trees}",
        f"model: {model}",
        f"min_sample_coverage: {min_sample_coverage}",
        f"min_sample_alignment_length: {min_sample_alignment_length}",
        f"min_alignment_length: {min_alignment_length}",
        f"redo: {redo}",
        f"force: {force}",
        f"windows_planned: {len(manifest)}",
        f"accepted_pending_tree: {counts.get('accepted_pending_tree', 0)}",
        f"polytomy_written: {counts.get('polytomy_written', 0)}",
        f"tree_completed: {counts.get('tree_completed', 0)}",
        f"tree_failed: {counts.get('tree_failed', 0)}",
        f"skipped_no_data: {counts.get('skipped_no_data', 0)}",
        f"skipped_short_alignment: {counts.get('skipped_short_alignment', 0)}",
        f"skipped_few_samples: {counts.get('skipped_few_samples', 0)}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_treeslider_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    window_size: int | None,
    slide_size: int | None,
    print_scaffold_table: bool,
    scaffolds: list[str] | None,
    min_sample_coverage: int,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    min_sample_alignment_length: int,
    min_alignment_length: int,
    threads: int | str,
    workers: int | str,
    bs_trees: int,
    model: str,
    raxml_ng_binary: Path | str | None,
    seed: int | None,
    force: bool,
    redo: bool,
    log_level: str = "INFO",
) -> None:
    """CLI entrypoint for checkpointed sequence-window tree inference."""
    data = Path(data).expanduser().absolute()
    outdir = Path(outdir).expanduser().absolute()
    window_size = _normalize_positive_int(window_size, "--window-size")
    slide_size = _normalize_positive_int(slide_size, "--slide-size")
    min_sample_coverage = _normalize_positive_int(min_sample_coverage, "--min-sample-coverage") or 4
    min_sample_alignment_length = _normalize_positive_int(
        min_sample_alignment_length,
        "--min-sample-alignment-length",
    ) or 1
    min_alignment_length = _normalize_positive_int(
        min_alignment_length,
        "--min-alignment-length",
    ) or 1
    bs_trees = _normalize_nonnegative_int(bs_trees, "--bs-trees") or 0

    if slide_size is not None and window_size is None:
        raise IPyradError("--slide-size requires --window-size.")
    if window_size is None:
        mode = "locus"
        if slide_size is not None:
            raise IPyradError("--slide-size is only valid with --window-size.")
    else:
        mode = "genomic"
        if slide_size is None:
            slide_size = window_size

    _require_sequence_hdf5(data)
    extracter = WindowExtracter(
        data=data,
        name=name,
        outdir=outdir,
        out_format="fa",
        windows=[],
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=1.0,
        exclude=exclude,
        include_reference=include_reference,
        imap=imap,
        minmap=minmap,
        stdout=False,
        force=True,
    )
    if print_scaffold_table:
        extracter.scaffold_table.to_csv(sys.stdout, sep="\t", index=False)
        return

    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / f"{name}.stats.tsv"
    stats_path = outdir / f"{name}.stats.txt"
    trees_path = outdir / f"{name}.trees.nex"
    stage_dir = outdir / f".{name}.stage"
    align_dir = stage_dir / "alignments"

    if force:
        for path in (manifest_path, stats_path, trees_path):
            path.unlink(missing_ok=True)
        shutil.rmtree(stage_dir, ignore_errors=True)

    phymap = _load_phymap(data)
    scaffold_names = extracter.scaffold_table["scaffold_name"].tolist()
    scaffold_lengths = extracter.scaffold_table["scaffold_length"].astype(int).tolist()
    selected_scaffolds = _resolve_scaffold_subset(extracter.scaffold_table, scaffolds)

    if mode == "genomic":
        specs = _plan_genomic_windows(
            phymap,
            scaffold_names,
            scaffold_lengths,
            selected_scaffolds,
            window_size=window_size,
            slide_size=slide_size,
        )
    else:
        specs = _plan_locus_windows(phymap, scaffold_names, selected_scaffolds)
    if not specs:
        raise IPyradError("No windows were planned from the requested scaffold selection.")

    spec_map = {spec.window_id: spec for spec in specs}
    if manifest_path.exists():
        manifest = _load_manifest(manifest_path)
        _validate_manifest_against_specs(manifest, specs)
    else:
        if any(path.exists() for path in (stats_path, trees_path)):
            raise IPyradError(
                f"Existing treeslider outputs were found in {outdir}. Use --force to restart from scratch."
            )
        manifest = _initialize_manifest(specs)
        _write_manifest(manifest, manifest_path)

    if redo:
        failed_ids = manifest.index[manifest["status"] == "tree_failed"]
        for window_id in failed_ids:
            manifest.at[window_id, "status"] = "planned"
            manifest.at[window_id, "status_detail"] = "Retrying failed tree inference."
            manifest.at[window_id, "alignment_path"] = ""
            _clear_tree_fields(manifest, window_id)
        if len(failed_ids):
            _write_manifest(manifest, manifest_path)

    align_dir.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        status = str(manifest.at[spec.window_id, "status"])
        if status in TERMINAL_STATUSES:
            continue
        alignment_path = align_dir / f"window_{spec.window_id:06d}.fa"
        if status == "accepted_pending_tree" and alignment_path.exists():
            continue

        result = _extract_window_alignment(
            data=data,
            extracter=extracter,
            spec=spec,
            min_sample_alignment_length=min_sample_alignment_length,
            min_alignment_length=min_alignment_length,
            alignment_path=alignment_path,
        )
        for key, value in result.items():
            if key in {"filtered_names", "filtered_alignment"}:
                continue
            manifest.at[spec.window_id, key] = value
        if result["status"] != "accepted_pending_tree":
            _clear_tree_fields(manifest, spec.window_id)
        _write_manifest(manifest, manifest_path)

    pending = manifest[manifest["status"] == "accepted_pending_tree"].copy()
    resolved_threads, resolved_workers = _resolve_parallelism(
        threads=threads,
        workers=workers,
        pending_jobs=len(pending),
    )
    binary = _resolve_binary(raxml_ng_binary)

    for window_id in pending.index.tolist():
        spec = spec_map[int(window_id)]
        alignment_path = Path(manifest.at[window_id, "alignment_path"])
        if not alignment_path.exists():
            result = _extract_window_alignment(
                data=data,
                extracter=extracter,
                spec=spec,
                min_sample_alignment_length=min_sample_alignment_length,
                min_alignment_length=min_alignment_length,
                alignment_path=alignment_path,
            )
            for key, value in result.items():
                if key in {"filtered_names", "filtered_alignment"}:
                    continue
                manifest.at[window_id, key] = value
            _write_manifest(manifest, manifest_path)

        if manifest.at[window_id, "status"] != "accepted_pending_tree":
            continue

        nvariants = int(manifest.at[window_id, "nvariants_after_filtering"])
        if nvariants == 0:
            names = [
                value
                for value in str(manifest.at[window_id, "retained_sample_names"]).split(",")
                if value
            ]
            manifest.at[window_id, "status"] = "polytomy_written"
            manifest.at[window_id, "status_detail"] = "No variable sites remained after filtering."
            manifest.at[window_id, "tree_newick"] = _build_polytomy(names)
            manifest.at[window_id, "tree_source"] = "polytomy"
            manifest.at[window_id, "tree_error"] = ""
            _cleanup_window_stage(manifest.at[window_id, "alignment_path"], "")
            manifest.at[window_id, "alignment_path"] = ""
            _write_manifest(manifest, manifest_path)

    runnable = manifest[manifest["status"] == "accepted_pending_tree"].copy()
    if len(runnable):
        with ThreadPoolExecutor(max_workers=resolved_workers) as pool:
            futures = {
                pool.submit(
                    _run_tree_job,
                    binary=binary,
                    spec=spec_map[int(window_id)],
                    alignment_path=Path(manifest.at[window_id, "alignment_path"]),
                    stage_dir=stage_dir,
                    threads=resolved_threads,
                    model=model,
                    bs_trees=bs_trees,
                    seed=seed,
                    redo=redo,
                ): int(window_id)
                for window_id in runnable.index
            }
            for future in as_completed(futures):
                result = future.result()
                window_id = int(result["window_id"])
                manifest.at[window_id, "status"] = result["status"]
                manifest.at[window_id, "status_detail"] = (
                    "Tree inference completed."
                    if result["status"] == "tree_completed"
                    else "Tree inference failed."
                )
                manifest.at[window_id, "tree_newick"] = result["tree_newick"]
                manifest.at[window_id, "tree_source"] = result["tree_source"]
                manifest.at[window_id, "tree_error"] = result["tree_error"]
                _cleanup_window_stage(result["alignment_path"], result["workdir"])
                manifest.at[window_id, "alignment_path"] = ""
                _write_manifest(manifest, manifest_path)

    _write_trees_nexus(manifest, trees_path)
    _write_stats_text(
        path=stats_path,
        data=data,
        name=name,
        outdir=outdir,
        mode=mode,
        selected_scaffolds=selected_scaffolds,
        binary=binary,
        threads=threads,
        workers=workers,
        resolved_threads=resolved_threads,
        resolved_workers=resolved_workers,
        bs_trees=bs_trees,
        model=model,
        min_sample_coverage=min_sample_coverage,
        min_sample_alignment_length=min_sample_alignment_length,
        min_alignment_length=min_alignment_length,
        redo=redo,
        force=force,
        manifest=manifest,
    )

    _prune_empty_dirs(stage_dir)

    logger.info("wrote treeslider window stats to {}", manifest_path)
    logger.info("wrote treeslider summary report to {}", stats_path)
    logger.info("wrote treeslider Nexus trees to {}", trees_path)
