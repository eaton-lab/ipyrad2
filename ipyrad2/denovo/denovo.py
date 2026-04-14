#!/usr/bin/env python

"""Build a denovo reference library by clustering reads and locus consensuses."""

from __future__ import annotations

import os
import pty
import shutil
import sys
import csv
import math
import random
import re
import select
import signal
import subprocess as sp
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from .align import AlignmentRunSummary, write_ordered_consensus_stream_to_file
from .cluster import build_sample_summary, concat_summaries
from .graph import GraphTableSummary, make_global_tables
from .common import (
    DENOVO_MAPPING_FILENAME,
    DENOVO_SAMPLE_GRAPH_SUMMARY_FILENAME,
    DENOVO_STATS_FILENAME,
    get_arm_boundary,
    infer_record_type,
)
from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.parallel import run_pipeline, run_with_pool
from ..utils.pops import expand_imap_patterns, parse_imap
from ..utils.progress import ProgressBar


WORKDIR_NAME = "_denovo_work"
DEFAULT_MAX_DENOVO_SAMPLES = 10

BIN = Path(sys.prefix) / "bin"
BIN_VSEARCH = str(BIN / "vsearch")
BIN_MAFFT = str(BIN / "mafft")


def _validate_runtime_args(
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    cores: int,
    threads: int,
    delim_idx: int,
) -> None:
    """Validate runner arguments for direct library use."""
    if not 0 < within_similarity <= 1:
        raise IPyradError("within_similarity must be > 0 and <= 1")
    if not 0 < across_similarity <= 1:
        raise IPyradError("across_similarity must be > 0 and <= 1")
    if min_derep_size < 1:
        raise IPyradError("min_derep_size must be >= 1")
    if min_length < 1:
        raise IPyradError("min_length must be >= 1")
    if min_merge_overlap < 1:
        raise IPyradError("min_merge_overlap must be >= 1")
    if max_merge_diffs < 0:
        raise IPyradError("max_merge_diffs must be >= 0")
    if cores < 1:
        raise IPyradError("cores must be >= 1")
    if threads < 1:
        raise IPyradError("threads must be >= 1")
    if threads > cores:
        raise IPyradError("threads cannot exceed cores")
    if delim_idx < 1:
        raise IPyradError("delim_idx must be >= 1")


def _validate_required_binaries() -> str:
    """Validate one required denovo runtime binary in the active env."""
    for binary, name in [(BIN_VSEARCH, "vsearch"), (BIN_MAFFT, "mafft")]:
        path = Path(binary)
        if not (path.exists() and path.is_file() and os.access(path, os.X_OK)):
            raise IPyradError(f"{name} binary is not executable: {path}")


def _validate_fastq_layout(
    fastq_dict: dict[str, tuple[Path, Path | None]],
) -> bool:
    """Validate parsed FASTQ layout and return whether the run is paired-end."""
    if not fastq_dict:
        raise IPyradError("No FASTQ inputs were parsed for denovo.")

    paired_states = {paths[1] is not None for paths in fastq_dict.values()}
    if len(paired_states) > 1:
        raise IPyradError(
            "denovo requires all inputs to be consistently single-end or paired-end; "
            "mixed input layouts were detected."
        )
    return paired_states.pop()


def _iter_denovo_outputs(outdir: Path) -> list[Path]:
    """Return the curated denovo output paths for one output directory."""
    return [
        outdir / "denovo_reference.fa",
        outdir / DENOVO_MAPPING_FILENAME,
        outdir / DENOVO_STATS_FILENAME,
        outdir / DENOVO_SAMPLE_GRAPH_SUMMARY_FILENAME,
        outdir / "denovo.stats.txt",
        outdir / "denovo.audit",
    ]


def _sample_total_fastq_size(paths: tuple[Path, Path | None]) -> int:
    """Return total on-disk FASTQ size for one parsed sample."""
    total = 0
    for path in paths:
        if path is None:
            continue
        try:
            total += path.stat().st_size
        except OSError:
            continue
    return total


def _select_denovo_samples(
    fastq_dict: dict[str, tuple[Path, Path | None]],
    *,
    imap_path: Path | None,
    use_all_samples: bool,
) -> tuple[dict[str, tuple[Path, Path | None]], str]:
    """Select the subset of samples used to build the denovo pseudoreference."""
    if use_all_samples:
        return dict(fastq_dict), "all"

    size_map = {
        sample: _sample_total_fastq_size(paths)
        for sample, paths in fastq_dict.items()
    }

    if imap_path is not None:
        imap = parse_imap(imap_path)
        imap, _unmatched = expand_imap_patterns(
            imap,
            fastq_dict,
            mapping_name="IMAP",
            available_name="the denovo inputs",
        )
        selected_names: set[str] = set()
        for group, names in sorted(imap.items()):
            logger.info(
                "denovo IMAP group '{}' -> matched {} sample(s)",
                group,
                len(names),
            )
            selected_names.update(names)
        if not selected_names:
            raise IPyradError("IMAP did not yield any denovo input samples.")
        selected = {
            name: fastq_dict[name]
            for name in fastq_dict
            if name in selected_names
        }
        return selected, "imap"

    if len(fastq_dict) <= DEFAULT_MAX_DENOVO_SAMPLES:
        return dict(fastq_dict), "all"

    ranked = sorted(
        fastq_dict,
        key=lambda sample: (-size_map[sample], sample),
    )
    eligible_count = max(1, (len(ranked) + 1) // 2)
    eligible = ranked[:eligible_count]
    if len(eligible) < DEFAULT_MAX_DENOVO_SAMPLES:
        shortfall = DEFAULT_MAX_DENOVO_SAMPLES - len(eligible)
        selected_names = eligible + ranked[eligible_count : eligible_count + shortfall]
    elif len(eligible) == DEFAULT_MAX_DENOVO_SAMPLES:
        selected_names = eligible
    else:
        rng = random.Random(0)
        selected_names = sorted(rng.sample(eligible, DEFAULT_MAX_DENOVO_SAMPLES))
    selected = {name: fastq_dict[name] for name in selected_names}
    return selected, "top-half-random"


def _log_sample_selection(
    full_fastq_dict: dict[str, tuple[Path, Path | None]],
    selected_fastq_dict: dict[str, tuple[Path, Path | None]],
    selection_mode: str,
) -> None:
    """Log the denovo sample subset used for pseudoreference construction."""
    dropped = sorted(set(full_fastq_dict).difference(selected_fastq_dict))
    if selection_mode == "all":
        logger.info("using all {} denovo input samples", len(selected_fastq_dict))
        return

    logger.info(
        "selected {} of {} denovo input samples using selection mode '{}'",
        len(selected_fastq_dict),
        len(full_fastq_dict),
        selection_mode,
    )
    logger.debug("selected samples: {}", ", ".join(sorted(selected_fastq_dict)))
    if dropped:
        logger.debug("excluded samples: {}", ", ".join(dropped))


def _cleanup_post_derep_work_files(
    sname: str,
    outdir: Path,
    *,
    paired: bool,
) -> None:
    """Remove large pre-derep sample files once dereplication succeeds."""
    paths = [
        outdir / f"{sname}.joined.fa",
        outdir / f"{sname}.cluster.fa",
    ]
    if paired:
        paths.extend(
            [
                outdir / f"{sname}.merged.fa",
                outdir / f"{sname}.unmerged_R1.fq",
                outdir / f"{sname}.unmerged_R2.fq",
            ]
        )
    for path in paths:
        path.unlink(missing_ok=True)


def _cleanup_post_summary_work_files(
    sname: str,
    outdir: Path,
) -> None:
    """Remove derep/cluster outputs once the sample summary exists."""
    for path in [
        outdir / f"{sname}.consensus.fa",
        outdir / f"{sname}.clusters.tsv",
        outdir / f"{sname}.derep.fa",
        outdir / f"{sname}.derep.sizesorted.fa",
    ]:
        path.unlink(missing_ok=True)


def _cleanup_denovo_workdir_stage_files(workdir: Path) -> None:
    """Remove cross-sample intermediates that are no longer needed."""
    for path in sorted(workdir.glob("*.summary.tsv")):
        path.unlink(missing_ok=True)
    for path in [
        workdir / "consensus.concat.fa",
        workdir / "global_hits.uc.tsv",
        workdir / "concat.summary.tsv",
    ]:
        path.unlink(missing_ok=True)


def _iter_status_records(buffer: str) -> tuple[list[str], str]:
    """Return complete status records split on CR/LF plus trailing partial text."""
    parts = re.split(r"\r\n|\r|\n", buffer)
    return [part for part in parts[:-1] if part], parts[-1]


def _extract_searching_percent(status: str) -> int | None:
    """Return the live Searching percent from one vsearch status line."""
    match = re.match(r"Searching\s+(\d+(?:\.\d+)?)%", status.strip())
    if match is None:
        return None
    return min(100, int(float(match.group(1))))


def _open_vsearch_process_with_stderr_stream(
    cmd: list[str],
) -> tuple[sp.Popen[bytes], int, bool]:
    """Start vsearch and return `(proc, stderr_fd, uses_pty)`."""
    try:
        master_fd, slave_fd = pty.openpty()
    except OSError:
        proc = sp.Popen(
            cmd,
            stdout=sp.DEVNULL,
            stderr=sp.PIPE,
            stdin=sp.DEVNULL,
            start_new_session=True,
        )
        assert proc.stderr is not None
        return proc, proc.stderr.fileno(), False

    try:
        proc = sp.Popen(
            cmd,
            stdout=sp.DEVNULL,
            stderr=slave_fd,
            stdin=sp.DEVNULL,
            start_new_session=True,
            close_fds=True,
        )
    finally:
        os.close(slave_fd)
    return proc, master_fd, True


def _iter_vsearch_stderr_chunks(stderr_fd: int) -> Iterator[str]:
    """Yield decoded stderr chunks live from one vsearch run."""
    while True:
        readable, _, _ = select.select([stderr_fd], [], [], 0.5)
        if not readable:
            continue
        try:
            chunk = os.read(stderr_fd, 4096)
        except OSError:
            break
        if not chunk:
            break
        yield chunk.decode("utf-8", "replace")


def _run_vsearch_with_progress(cmd: list[str], *, message: str) -> None:
    """Run one vsearch command while mirroring live Searching progress."""
    progress = ProgressBar(100, None, message)
    progress.update()
    proc: sp.Popen[bytes] | None = None
    stderr_fd: int | None = None
    uses_pty = False
    last_percent = 0
    stderr_parts: list[str] = []
    try:
        proc, stderr_fd, uses_pty = _open_vsearch_process_with_stderr_stream(cmd)
        # vsearch emits incremental carriage-return progress only on a terminal-like stderr.
        pending = ""
        for raw_text in _iter_vsearch_stderr_chunks(stderr_fd):
            stderr_parts.append(raw_text)
            complete, pending = _iter_status_records(pending + raw_text)
            for status in complete:
                percent = _extract_searching_percent(status)
                if percent is None or percent <= last_percent:
                    continue
                last_percent = percent
                progress.finished = percent
                progress.update()
        if pending:
            percent = _extract_searching_percent(pending)
            if percent is not None and percent > last_percent:
                last_percent = percent
                progress.finished = percent
                progress.update()
        rc = proc.wait()
    except KeyboardInterrupt:
        progress.close()
        logger.warning("interrupted by user. Cleaning up.")
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except Exception:
                pass
        raise SystemExit(130)
    finally:
        if proc is not None and proc.poll() == 0 and last_percent < 100:
            progress.finished = 100
            progress.update()
        if uses_pty and stderr_fd is not None:
            try:
                os.close(stderr_fd)
            except OSError:
                pass
        progress.close()

    if rc != 0:
        stderr_text = "".join(stderr_parts)
        raise RuntimeError(f"command failed rc={rc}: {' '.join(cmd)} stderr={stderr_text.strip()}")


def _prepare_output_paths(
    outdir: Path,
    force: bool,
) -> tuple[Path, dict[str, Path]]:
    """Prepare curated denovo outputs and the internal working directory."""
    outdir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "reference": outdir / "denovo_reference.fa",
        "mapping": outdir / DENOVO_MAPPING_FILENAME,
        "loci_stats": outdir / DENOVO_STATS_FILENAME,
        "sample_graph_summary": outdir / DENOVO_SAMPLE_GRAPH_SUMMARY_FILENAME,
        "run_stats": outdir / "denovo.stats.txt",
        "audit_dir": outdir / "denovo.audit",
        "workdir": outdir / WORKDIR_NAME,
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing and not force:
        joined = ", ".join(path.name for path in existing)
        raise IPyradError(
            f"denovo outputs already exist in {outdir}. Use --force to overwrite: {joined}"
        )

    if force:
        for path in _iter_denovo_outputs(outdir):
            if path.exists():
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        if outputs["workdir"].exists():
            shutil.rmtree(outputs["workdir"])

    outputs["workdir"].mkdir(parents=True, exist_ok=True)
    return outputs["workdir"], outputs


@dataclass(frozen=True, slots=True)
class SelectedSampleQc:
    """Compact per-selected-sample denovo burden metrics."""

    sample: str
    consensus_records: int
    n_reads_sum: int
    joined_records: int
    merged_records: int
    single_records: int


@dataclass(frozen=True, slots=True)
class DenovoQcSummary:
    """Streamed QC metrics derived from final denovo outputs."""

    selected_sample_count: int
    total_input_sample_count: int
    consensus_records: int
    loci_written: int
    singleton_loci: int
    singleton_locus_fraction: float
    loci_with_2plus_samples: int
    loci_with_half_or_more_selected_samples: int
    loci_with_all_selected_samples: int
    mean_samples_per_locus: float
    median_samples_per_locus: float
    mean_cores_per_locus: float
    median_cores_per_locus: float
    max_samples_per_locus: int
    max_cores_per_locus: int
    multi_core_single_sample_loci: int
    duplicated_component_loci: int
    reconciled_loci: int
    audited_components: int
    processed_components: int
    oversize_unsplit_components: int
    largest_component_nodes: int
    component_input_nodes_p50: float
    component_input_nodes_p90: float
    component_input_nodes_p99: float
    component_input_nodes_max: int
    component_contracted_nodes_p50: float
    component_contracted_nodes_p90: float
    component_contracted_nodes_p99: float
    component_contracted_nodes_max: int
    occupancy_counts: tuple[tuple[int, int], ...]
    selected_sample_rows: tuple[SelectedSampleQc, ...]


def _parse_bool_text(value: str | None) -> bool:
    """Return a boolean parsed from TSV text."""
    return str(value).strip().lower() in {"1", "true", "yes"}


def _counter_value_at_rank(counts: Counter[int], rank: int) -> int:
    """Return the 0-based ranked value from a histogram counter."""
    seen = 0
    for value, count in sorted(counts.items()):
        seen += count
        if rank < seen:
            return int(value)
    return 0


def _counter_quantile(counts: Counter[int], total: int, quantile: float) -> float:
    """Return one linear-interpolated quantile from histogram counts."""
    if total <= 0:
        return 0.0
    if total == 1:
        return float(next(iter(counts)))
    position = (total - 1) * quantile
    low_rank = int(math.floor(position))
    high_rank = int(math.ceil(position))
    low_value = _counter_value_at_rank(counts, low_rank)
    high_value = _counter_value_at_rank(counts, high_rank)
    if low_rank == high_rank:
        return float(low_value)
    fraction = position - low_rank
    return float(low_value + (high_value - low_value) * fraction)


def _collect_denovo_qc(
    *,
    selected_fastq_dict: dict[str, tuple[Path, Path | None]],
    total_input_sample_count: int,
    graph_summary: GraphTableSummary,
    outputs: dict[str, Path],
) -> DenovoQcSummary:
    """Collect low-memory QC metrics from final denovo output tables."""
    selected_sample_count = len(selected_fastq_dict)
    half_selected_cutoff = max(1, math.ceil(selected_sample_count / 2))

    loci_total = 0
    singleton_loci = 0
    loci_with_2plus_samples = 0
    loci_with_half_or_more_selected_samples = 0
    loci_with_all_selected_samples = 0
    multi_core_single_sample_loci = 0
    duplicated_component_loci = 0
    reconciled_loci = 0
    n_samples_sum = 0
    n_cores_sum = 0
    max_samples_per_locus = 0
    max_cores_per_locus = 0
    occupancy_counts: Counter[int] = Counter()
    n_cores_counts: Counter[int] = Counter()

    with open(outputs["loci_stats"], "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            n_samples = int(row["n_samples"])
            n_cores = int(row["n_cores"])
            loci_total += 1
            occupancy_counts[n_samples] += 1
            n_cores_counts[n_cores] += 1
            n_samples_sum += n_samples
            n_cores_sum += n_cores
            max_samples_per_locus = max(max_samples_per_locus, n_samples)
            max_cores_per_locus = max(max_cores_per_locus, n_cores)
            singleton_loci += int(n_samples == 1)
            loci_with_2plus_samples += int(n_samples >= 2)
            loci_with_half_or_more_selected_samples += int(n_samples >= half_selected_cutoff)
            loci_with_all_selected_samples += int(n_samples >= selected_sample_count)
            multi_core_single_sample_loci += int(n_samples == 1 and n_cores > 1)
            duplicated_component_loci += int(_parse_bool_text(row.get("duplicated_component")))
            reconciled_loci += int(_parse_bool_text(row.get("used_reconciliation")))

    audited_components = 0
    processed_components = 0
    oversize_unsplit_components = 0
    largest_component_nodes = 0
    component_input_counts: Counter[int] = Counter()
    component_contracted_counts: Counter[int] = Counter()
    components_summary_path = outputs["audit_dir"] / "components.summary.tsv"
    if components_summary_path.exists():
        with open(components_summary_path, "rt", encoding="utf-8", newline="") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                audited_components += 1
                n_input_nodes = int(row["n_input_nodes"])
                n_contracted_nodes = int(row["n_contracted_nodes"])
                component_input_counts[n_input_nodes] += 1
                component_contracted_counts[n_contracted_nodes] += 1
                largest_component_nodes = max(largest_component_nodes, n_input_nodes)
                status = str(row.get("status", ""))
                processed_components += int(status == "processed")
                oversize_unsplit_components += int(status == "oversize_unsplit")

    sample_counts = {
        sample: {
            "consensus_records": 0,
            "n_reads_sum": 0,
            "joined_records": 0,
            "merged_records": 0,
            "single_records": 0,
        }
        for sample in selected_fastq_dict
    }
    concat_summary_tsv = outputs["workdir"] / "concat.summary.tsv"
    if concat_summary_tsv.exists():
        with open(concat_summary_tsv, "rt", encoding="utf-8", newline="") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                sample = str(row["sample"])
                counts = sample_counts.setdefault(
                    sample,
                    {
                        "consensus_records": 0,
                        "n_reads_sum": 0,
                        "joined_records": 0,
                        "merged_records": 0,
                        "single_records": 0,
                    },
                )
                counts["consensus_records"] += 1
                counts["n_reads_sum"] += int(row["n_reads"])
                record_type = str(row["record_type"])
                if record_type == "joined":
                    counts["joined_records"] += 1
                elif record_type == "merged":
                    counts["merged_records"] += 1
                elif record_type == "single":
                    counts["single_records"] += 1

    selected_sample_rows = tuple(
        SelectedSampleQc(
            sample=sample,
            consensus_records=int(counts["consensus_records"]),
            n_reads_sum=int(counts["n_reads_sum"]),
            joined_records=int(counts["joined_records"]),
            merged_records=int(counts["merged_records"]),
            single_records=int(counts["single_records"]),
        )
        for sample, counts in sorted(
            sample_counts.items(),
            key=lambda item: (-int(item[1]["consensus_records"]), item[0]),
        )
    )

    return DenovoQcSummary(
        selected_sample_count=selected_sample_count,
        total_input_sample_count=total_input_sample_count,
        consensus_records=int(graph_summary.consensus_records),
        loci_written=int(graph_summary.loci_written),
        singleton_loci=singleton_loci,
        singleton_locus_fraction=(float(singleton_loci / loci_total) if loci_total else 0.0),
        loci_with_2plus_samples=loci_with_2plus_samples,
        loci_with_half_or_more_selected_samples=loci_with_half_or_more_selected_samples,
        loci_with_all_selected_samples=loci_with_all_selected_samples,
        mean_samples_per_locus=(float(n_samples_sum / loci_total) if loci_total else 0.0),
        median_samples_per_locus=_counter_quantile(occupancy_counts, loci_total, 0.5),
        mean_cores_per_locus=(float(n_cores_sum / loci_total) if loci_total else 0.0),
        median_cores_per_locus=_counter_quantile(n_cores_counts, loci_total, 0.5),
        max_samples_per_locus=max_samples_per_locus,
        max_cores_per_locus=max_cores_per_locus,
        multi_core_single_sample_loci=multi_core_single_sample_loci,
        duplicated_component_loci=duplicated_component_loci,
        reconciled_loci=reconciled_loci,
        audited_components=audited_components,
        processed_components=processed_components,
        oversize_unsplit_components=oversize_unsplit_components,
        largest_component_nodes=largest_component_nodes,
        component_input_nodes_p50=_counter_quantile(component_input_counts, audited_components, 0.5),
        component_input_nodes_p90=_counter_quantile(component_input_counts, audited_components, 0.9),
        component_input_nodes_p99=_counter_quantile(component_input_counts, audited_components, 0.99),
        component_input_nodes_max=max(component_input_counts, default=0),
        component_contracted_nodes_p50=_counter_quantile(component_contracted_counts, audited_components, 0.5),
        component_contracted_nodes_p90=_counter_quantile(component_contracted_counts, audited_components, 0.9),
        component_contracted_nodes_p99=_counter_quantile(component_contracted_counts, audited_components, 0.99),
        component_contracted_nodes_max=max(component_contracted_counts, default=0),
        occupancy_counts=tuple(sorted((int(k), int(v)) for k, v in occupancy_counts.items())),
        selected_sample_rows=selected_sample_rows,
    )


_DENOVO_SAMPLE_GRAPH_SUMMARY_FIELDS = [
    "sample",
    "components_seen",
    "split_components_seen",
    "prop_split_components_seen",
    "multi_subcomponent_components",
    "prop_multi_subcomponent_components",
    "duplicated_components_seen",
    "prop_duplicated_components_seen",
    "reconciled_components_seen",
    "prop_reconciled_components_seen",
]


def _validate_tsv_columns(
    fieldnames: Iterable[str] | None,
    *,
    required: set[str],
    label: str,
) -> None:
    """Validate that one TSV contains the required columns."""
    present = set(fieldnames or ())
    missing = sorted(required.difference(present))
    if missing:
        raise IPyradError(
            f"{label} is missing required columns: {', '.join(missing)}"
        )


def _build_denovo_sample_graph_summary_rows(
    *,
    mapping_tsv: Path,
    loci_stats_tsv: Path,
    sample_names: Iterable[str] | None = None,
) -> list[dict[str, str]]:
    """Summarize per-sample graph-splitting burden from final denovo tables."""
    component_subcomponents: dict[int, set[int]] = defaultdict(set)
    component_duplicated: dict[int, bool] = {}
    component_reconciled: dict[int, bool] = {}

    with open(loci_stats_tsv, "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        _validate_tsv_columns(
            reader.fieldnames,
            required={
                "component_id",
                "subcomponent_id",
                "duplicated_component",
                "used_reconciliation",
            },
            label=loci_stats_tsv.name,
        )
        for row in reader:
            component_id = int(row["component_id"])
            component_subcomponents[component_id].add(int(row["subcomponent_id"]))
            component_duplicated[component_id] = _parse_bool_text(
                row.get("duplicated_component")
            )
            component_reconciled[component_id] = _parse_bool_text(
                row.get("used_reconciliation")
            )

    sample_component_subcomponents: dict[str, dict[int, set[int]]] = defaultdict(
        lambda: defaultdict(set)
    )
    with open(mapping_tsv, "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        _validate_tsv_columns(
            reader.fieldnames,
            required={"sample", "component_id", "subcomponent_id"},
            label=mapping_tsv.name,
        )
        for row in reader:
            sample_component_subcomponents[str(row["sample"])][int(row["component_id"])].add(
                int(row["subcomponent_id"])
            )

    all_samples = set(sample_component_subcomponents)
    if sample_names is not None:
        all_samples.update(str(name) for name in sample_names)

    rows: list[dict[str, str | int | float]] = []
    for sample in all_samples:
        component_map = sample_component_subcomponents.get(sample, {})
        components_seen = len(component_map)
        split_components_seen = sum(
            1
            for component_id in component_map
            if len(component_subcomponents.get(component_id, {0})) > 1
        )
        multi_subcomponent_components = sum(
            1
            for subcomponent_ids in component_map.values()
            if len(subcomponent_ids) > 1
        )
        duplicated_components_seen = sum(
            1
            for component_id in component_map
            if component_duplicated.get(component_id, False)
        )
        reconciled_components_seen = sum(
            1
            for component_id in component_map
            if component_reconciled.get(component_id, False)
        )
        rows.append(
            {
                "sample": sample,
                "components_seen": components_seen,
                "split_components_seen": split_components_seen,
                "prop_split_components_seen": _safe_fraction(
                    split_components_seen, components_seen
                ),
                "multi_subcomponent_components": multi_subcomponent_components,
                "prop_multi_subcomponent_components": _safe_fraction(
                    multi_subcomponent_components, components_seen
                ),
                "duplicated_components_seen": duplicated_components_seen,
                "prop_duplicated_components_seen": _safe_fraction(
                    duplicated_components_seen, components_seen
                ),
                "reconciled_components_seen": reconciled_components_seen,
                "prop_reconciled_components_seen": _safe_fraction(
                    reconciled_components_seen, components_seen
                ),
            }
        )

    rows.sort(
        key=lambda row: (
            -float(row["prop_multi_subcomponent_components"]),
            -float(row["prop_split_components_seen"]),
            str(row["sample"]),
        )
    )
    return [
        {
            "sample": str(row["sample"]),
            "components_seen": str(int(row["components_seen"])),
            "split_components_seen": str(int(row["split_components_seen"])),
            "prop_split_components_seen": _format_report_fraction(
                float(row["prop_split_components_seen"])
            ),
            "multi_subcomponent_components": str(
                int(row["multi_subcomponent_components"])
            ),
            "prop_multi_subcomponent_components": _format_report_fraction(
                float(row["prop_multi_subcomponent_components"])
            ),
            "duplicated_components_seen": str(
                int(row["duplicated_components_seen"])
            ),
            "prop_duplicated_components_seen": _format_report_fraction(
                float(row["prop_duplicated_components_seen"])
            ),
            "reconciled_components_seen": str(
                int(row["reconciled_components_seen"])
            ),
            "prop_reconciled_components_seen": _format_report_fraction(
                float(row["prop_reconciled_components_seen"])
            ),
        }
        for row in rows
    ]


def _write_denovo_sample_graph_summary(
    outpath: Path,
    *,
    mapping_tsv: Path,
    loci_stats_tsv: Path,
    sample_names: Iterable[str] | None = None,
) -> Path:
    """Write one per-sample denovo graph burden summary TSV."""
    rows = _build_denovo_sample_graph_summary_rows(
        mapping_tsv=mapping_tsv,
        loci_stats_tsv=loci_stats_tsv,
        sample_names=sample_names,
    )
    with open(outpath, "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=_DENOVO_SAMPLE_GRAPH_SUMMARY_FIELDS,
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"wrote denovo sample graph summary to {outpath}")
    return outpath


def write_denovo_sample_graph_summary(
    outdir: Path,
    *,
    sample_names: Iterable[str] | None = None,
) -> Path:
    """Write the denovo sample graph summary for one completed output dir."""
    outdir = Path(outdir)
    return _write_denovo_sample_graph_summary(
        outdir / DENOVO_SAMPLE_GRAPH_SUMMARY_FILENAME,
        mapping_tsv=outdir / DENOVO_MAPPING_FILENAME,
        loci_stats_tsv=outdir / DENOVO_STATS_FILENAME,
        sample_names=sample_names,
    )


def _safe_fraction(numer: int, denom: int) -> float:
    """Return `numer / denom` or 0.0 when the denominator is empty."""
    if denom <= 0:
        return 0.0
    return float(numer / denom)


def _format_report_count(value: int) -> str:
    """Format integer counts for the human-readable denovo report."""
    return f"{int(value):,}"


def _format_report_float(value: float, digits: int = 3) -> str:
    """Format floating-point summary values for the human-readable report."""
    return f"{float(value):.{digits}f}"


def _format_report_fraction(value: float) -> str:
    """Format fraction-like values consistently in the human-readable report."""
    return f"{float(value):.6f}"


def _format_report_quantile(value: int | float) -> str:
    """Format component-node quantiles without forcing trailing decimals."""
    if float(value).is_integer():
        return _format_report_count(int(value))
    return _format_report_float(float(value))


def _append_report_key_value_section(
    lines: list[str],
    title: str,
    rows: list[tuple[str, str]],
) -> None:
    """Append one assemble-style aligned key/value section."""
    lines.append(f"# {title}")
    if rows:
        width = max(len(key) for key, _value in rows)
        for key, value in rows:
            lines.append(f"{key.ljust(width)}  {value}")
    lines.append("")


def _append_report_table_section(
    lines: list[str],
    title: str,
    headers: list[str],
    rows: list[list[str]],
) -> None:
    """Append one assemble-style whitespace-aligned table section."""
    lines.append(f"# {title}")
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    lines.append("  ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)))
    for row in rows:
        lines.append("  ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)))
    lines.append("")


def _write_denovo_stats(
    outpath: Path,
    *,
    all_fastq_dict: dict[str, tuple[Path, Path | None]],
    selected_fastq_dict: dict[str, tuple[Path, Path | None]],
    selection_mode: str,
    paired: bool,
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    allow_reverse_complement: bool,
    cores: int,
    threads: int,
    workers: int,
    alignment_summary: AlignmentRunSummary,
    keep_intermediates: bool,
    workdir: Path,
    graph_summary: GraphTableSummary,
    qc_summary: DenovoQcSummary,
    outputs: dict[str, Path],
) -> None:
    """Write a human-readable summary of the denovo run."""
    occupancy_lookup = dict(qc_summary.occupancy_counts)

    inputs_rows = [
        (
            "fastq_count",
            _format_report_count(
                sum(2 if paths[1] else 1 for paths in selected_fastq_dict.values())
            ),
        ),
        ("selected_sample_count", _format_report_count(len(selected_fastq_dict))),
        ("total_input_sample_count", _format_report_count(len(all_fastq_dict))),
        ("sample_selection_mode", selection_mode),
        ("paired_mode", "paired-end" if paired else "single-end"),
    ]
    clustering_rows = [
        ("within_similarity", _format_report_fraction(within_similarity)),
        ("across_similarity", _format_report_fraction(across_similarity)),
        ("min_derep_size", _format_report_count(min_derep_size)),
        ("min_length", _format_report_count(min_length)),
        ("min_merge_overlap", _format_report_count(min_merge_overlap)),
        ("max_merge_diffs", _format_report_count(max_merge_diffs)),
        ("allow_reverse_complement", str(bool(allow_reverse_complement))),
    ]
    summary_rows = [
        ("consensus_records", _format_report_count(graph_summary.consensus_records)),
        ("loci_written", _format_report_count(graph_summary.loci_written)),
        ("single_sequence_loci", _format_report_count(alignment_summary.single_sequence_loci)),
        ("identical_sequence_loci", _format_report_count(alignment_summary.identical_sequence_loci)),
        ("mafft_required_loci", _format_report_count(alignment_summary.mafft_required_loci)),
        ("joined_spacer_loci", _format_report_count(alignment_summary.joined_spacer_loci)),
        (
            "mixed_reconciled_spacer_loci",
            _format_report_count(alignment_summary.mixed_reconciled_spacer_loci),
        ),
        ("stripped_output_loci", _format_report_count(alignment_summary.stripped_output_loci)),
        ("duplicated_components_seen", _format_report_count(graph_summary.duplicated_components_seen)),
        (
            "same_sample_reconciliation_attempted",
            _format_report_count(graph_summary.same_sample_reconciliation_attempted),
        ),
        ("components_reconciled", _format_report_count(graph_summary.components_reconciled)),
        (
            "joined_only_reconciled_loci",
            _format_report_count(graph_summary.joined_only_reconciled_loci),
        ),
        ("mixed_reconciled_loci", _format_report_count(graph_summary.mixed_reconciled_loci)),
        ("mixed_reconciled_groups", _format_report_count(graph_summary.mixed_reconciled_groups)),
    ]
    locus_qc_rows = [
        ("singleton_loci", _format_report_count(qc_summary.singleton_loci)),
        (
            "singleton_locus_fraction",
            _format_report_fraction(qc_summary.singleton_locus_fraction),
        ),
        ("loci_with_2plus_samples", _format_report_count(qc_summary.loci_with_2plus_samples)),
        (
            "loci_with_half_or_more_selected_samples",
            _format_report_count(qc_summary.loci_with_half_or_more_selected_samples),
        ),
        (
            "loci_with_all_selected_samples",
            _format_report_count(qc_summary.loci_with_all_selected_samples),
        ),
        ("mean_samples_per_locus", _format_report_float(qc_summary.mean_samples_per_locus)),
        ("median_samples_per_locus", _format_report_float(qc_summary.median_samples_per_locus)),
        ("max_samples_per_locus", _format_report_count(qc_summary.max_samples_per_locus)),
        ("mean_cores_per_locus", _format_report_float(qc_summary.mean_cores_per_locus)),
        ("median_cores_per_locus", _format_report_float(qc_summary.median_cores_per_locus)),
        ("max_cores_per_locus", _format_report_count(qc_summary.max_cores_per_locus)),
        (
            "multi_core_single_sample_loci",
            _format_report_count(qc_summary.multi_core_single_sample_loci),
        ),
        (
            "duplicated_component_loci",
            _format_report_count(qc_summary.duplicated_component_loci),
        ),
        ("reconciled_loci", _format_report_count(qc_summary.reconciled_loci)),
    ]
    component_qc_rows = [
        ("audited_components", _format_report_count(qc_summary.audited_components)),
        ("processed_components", _format_report_count(qc_summary.processed_components)),
        (
            "oversize_unsplit_components",
            _format_report_count(qc_summary.oversize_unsplit_components),
        ),
        ("largest_component_nodes", _format_report_count(qc_summary.largest_component_nodes)),
    ]
    component_node_headers = ["quantile", "input_nodes", "contracted_nodes"]
    component_node_rows = [
        [
            "p50",
            _format_report_quantile(qc_summary.component_input_nodes_p50),
            _format_report_quantile(qc_summary.component_contracted_nodes_p50),
        ],
        [
            "p90",
            _format_report_quantile(qc_summary.component_input_nodes_p90),
            _format_report_quantile(qc_summary.component_contracted_nodes_p90),
        ],
        [
            "p99",
            _format_report_quantile(qc_summary.component_input_nodes_p99),
            _format_report_quantile(qc_summary.component_contracted_nodes_p99),
        ],
        [
            "max",
            _format_report_count(qc_summary.component_input_nodes_max),
            _format_report_count(qc_summary.component_contracted_nodes_max),
        ],
    ]
    selected_sample_headers = [
        "sample",
        "consensus_records",
        "n_reads_sum",
        "joined_records",
        "merged_records",
        "single_records",
    ]
    selected_sample_rows = [
        [
            row.sample,
            _format_report_count(row.consensus_records),
            _format_report_count(row.n_reads_sum),
            _format_report_count(row.joined_records),
            _format_report_count(row.merged_records),
            _format_report_count(row.single_records),
        ]
        for row in qc_summary.selected_sample_rows
    ]
    occupancy_headers = ["samples_with_data", "loci", "fraction_of_final_loci"]
    occupancy_rows = [
        [
            _format_report_count(n_samples),
            _format_report_count(int(occupancy_lookup.get(n_samples, 0))),
            _format_report_fraction(
                _safe_fraction(int(occupancy_lookup.get(n_samples, 0)), graph_summary.loci_written)
            ),
        ]
        for n_samples in range(qc_summary.selected_sample_count + 1)
    ]
    runtime_rows = [
        ("cores", _format_report_count(cores)),
        ("vsearch_threads_per_job", _format_report_count(threads)),
        ("vsearch_worker_processes", _format_report_count(workers)),
        ("across_vsearch_threads", _format_report_count(cores)),
        ("mafft_threads_per_job", _format_report_count(alignment_summary.mafft_threads_per_job)),
        (
            "mafft_worker_processes",
            _format_report_count(alignment_summary.mafft_worker_processes),
        ),
        ("alignment_mode", alignment_summary.alignment_mode),
        (
            "mafft_timeout_seconds",
            _format_report_count(alignment_summary.mafft_timeout_seconds),
        ),
        ("duplicated_component_reconciliation", "same-sample graph"),
        ("cluster_spacer_mode", "stripped"),
        ("output_spacer_length", _format_report_count(alignment_summary.output_spacer_length)),
        ("keep_intermediates", str(bool(keep_intermediates))),
        ("vsearch_binary", BIN_VSEARCH),
        ("mafft_binary", BIN_MAFFT),
        ("workdir", str(workdir)),
    ]
    output_rows = [
        ("reference", str(outputs["reference"])),
        ("mapping", str(outputs["mapping"])),
        ("loci_stats", str(outputs["loci_stats"])),
        ("sample_graph_summary", str(outputs["sample_graph_summary"])),
        ("run_stats", str(outputs["run_stats"])),
        ("audit_dir", str(outputs["audit_dir"])),
        (
            "intermediates",
            "retained in workdir" if keep_intermediates else "cleaned on success",
        ),
    ]

    lines: list[str] = []
    _append_report_key_value_section(lines, "Inputs", inputs_rows)
    _append_report_key_value_section(lines, "Clustering Parameters", clustering_rows)
    _append_report_key_value_section(lines, "Denovo Summary", summary_rows)
    _append_report_key_value_section(lines, "Locus QC", locus_qc_rows)
    _append_report_key_value_section(lines, "Component QC", component_qc_rows)
    _append_report_table_section(lines, "Component Node Summary", component_node_headers, component_node_rows)
    _append_report_table_section(lines, "Selected Sample Summary", selected_sample_headers, selected_sample_rows)
    _append_report_table_section(lines, "Locus Occupancy", occupancy_headers, occupancy_rows)
    _append_report_key_value_section(lines, "Runtime", runtime_rows)
    _append_report_key_value_section(lines, "Outputs", output_rows)

    outpath.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    logger.info(f"wrote denovo run stats to {outpath}")


def _write_stripped_clustering_fasta(
    out_fasta: Path,
    input_fastas: list[Path],
) -> dict[str, tuple[str, int]]:
    """Write one FASTA with joined spacers stripped before within-sample clustering."""
    seed_to_meta: dict[str, tuple[str, int]] = {}
    with open(out_fasta, "wt", encoding="utf-8") as out:
        for fasta_path in input_fastas:
            if not fasta_path.exists():
                continue
            header: str | None = None
            with open(fasta_path, "rt", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        header = line[1:].strip()
                        out.write(f">{header}\n")
                        continue
                    if header is None:
                        continue
                    record_type = infer_record_type(header)
                    cluster_sequence, arm_boundary = get_arm_boundary(line)
                    seed_to_meta[header] = (record_type, int(arm_boundary))
                    out.write(f"{cluster_sequence}\n")
    return seed_to_meta


def vsearch_pairs(
    sname: str,
    r1: Path,
    r2: Path | None,
    outdir: Path,
    min_derep_size: int,
    min_merge_overlap: int,
    min_length: int,
    max_merge_diffs: int,
    allow_reverse_complement: bool,
    within_similarity: float,
    threads: int,
    keep_intermediates: bool,
    paired: bool = False,
) -> None:
    """Run the within-sample vsearch workflow for one sample."""
    unmerged_r1 = outdir / f"{sname}.unmerged_R1.fq"
    unmerged_r2 = outdir / f"{sname}.unmerged_R2.fq"
    merged = outdir / f"{sname}.merged.fa"
    joined = outdir / f"{sname}.joined.fa"
    cluster_fa = outdir / f"{sname}.cluster.fa"
    derep_unsorted = outdir / f"{sname}.derep.fa"
    derep = outdir / f"{sname}.derep.sizesorted.fa"
    consensus = outdir / f"{sname}.consensus.fa"
    clusters = outdir / f"{sname}.clusters.tsv"

    # get PE or SE data ready for dereplication
    cluster_inputs: list[Path] = []
    if paired:
        if r2 is None:
            raise IPyradError(f"Missing R2 FASTQ for paired sample: {sname}")
        cmd1 = [
            BIN_VSEARCH, "--fastq_mergepairs", str(r1),
            "--reverse", str(r2),
            "--fastq_minovlen", str(min_merge_overlap),
            "--fastq_maxdiffs", str(max_merge_diffs),
            "--fastq_minlen", str(min_length),
            "--fastq_allowmergestagger",
            "--fasta_width", "0",
            "--fastqout_notmerged_fwd", str(unmerged_r1),
            "--fastqout_notmerged_rev", str(unmerged_r2),
            "--relabel", f"{sname};M",
            "--fastaout", str(merged),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])

        cmd1 = [
            BIN_VSEARCH, "--fastq_join", str(unmerged_r1),
            "--reverse", str(unmerged_r2),
            "--join_padgap", "N" * 24,
            "--join_padgapq", "I" * 24,
            "--fasta_width", "0",
            "--relabel", f"{sname};J",
            "--fastaout", str(joined),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])
        cluster_inputs.extend([joined, merged])
    else:
        cmd1 = [
            BIN_VSEARCH,
            "--fastx_subsample", str(r1),
            "--sample_pct", "100",
            "--relabel", f"{sname};S",
            "--fastaout", str(joined),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])
        cluster_inputs.append(joined)

    seed_to_meta = _write_stripped_clustering_fasta(cluster_fa, cluster_inputs)

    # run dereplication, then sort dereplicated records by length
    cmd1 = [
        BIN_VSEARCH,
        "--fastx_uniques", str(cluster_fa),
        "--minuniquesize", str(min_derep_size),
        "--strand", "both" if allow_reverse_complement else "plus",
        "--fasta_width", "0",
        "--sizeout",
        "--relabel_keep",
        "--fastaout", str(derep_unsorted),
    ]
    logger.debug(" ".join(cmd1))
    run_pipeline([cmd1])
    if not keep_intermediates:
        _cleanup_post_derep_work_files(sname, outdir, paired=paired)

    cmd1 = [
        BIN_VSEARCH,
        "--sortbylength", str(derep_unsorted),
        "--sizein",
        "--sizeout",
        "--fasta_width", "0",
        "--output", str(derep),
    ]
    logger.debug(" ".join(cmd1))
    run_pipeline([cmd1])

    # run within-sample clustering
    cmd1 = [
        BIN_VSEARCH,
        "--cluster_fast", str(derep),
        "--id", str(within_similarity),
        "--strand", "both" if allow_reverse_complement else "plus",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        "--query_cov", "0.75",
        "--fasta_width", "0",
        "--qmask", "none",
        "--consout", str(consensus),
        "--uc", str(clusters),
        "--threads", str(threads),
    ]
    logger.debug(" ".join(cmd1))
    run_pipeline([cmd1])

    # write a '{sname}.summary.tsv' and remove intermediate files
    build_sample_summary(sname, outdir, seed_to_meta=seed_to_meta)
    if not keep_intermediates:
        _cleanup_post_summary_work_files(sname, outdir)


def _write_cluster_sequence_fasta(
    summary_tsv: Path,
    out_fasta: Path,
) -> Path:
    """Write a clustering FASTA from spacer-stripped summary sequences."""
    with open(summary_tsv, "rt", encoding="utf-8", newline="") as infile, open(
        out_fasta,
        "wt",
        encoding="utf-8",
    ) as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        required = {"seed", "cluster_sequence"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            joined = ", ".join(sorted(required))
            raise RuntimeError(f"concat.summary.tsv is missing required columns: {joined}")
        for row in reader:
            seed = str(row["seed"])
            seq = str(row["cluster_sequence"]).upper()
            outfile.write(f">{seed}\n{seq}\n")
    return out_fasta


def vsearch_cluster_across(
    outdir: Path,
    summary_tsv: Path,
    across_similarity: float,
    threads: int,
) -> None:
    """Cluster all sample-level consensus sequences across samples."""

    # write concatenated consensus seqs to
    consensus_concat = outdir / "consensus.concat.fa"
    _write_cluster_sequence_fasta(summary_tsv, consensus_concat)

    # tabular cluster hits will be written here
    cluster_table = outdir / "global_hits.uc.tsv"

    cmd1 = [
        BIN_VSEARCH,
        "--usearch_global", str(consensus_concat),
        "--db", str(consensus_concat),
        "--id", str(across_similarity),
        "--userout", str(cluster_table),
        "--userfields", "query+target+id+qstrand+qcov+ql+tl",
        "--maxaccepts", "0",
        "--maxrejects", "0",
        "--query_cov", "0.75",
        "--self",
        "--qmask", "none",
        "--notmatched", "/dev/null",
        "--fasta_width", "0",
        "--threads", str(threads),
    ]
    _run_vsearch_with_progress(cmd1, message="Across-sample clustering")


def run_denovo(
    fastqs: list[Path],
    outdir: Path,
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    delim_str: str | None,
    delim_idx: int,
    allow_reverse_complement: bool,
    cores: int,
    threads: int,
    no_alignment: bool,
    force: bool,
    imap: Path | None,
    use_all_samples: bool,
    keep_intermediates: bool,
    log_level: str,
) -> None:
    """Run the denovo reference construction workflow."""
    # check cli arg values
    _validate_runtime_args(
        within_similarity=within_similarity,
        across_similarity=across_similarity,
        min_derep_size=min_derep_size,
        min_length=min_length,
        min_merge_overlap=min_merge_overlap,
        max_merge_diffs=max_merge_diffs,
        cores=cores,
        threads=threads,
        delim_idx=delim_idx,
    )

    # create outdir and clear old files if --force
    outdir = outdir.expanduser().absolute()
    workdir, outputs = _prepare_output_paths(outdir, force=force)

    # parse input fastqs/pairs, subselect, and report
    full_fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)
    _validate_fastq_layout(full_fastq_dict)
    fastq_dict, selection_mode = _select_denovo_samples(
        full_fastq_dict,
        imap_path=imap,
        use_all_samples=use_all_samples,
    )
    _log_sample_selection(full_fastq_dict, fastq_dict, selection_mode)
    is_paired = _validate_fastq_layout(fastq_dict)
    _validate_required_binaries()
    workers = max(1, cores // threads)

    # I REMOVED THIS ON PURPOSE. DO NOT ADD IT BACK IN.
    # vsearch_path = _validate_required_binary(BIN_VSEARCH, "vsearch")
    # logger.info("paired mode: {}", "paired-end" if is_paired else "single-end")
    # logger.info("using vsearch binary: {}", vsearch_path)
    # logger.info("using mafft binary: {}", mafft_path)

    # perform within-samples operations to get within-sample consensus loci
    msg = "Joining/merging pairs, dereplicating, and clustering" if is_paired else "Dereplicating and clustering"
    jobs: dict[str, tuple[object, dict[str, object]]] = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            sname=sname,
            r1=fastq_tuple[0],
            r2=fastq_tuple[1],
            outdir=workdir,
            min_derep_size=min_derep_size,
            min_length=min_length,
            min_merge_overlap=min_merge_overlap,
            max_merge_diffs=max_merge_diffs,
            allow_reverse_complement=allow_reverse_complement,
            within_similarity=within_similarity,
            threads=threads,
            keep_intermediates=keep_intermediates,
            paired=is_paired,
        )
        jobs[sname] = (vsearch_pairs, kwargs)
    run_with_pool(jobs, log_level, workers, msg=msg)

    # write '.concat.summary.tsv' (this is a tmp file for debugging)
    concat_summaries(workdir)

    # perform across-sample clustering and write a summary
    logger.info("Clustering consensus sequences across samples")
    vsearch_cluster_across(
        outdir=workdir,
        summary_tsv=workdir / "concat.summary.tsv",
        across_similarity=across_similarity,
        threads=cores,
    )

    # perform graph splitting on the across samples clusters
    logger.info("Splitting global clusters and writing locus tables")
    graph_summary = make_global_tables(
        workdir,
        cores=cores,
        log_level=log_level,
        within_similarity=within_similarity,
    )

    # optionally perform alignment on the across-sample clusters
    if no_alignment:
        logger.info("Selecting longest locus representatives and writing denovo reference")
    else:
        logger.info("Aligning locus consensuses and writing denovo reference")
    alignment_summary = write_ordered_consensus_stream_to_file(
        mapping_tsv=outputs["mapping"],
        summary_tsv=workdir / "concat.summary.tsv",
        out_fa=outputs["reference"],
        mafft_binary=BIN_MAFFT,
        cores=cores,
        alignment_mode="none" if no_alignment else "mafft",
    )

    qc_summary = _collect_denovo_qc(
        selected_fastq_dict=fastq_dict,
        total_input_sample_count=len(full_fastq_dict),
        graph_summary=graph_summary,
        outputs=outputs,
    )
    _write_denovo_sample_graph_summary(
        outputs["sample_graph_summary"],
        mapping_tsv=outputs["mapping"],
        loci_stats_tsv=outputs["loci_stats"],
        sample_names=sorted(fastq_dict),
    )

    # collect and write stats on the denovo assembly
    _write_denovo_stats(
        outputs["run_stats"],
        all_fastq_dict=full_fastq_dict,
        selected_fastq_dict=fastq_dict,
        selection_mode=selection_mode,
        paired=is_paired,
        within_similarity=within_similarity,
        across_similarity=across_similarity,
        min_derep_size=min_derep_size,
        min_length=min_length,
        min_merge_overlap=min_merge_overlap,
        max_merge_diffs=max_merge_diffs,
        allow_reverse_complement=allow_reverse_complement,
        cores=cores,
        threads=threads,
        workers=workers,
        alignment_summary=alignment_summary,
        keep_intermediates=keep_intermediates,
        workdir=workdir,
        graph_summary=graph_summary,
        qc_summary=qc_summary,
        outputs=outputs,
    )

    if not keep_intermediates:
        _cleanup_denovo_workdir_stage_files(workdir)
        shutil.rmtree(workdir)


if __name__ == "__main__":
    pass
