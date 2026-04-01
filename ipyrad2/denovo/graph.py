#!/usr/bin/env python

"""Build, contract, and split denovo across-sample graphs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from collections import Counter
from itertools import combinations
from statistics import median

from loguru import logger
import pandas as pd

from .align import mafft_align_one
from .common import infer_record_type, split_joined_sequence
from .splitters import get_splitter
from .splitters.common import connected_components, sort_components
from ..utils.parallel import run_with_pool_iter


RECONCILE_QCOV = 0.75
RECONCILE_MIN_BOUNDARY_SPAN = 12


@dataclass(frozen=True)
class ContractedComponent:
    """One post-contraction graph component ready for splitter dispatch."""

    component_id: int
    duplicated_component: bool
    aligned_for_reconciliation: bool
    arm_boundary_found: bool
    used_reconciliation: bool
    node_members: dict[str, tuple[str, ...]]
    node_samples: dict[str, frozenset[str]]
    node_order: dict[str, int]
    node_modes: dict[str, str]
    edges: dict[tuple[str, str], tuple[float, float]]


@dataclass(frozen=True)
class SummaryRecord:
    """Minimal summary fields needed for graph contraction and table writing."""

    seed: str
    sample: str
    n_reads: int
    n_unique: int
    length: int
    cluster_length: int
    merged: bool
    record_type: str
    cluster_id: int
    cluster_sequence: str
    left_arm: str
    right_arm: str


@dataclass(frozen=True)
class ComponentPart:
    """One final split result returned from one component worker."""

    component_id: int
    subcomponent_id: int
    mapping_rows: tuple[dict[str, object], ...]
    stats_row: dict[str, object]


@dataclass(frozen=True)
class ComponentResult:
    """One processed component plus audit artifacts."""

    component_id: int
    parts: tuple[ComponentPart, ...]
    audit_summary: dict[str, object]
    audit_rows: tuple[dict[str, object], ...]
    audit_fasta: tuple[tuple[str, str], ...]


def get_edges_dict(outdir: Path) -> dict[tuple[str, str], tuple[float, float]]:
    """Return graph edges mapped to their strongest `(pid, qcov)` score."""
    uc_path = outdir / "global_hits.uc.tsv"
    edges: dict[tuple[str, str], tuple[float, float]] = {}
    with open(uc_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            query, target, pid, _qstrand, qcov, _qlen, _tlen = line.rstrip().split("\t")
            if query == target:
                continue
            edge = tuple(sorted((query, target)))
            score = (float(pid) / 100.0, float(qcov) / 100.0)
            prev = edges.get(edge)
            if prev is None or score > prev:
                edges[edge] = score
    return edges


def get_summary_df(outdir: Path) -> pd.DataFrame:
    """Load the concatenated denovo summary table."""
    tsv = outdir / "concat.summary.tsv"
    df = pd.read_csv(tsv, sep="\t").copy()

    if "record_type" not in df.columns:
        df["record_type"] = df["seed"].map(infer_record_type)
    if "cluster_sequence" not in df.columns or "left_arm" not in df.columns or "right_arm" not in df.columns:
        cluster_values = df["consensus"].map(split_joined_sequence)
        df["cluster_sequence"] = [value[0] for value in cluster_values]
        df["left_arm"] = [value[1] for value in cluster_values]
        df["right_arm"] = [value[2] for value in cluster_values]
    if "cluster_length" not in df.columns:
        df["cluster_length"] = df["cluster_sequence"].map(len)
    return df


def _get_summary_records(df: pd.DataFrame) -> dict[str, SummaryRecord]:
    """Return minimal per-seed summary records for graph processing."""
    return {
        str(row["seed"]): SummaryRecord(
            seed=str(row["seed"]),
            sample=str(row["sample"]),
            n_reads=int(row["n_reads"]),
            n_unique=int(row["n_unique"]),
            length=int(row["length"]),
            cluster_length=int(row["cluster_length"]),
            merged=bool(row["merged"]),
            record_type=str(row["record_type"]),
            cluster_id=int(row["cluster_id"]),
            cluster_sequence=str(row["cluster_sequence"]).upper(),
            left_arm=str(row["left_arm"]).upper(),
            right_arm=str(row["right_arm"]).upper(),
        )
        for row in df.to_dict(orient="records")
    }


def _seed_order(df: pd.DataFrame) -> dict[str, int]:
    """Return stable seed ordering from concat.summary.tsv row order."""
    return {str(seed): idx for idx, seed in enumerate(df["seed"].tolist(), start=1)}


def _component_has_duplicate_samples(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
) -> bool:
    """Return True when one component contains repeated samples."""
    counts = Counter(summary_records[node].sample for node in component)
    return any(count > 1 for count in counts.values())


def _is_joined_record(record: SummaryRecord) -> bool:
    """Return True when one summary record is a joined pair with both arms."""
    return record.record_type == "joined" and bool(record.right_arm)


def _component_needs_reconciliation(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
) -> bool:
    """Return True when one duplicated component should enter reconciliation."""
    return _component_has_duplicate_samples(component, summary_records)


def _aligned_breakpoint_column(aligned_seq: str, left_len: int) -> int:
    """Return the aligned column index immediately after the left arm."""
    if left_len <= 0:
        return 0
    seen = 0
    for idx, char in enumerate(aligned_seq):
        if char != "-":
            seen += 1
            if seen >= left_len:
                return idx + 1
    return len(aligned_seq)


def _infer_component_arm_boundary(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    aligned_map: dict[str, str],
) -> int | None:
    """Return a coherent aligned arm boundary for joined records, or None."""
    boundaries = [
        _aligned_breakpoint_column(aligned_map[node], len(summary_records[node].left_arm))
        for node in component
        if _is_joined_record(summary_records[node])
    ]
    if not boundaries:
        return None
    boundary = int(round(median(boundaries)))
    align_len = len(next(iter(aligned_map.values()))) if aligned_map else 0
    tolerance = max(
        RECONCILE_MIN_BOUNDARY_SPAN,
        int(round(align_len * 0.05)),
    )
    if any(abs(value - boundary) > tolerance for value in boundaries):
        return None
    return boundary


def _align_component_sequences(
    component_id: int,
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    mafft_binary: str,
) -> dict[str, str]:
    """Return aligned cluster sequences for one component."""
    record = [
        (core, summary_records[core].cluster_sequence)
        for core in sorted(component, key=seed_order.get)
    ]
    return dict(
        mafft_align_one(
            record,
            mafft_binary=mafft_binary,
            threads=1,
            locus_id=component_id,
        )
    )


def _informative_length(seq: str) -> int:
    """Return the number of informative ACGT sites in one stripped sequence."""
    return sum(char.upper() in "ACGT" for char in seq)


def _aligned_pair_score(
    left: str,
    right: str,
    summary_records: dict[str, SummaryRecord],
    aligned_map: dict[str, str],
) -> tuple[float, float] | None:
    """Return one symmetric `(identity, overlap)` score from an aligned pair."""
    matches = 0
    informative = 0
    left_seq = aligned_map[left]
    right_seq = aligned_map[right]
    for left_char, right_char in zip(left_seq, right_seq):
        left_base = left_char.upper()
        right_base = right_char.upper()
        if left_base not in "ACGT" or right_base not in "ACGT":
            continue
        informative += 1
        if left_base == right_base:
            matches += 1
    if informative == 0:
        return None
    left_len = _informative_length(summary_records[left].cluster_sequence)
    right_len = _informative_length(summary_records[right].cluster_sequence)
    denom = max(left_len, right_len, 1)
    return matches / informative, informative / denom


def _passes_reconcile_threshold(
    score: tuple[float, float] | None,
    across_similarity: float,
) -> bool:
    """Return True when one aligned pair passes reconciliation thresholds."""
    if score is None:
        return False
    identity, overlap = score
    return identity >= across_similarity and overlap >= RECONCILE_QCOV


def _get_reconcile_scores(
    component_id: int,
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    across_similarity: float,
    mafft_binary: str,
) -> tuple[bool, bool, dict[tuple[str, str], tuple[float, float]]]:
    """Return aligned compatibility scores for one duplicated component."""
    if not _component_needs_reconciliation(component, summary_records):
        return False, False, {}

    aligned_map = _align_component_sequences(
        component_id=component_id,
        component=component,
        summary_records=summary_records,
        seed_order=seed_order,
        mafft_binary=mafft_binary,
    )
    arm_boundary_found = (
        _infer_component_arm_boundary(component, summary_records, aligned_map) is not None
    )

    compatible_scores: dict[tuple[str, str], tuple[float, float]] = {}
    for left, right in combinations(sorted(component, key=seed_order.get), 2):
        score = _aligned_pair_score(left, right, summary_records, aligned_map)
        if not _passes_reconcile_threshold(score, across_similarity):
            continue
        edge = tuple(sorted((left, right)))
        compatible_scores[edge] = score
    return True, arm_boundary_found, compatible_scores


def _same_sample_groups(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    compatible_scores: dict[tuple[str, str], tuple[float, float]],
) -> list[tuple[str, tuple[str, ...], str]]:
    """Return same-sample technical-duplicate groups eligible for contraction."""
    out: list[tuple[str, tuple[str, ...], str]] = []
    samples = sorted({summary_records[node].sample for node in component})
    for sample in samples:
        sample_nodes = [
            node
            for node in component
            if summary_records[node].sample == sample
        ]
        if len(sample_nodes) <= 1:
            continue
        sample_edges: dict[tuple[str, str], tuple[float, float]] = {}
        for left_idx, left in enumerate(sample_nodes):
            for right in sample_nodes[left_idx + 1:]:
                edge = tuple(sorted((left, right)))
                score = compatible_scores.get(edge)
                if score is not None:
                    sample_edges[edge] = score
        for subcomponent in sort_components(
            connected_components(sample_nodes, sample_edges),
            seed_order,
        ):
            if len(subcomponent) <= 1:
                continue
            ordered = tuple(sorted(subcomponent, key=seed_order.get))
            mode = (
                "mixed"
                if any(summary_records[node].record_type != "joined" for node in ordered)
                else "joined_only"
            )
            out.append((sample, ordered, mode))
    return out


def _build_contracted_component(
    component_id: int,
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    edges: dict[tuple[str, str], tuple[float, float]],
    seed_order: dict[str, int],
    across_similarity: float,
    mafft_binary: str,
) -> ContractedComponent:
    """Contract technical duplicates inside one component."""
    duplicated_component = _component_needs_reconciliation(component, summary_records)
    aligned_for_reconciliation, arm_boundary_found, compatible_scores = _get_reconcile_scores(
        component_id=component_id,
        component=component,
        summary_records=summary_records,
        seed_order=seed_order,
        across_similarity=across_similarity,
        mafft_binary=mafft_binary,
    )

    groups: list[tuple[str, tuple[str, ...], str]] = []
    used: set[str] = set()

    if aligned_for_reconciliation:
        for sample, members, mode in _same_sample_groups(
            component,
            summary_records,
            seed_order,
            compatible_scores,
        ):
            groups.append((sample, members, mode))
            used.update(members)

    used_reconciliation = bool(groups)

    for core in sorted(component, key=seed_order.get):
        if core in used:
            continue
        sample = summary_records[core].sample
        groups.append((sample, (core,), "none"))

    groups.sort(key=lambda item: min(seed_order[core] for core in item[1]))
    node_members: dict[str, tuple[str, ...]] = {}
    node_samples: dict[str, frozenset[str]] = {}
    node_order: dict[str, int] = {}
    node_modes: dict[str, str] = {}

    for local_idx, (sample, members, mode) in enumerate(groups, start=1):
        node_id = f"contract_{component_id}_{local_idx}"
        ordered_members = tuple(sorted(members, key=seed_order.get))
        node_members[node_id] = ordered_members
        node_samples[node_id] = frozenset({sample})
        node_order[node_id] = min(seed_order[core] for core in ordered_members)
        node_modes[node_id] = mode

    node_ids = sorted(node_members, key=node_order.get)
    contracted_edges: dict[tuple[str, str], tuple[float, float]] = {}
    for left_idx, left in enumerate(node_ids):
        left_members = node_members[left]
        for right in node_ids[left_idx + 1:]:
            best: tuple[float, float] | None = None
            for left_core in left_members:
                for right_core in node_members[right]:
                    edge = tuple(sorted((left_core, right_core)))
                    score = compatible_scores.get(edge) if aligned_for_reconciliation else edges.get(edge)
                    if score is not None and (best is None or score > best):
                        best = score
            if best is not None:
                contracted_edges[(left, right)] = best

    return ContractedComponent(
        component_id=component_id,
        duplicated_component=duplicated_component,
        aligned_for_reconciliation=aligned_for_reconciliation,
        arm_boundary_found=arm_boundary_found,
        used_reconciliation=used_reconciliation,
        node_members=node_members,
        node_samples=node_samples,
        node_order=node_order,
        node_modes=node_modes,
        edges=contracted_edges,
    )


def _validate_sample_uniqueness(
    nodes: Iterable[str],
    node_samples: dict[str, frozenset[str]],
) -> None:
    """Raise if a final split still contains duplicate samples."""
    samples: list[str] = []
    for node in nodes:
        samples.extend(node_samples[node])
    if len(samples) != len(set(samples)):
        raise RuntimeError("denovo splitter emitted a final locus with duplicate samples")


def _bucket_component_edges(
    components: list[set[str]],
    edges: dict[tuple[str, str], tuple[float, float]],
) -> dict[int, dict[tuple[str, str], tuple[float, float]]]:
    """Return per-component edge dictionaries in one pass over the global edges."""
    node_to_component: dict[str, int] = {}
    for component_id, component in enumerate(components, start=1):
        for node in component:
            node_to_component[node] = component_id

    component_edges = {component_id: {} for component_id in range(1, len(components) + 1)}
    for edge, score in edges.items():
        component_id = node_to_component[edge[0]]
        component_edges[component_id][edge] = score
    return component_edges


def _std_or_nan(values: list[float | int]) -> float:
    """Return pandas-compatible sample std, including NaN for one value."""
    return float(pd.Series(values, dtype="float64").std())


def _part_reconcile_mode(
    *,
    used_reconciliation: bool,
    has_joined: bool,
    has_merged: bool,
) -> str:
    """Return the reconcile mode label for one final locus part."""
    if not used_reconciliation:
        return "none"
    if has_joined and has_merged:
        return "mixed"
    if has_joined:
        return "joined_only"
    return "none"


def _part_output_form(
    *,
    mapping_rows: list[dict[str, object]],
    reconcile_mode: str,
    arm_boundary_found: bool,
) -> str:
    """Return the output form for one final locus part."""
    if all(row["record_type"] == "joined" for row in mapping_rows):
        return "spaced"
    if reconcile_mode in {"joined_only", "mixed"} and arm_boundary_found:
        return "spaced"
    return "stripped"


def _process_component(
    component_id: int,
    component: tuple[str, ...],
    summary_records: dict[str, SummaryRecord],
    component_edges: dict[tuple[str, str], tuple[float, float]],
    seed_order: dict[str, int],
    graph_splitter: str,
    across_similarity: float,
    mafft_binary: str,
) -> tuple[int, ComponentResult]:
    """Contract, split, and expand one connected component."""
    splitter = get_splitter(graph_splitter)
    contracted = _build_contracted_component(
        component_id,
        set(component),
        summary_records,
        component_edges,
        seed_order,
        across_similarity=across_similarity,
        mafft_binary=mafft_binary,
    )
    subcomponents = splitter(
        set(contracted.node_members),
        contracted.edges,
        contracted.node_samples,
        contracted.node_order,
    )

    parts: list[ComponentPart] = []
    audit_rows: list[dict[str, object]] = []
    for subcomponent_id, nodes in enumerate(subcomponents, start=1):
        _validate_sample_uniqueness(nodes, contracted.node_samples)
        ordered_nodes = tuple(sorted(nodes, key=contracted.node_order.get))
        ordered_cores: list[tuple[str, str]] = []
        for node in ordered_nodes:
            for core in contracted.node_members[node]:
                ordered_cores.append((node, core))
        ordered_cores.sort(key=lambda item: seed_order[item[1]])

        mapping_rows: list[dict[str, object]] = []
        n_reads_values: list[int] = []
        length_values: list[int] = []
        merged_values: list[int] = []
        samples: list[str] = []
        contract_group_counts: Counter[str] = Counter()
        node_modes: dict[str, str] = {}

        for node in ordered_nodes:
            samples.append(next(iter(contracted.node_samples[node])))
            node_modes[node] = contracted.node_modes[node]

        for contract_group, core in ordered_cores:
            info = summary_records[core]
            contract_group_counts[contract_group] += 1
            n_reads_values.append(info.n_reads)
            length_values.append(info.length)
            merged_values.append(int(info.merged))
            mapping_rows.append({
                "component_id": int(component_id),
                "subcomponent_id": int(subcomponent_id),
                "contract_group": contract_group,
                "sample": info.sample,
                "n_reads": info.n_reads,
                "n_unique": info.n_unique,
                "length": info.length,
                "cluster_length": info.cluster_length,
                "merged": int(info.merged),
                "record_type": info.record_type,
                "cluster_id": info.cluster_id,
                "core": core,
            })

        has_joined = any(row["record_type"] == "joined" for row in mapping_rows)
        has_merged = any(row["record_type"] == "merged" for row in mapping_rows)
        reconcile_mode = _part_reconcile_mode(
            used_reconciliation=contracted.used_reconciliation,
            has_joined=has_joined,
            has_merged=has_merged,
        )
        output_form = _part_output_form(
            mapping_rows=mapping_rows,
            reconcile_mode=reconcile_mode,
            arm_boundary_found=contracted.arm_boundary_found,
        )
        n_reconciled_groups = (
            len(set(contract_group_counts))
            if contracted.used_reconciliation
            else 0
        )
        n_mixed_records = len(mapping_rows) if reconcile_mode == "mixed" else 0
        for row in mapping_rows:
            row["reconcile_mode"] = reconcile_mode
            row["reconciled_group"] = row["contract_group"] if contracted.used_reconciliation else ""
            row["output_form"] = output_form
            audit_rows.append({
                "component_id": int(component_id),
                "subcomponent_id": int(subcomponent_id),
                "locus_name": f"locus_{component_id}_{subcomponent_id}",
                "sample": row["sample"],
                "core": row["core"],
                "record_type": row["record_type"],
                "cluster_id": row["cluster_id"],
                "contract_group": row["contract_group"],
                "reconcile_mode": reconcile_mode,
                "reconciled_group": row["reconciled_group"],
                "output_form": output_form,
            })

        stats_row = {
            "component_id": int(component_id),
            "subcomponent_id": int(subcomponent_id),
            "n_samples": len(set(samples)),
            "n_cores": len(mapping_rows),
            "n_contracted_groups": sum(1 for count in contract_group_counts.values() if count > 1),
            "n_reconciled_groups": int(n_reconciled_groups),
            "n_mixed_records": int(n_mixed_records),
            "n_reads_sum": int(sum(n_reads_values)),
            "n_reads_mean": float(sum(n_reads_values) / len(n_reads_values)),
            "n_reads_std": _std_or_nan(n_reads_values),
            "length_mean": float(sum(length_values) / len(length_values)),
            "length_std": _std_or_nan(length_values),
            "merged_freq": float(sum(merged_values) / len(merged_values)),
            "duplicated_component": bool(contracted.duplicated_component),
            "aligned_for_reconciliation": bool(contracted.aligned_for_reconciliation),
            "arm_boundary_found": bool(contracted.arm_boundary_found),
            "used_reconciliation": bool(contracted.used_reconciliation),
            "reconcile_mode": reconcile_mode,
            "used_joined_only_reconciliation": bool(reconcile_mode == "joined_only"),
            "used_mixed_reconciliation": bool(reconcile_mode == "mixed"),
            "output_form": output_form,
            "samples": ",".join(samples),
        }
        parts.append(
            ComponentPart(
                component_id=component_id,
                subcomponent_id=subcomponent_id,
                mapping_rows=tuple(mapping_rows),
                stats_row=stats_row,
            )
        )

    counts = Counter(summary_records[node].sample for node in component)
    component_has_joined = any(summary_records[node].record_type == "joined" for node in component)
    component_has_merged = any(summary_records[node].record_type == "merged" for node in component)
    component_mode = _part_reconcile_mode(
        used_reconciliation=contracted.used_reconciliation,
        has_joined=component_has_joined,
        has_merged=component_has_merged,
    )
    audit_summary = {
        "component_id": int(component_id),
        "n_input_nodes": int(len(component)),
        "n_input_samples": int(len(counts)),
        "n_duplicate_samples": int(sum(1 for value in counts.values() if value > 1)),
        "has_joined": bool(component_has_joined),
        "has_merged": bool(component_has_merged),
        "aligned_for_reconciliation": bool(contracted.aligned_for_reconciliation),
        "arm_boundary_found": bool(contracted.arm_boundary_found),
        "used_reconciliation": bool(contracted.used_reconciliation),
        "reconcile_mode": component_mode,
        "graph_splitter": graph_splitter,
        "n_final_loci": int(len(parts)),
    }
    audit_fasta = tuple(
        (core, summary_records[core].cluster_sequence)
        for core in sorted(component, key=seed_order.get)
    )
    return component_id, ComponentResult(
        component_id=component_id,
        parts=tuple(parts),
        audit_summary=audit_summary,
        audit_rows=tuple(audit_rows),
        audit_fasta=audit_fasta,
    )


def _write_component_audits(
    audit_dir: Path,
    results: list[ComponentResult],
) -> None:
    """Write compact audit outputs for duplicated components."""
    if audit_dir.exists():
        for path in audit_dir.iterdir():
            if path.is_dir():
                continue
            path.unlink()
    else:
        audit_dir.mkdir(parents=True, exist_ok=True)

    summaries = [
        result.audit_summary
        for result in results
        if result.audit_summary["n_duplicate_samples"] > 0
    ]
    summary_df = pd.DataFrame(summaries)
    summary_path = audit_dir / "components.summary.tsv"
    if summary_df.empty:
        summary_df = pd.DataFrame(
            columns=[
                "component_id",
                "n_input_nodes",
                "n_input_samples",
                "n_duplicate_samples",
                "has_joined",
                "has_merged",
                "aligned_for_reconciliation",
                "arm_boundary_found",
                "used_reconciliation",
                "reconcile_mode",
                "graph_splitter",
                "n_final_loci",
            ]
        )
    summary_df.to_csv(summary_path, sep="\t", index=False)

    for result in results:
        if result.audit_summary["n_duplicate_samples"] <= 0:
            continue
        members_df = pd.DataFrame(result.audit_rows)
        members_df.to_csv(
            audit_dir / f"component_{result.component_id}.members.tsv",
            sep="\t",
            index=False,
        )
        fasta_path = audit_dir / f"component_{result.component_id}.fa"
        with open(fasta_path, "wt", encoding="utf-8") as out:
            for core, seq in result.audit_fasta:
                out.write(f">{core}\n{seq}\n")


def make_global_tables(
    outdir: Path,
    graph_splitter: str = "constrained",
    *,
    cores: int = 1,
    log_level: str = "INFO",
    across_similarity: float = 0.85,
    mafft_binary: str = "mafft",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build locus mapping/stats tables from the across-sample denovo graph."""
    if cores < 1:
        raise ValueError("cores must be >= 1")
    out_mapping_tsv = outdir.parent / "loci.mapping.tsv"
    out_stats_tsv = outdir.parent / "loci.stats.tsv"
    audit_dir = outdir.parent / "denovo.audit"

    df = get_summary_df(outdir)
    summary_records = _get_summary_records(df)
    seed_order = _seed_order(df)
    edges = get_edges_dict(outdir)

    components = sort_components(connected_components(df["seed"].tolist(), edges), seed_order)
    component_edges = _bucket_component_edges(components, edges)
    final_parts: list[ComponentPart] = []
    component_results: list[ComponentResult] = []

    if components:
        jobs = {
            component_id: (
                _process_component,
                {
                    "component_id": component_id,
                    "component": tuple(sorted(component, key=seed_order.get)),
                    "summary_records": {seed: summary_records[seed] for seed in component},
                    "component_edges": component_edges[component_id],
                    "seed_order": {seed: seed_order[seed] for seed in component},
                    "graph_splitter": graph_splitter,
                    "across_similarity": across_similarity,
                    "mafft_binary": mafft_binary,
                },
            )
            for component_id, component in enumerate(components, start=1)
        }
        if cores == 1 or len(jobs) == 1:
            for component_id in sorted(jobs):
                _key, result = _process_component(**jobs[component_id][1])
                component_results.append(result)
                final_parts.extend(result.parts)
        else:
            for _key, (_component_id, result) in run_with_pool_iter(
                jobs.items(),
                log_level=log_level,
                max_workers=min(cores, len(jobs)),
                msg="Splitting global clusters",
                njobs=len(jobs),
            ):
                component_results.append(result)
                final_parts.extend(result.parts)
        component_results.sort(key=lambda result: result.component_id)
        final_parts.sort(key=lambda part: (part.component_id, part.subcomponent_id))

    logger.info(
        f"split {len(components)} clusters into {len(final_parts)} non-duplicated subclusters "
        f"using graph splitter '{graph_splitter}'"
    )

    mapping_rows: list[dict[str, object]] = []
    stats_rows: list[dict[str, object]] = []
    for locus_idx, part in enumerate(final_parts, start=1):
        locus_name = f"locus_{part.component_id}_{part.subcomponent_id}"
        for row in part.mapping_rows:
            mapping_rows.append({
                "locus": int(locus_idx),
                "locus_name": locus_name,
                **row,
            })
        stats_rows.append({
            "locus": int(locus_idx),
            "locus_name": locus_name,
            **part.stats_row,
        })

    mapping_df = pd.DataFrame(mapping_rows)
    stats_df = pd.DataFrame(stats_rows).sort_values("locus").reset_index(drop=True)
    mapping_df.to_csv(out_mapping_tsv, sep="\t", float_format="%12.6f", index=False)
    stats_df.to_csv(out_stats_tsv, sep="\t", float_format="%12.6f", index=False)
    _write_component_audits(audit_dir, component_results)
    logger.info(f"wrote locus mapping to {out_mapping_tsv}")
    logger.info(f"wrote locus stats to {out_stats_tsv}")
    logger.info(f"wrote denovo component audits to {audit_dir}")
    return mapping_df, stats_df
