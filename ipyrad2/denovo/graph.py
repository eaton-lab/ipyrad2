#!/usr/bin/env python

"""Build denovo locus tables from across-sample clustering graphs.

This stage turns sample-level consensus records plus across-sample hits into
final denovo loci in five phases:

1. load per-seed summary rows and strongest edge scores
2. find raw connected components in the across-sample hit graph
3. contract same-sample technical duplicates inside each component
4. split the contracted graph into final loci, or emit an oversize placeholder
5. expand each final locus into mapping, stats, and audit outputs
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import time
from typing import Iterable

from loguru import logger
import pandas as pd
from .common import DENOVO_MAPPING_FILENAME, DENOVO_STATS_FILENAME
from .graph_split import (
    connected_components,
    sort_components,
    split_component as split_constrained_component,
    subset_edges_above_threshold,
)
from ..utils.parallel import run_with_pool_iter


RECONCILE_QCOV = 0.75
OVERSIZE_COMPONENT_FACTOR = 10
FLUSH_STALL_LOG_SECONDS = 30.0
AUDIT_SUMMARY_FIELDS = [
    "component_id",
    "n_input_nodes",
    "n_contracted_nodes",
    "n_input_samples",
    "n_duplicate_samples",
    "has_joined",
    "has_merged",
    "same_sample_reconciliation_attempted",
    "used_reconciliation",
    "reconcile_mode",
    "used_oversize_rescue",
    "used_residue_cleanup",
    "status",
    "discard_reason",
    "discard_limit_nodes",
    "n_final_loci",
]


@dataclass(frozen=True, slots=True)
class ContractedComponent:
    """One post-contraction graph component ready for splitter dispatch."""

    component_id: int
    duplicated_component: bool
    same_sample_reconciliation_attempted: bool
    used_reconciliation: bool
    node_members: dict[str, tuple[str, ...]]
    node_samples: dict[str, frozenset[str]]
    node_order: dict[str, int]
    node_modes: dict[str, str]
    edges: dict[tuple[str, str], tuple[float, float]]


@dataclass(frozen=True, slots=True)
class SummaryRecord:
    """Minimal summary fields needed for graph contraction and table writing."""

    sample: str
    n_reads: int
    n_unique: int
    length: int
    cluster_length: int
    record_type: str
    cluster_id: int
    cluster_sequence: str


@dataclass(frozen=True, slots=True)
class ComponentPart:
    """One final split result returned from one component worker."""

    component_id: int
    subcomponent_id: int
    mapping_rows: tuple[dict[str, object], ...]
    stats_row: dict[str, object]


@dataclass(frozen=True, slots=True)
class ComponentResult:
    """One processed component plus audit artifacts."""

    component_id: int
    parts: tuple[ComponentPart, ...]
    audit_summary: dict[str, object]
    audit_rows: tuple[dict[str, object], ...]
    audit_fasta: tuple[tuple[str, str], ...]
    raw_oversize_fast_path: bool = False


@dataclass(frozen=True, slots=True)
class SplitExecution:
    """Splitter output plus audit flags for one contracted component."""

    subcomponents: tuple[set[str], ...]
    used_oversize_rescue: bool
    used_residue_cleanup: bool


@dataclass(frozen=True, slots=True)
class GraphTableSummary:
    """Compact counters describing the written denovo graph tables."""

    loci_written: int = 0
    consensus_records: int = 0
    duplicated_components_seen: int = 0
    same_sample_reconciliation_attempted: int = 0
    components_reconciled: int = 0
    joined_only_reconciled_loci: int = 0
    mixed_reconciled_loci: int = 0
    mixed_reconciled_groups: int = 0
    rescued_oversize_components: int = 0
    raw_oversize_placeholder_components: int = 0
    post_contraction_oversize_placeholder_components: int = 0


def get_edges_dict(outdir: Path) -> dict[tuple[str, str], tuple[float, float]]:
    """Return graph edges mapped to their strongest `(pid, qcov)` score."""
    uc_path = outdir / "global_hits.uc.tsv"
    edges: dict[tuple[str, str], tuple[float, float]] = {}
    with open(uc_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            query, target, pid, _qstrand, qcov, _qlen, _tlen = line.rstrip().split("	")
            if query == target:
                continue
            edge = tuple(sorted((query, target)))
            score = (float(pid) / 100.0, float(qcov) / 100.0)
            # Multiple hit rows can collapse onto the same unordered seed pair.
            # Keep the strongest observed support because downstream contraction
            # and splitting only need one best edge per pair.
            prev = edges.get(edge)
            if prev is None or score > prev:
                edges[edge] = score
    return edges


def _load_summary_records(
    outdir: Path,
) -> tuple[list[str], dict[str, SummaryRecord], dict[str, int], int]:
    """Load the graph-stage summary TSV into compact per-seed structures."""
    tsv = outdir / "concat.summary.tsv"
    # The summary schema is intentionally strict here. This table is the sole
    # interchange format between within-sample clustering and graph splitting.
    required = {
        "sample",
        "cluster_id",
        "seed",
        "length",
        "cluster_length",
        "n_unique",
        "n_reads",
        "record_type",
        "cluster_sequence",
        "arm_boundary",
    }
    with open(tsv, "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        missing = sorted(required.difference(fieldnames))
        if missing:
            joined = ", ".join(missing)
            raise RuntimeError(f"concat.summary.tsv is missing required columns: {joined}")

        ordered_seeds: list[str] = []
        summary_records: dict[str, SummaryRecord] = {}
        seed_order: dict[str, int] = {}
        samples: set[str] = set()
        for idx, row in enumerate(reader, start=1):
            seed = str(row["seed"])
            sample = str(row["sample"])
            ordered_seeds.append(seed)
            seed_order[seed] = idx
            samples.add(sample)
            summary_records[seed] = SummaryRecord(
                sample=sample,
                n_reads=int(row["n_reads"]),
                n_unique=int(row["n_unique"]),
                length=int(row["length"]),
                cluster_length=int(row["cluster_length"]),
                record_type=str(row["record_type"]),
                cluster_id=int(row["cluster_id"]),
                cluster_sequence=str(row["cluster_sequence"]).upper(),
            )
    return ordered_seeds, summary_records, seed_order, len(samples)


def _component_has_duplicate_samples(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
) -> bool:
    """Return True when one component contains repeated samples."""
    counts = Counter(summary_records[node].sample for node in component)
    return any(count > 1 for count in counts.values())


def _component_input_sample_count(
    component: Iterable[str],
    summary_records: dict[str, SummaryRecord],
) -> int:
    """Return the number of distinct samples represented in one component."""
    return len({summary_records[node].sample for node in component})


def _passes_same_sample_reconcile_threshold(
    score: tuple[float, float] | None,
    within_similarity: float,
) -> bool:
    """Return True when one direct same-sample edge passes contraction filters."""
    if score is None:
        return False
    identity, overlap = score
    return identity >= within_similarity and overlap >= RECONCILE_QCOV


def _same_sample_groups(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    edges: dict[tuple[str, str], tuple[float, float]],
    within_similarity: float,
) -> list[tuple[str, tuple[str, ...], str]]:
    """Return same-sample duplicate groups eligible for contraction."""
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
        # Build a same-sample subgraph containing only edges strong enough to
        # justify technical-duplicate reconciliation. Identical stripped
        # sequences are promoted to perfect support even if the global edge
        # table does not include an explicit same-sample hit row.
        reconcile_edges: dict[tuple[str, str], tuple[float, float]] = {}
        for left_idx, left in enumerate(sample_nodes):
            for right in sample_nodes[left_idx + 1:]:
                edge = tuple(sorted((left, right)))
                left_record = summary_records[left]
                right_record = summary_records[right]
                if left_record.cluster_sequence == right_record.cluster_sequence:
                    reconcile_edges[edge] = (1.0, 1.0)
                    continue
                score = edges.get(edge)
                if _passes_same_sample_reconcile_threshold(score, within_similarity):
                    reconcile_edges[edge] = score
        for subcomponent in sort_components(
            connected_components(sample_nodes, reconcile_edges),
            seed_order,
        ):
            if len(subcomponent) <= 1:
                continue
            ordered = tuple(sorted(subcomponent, key=seed_order.get))
            # Contracted groups keep a mode tag because later output-form
            # decisions need to know whether reconciliation touched only joined
            # records or mixed joined/merged records.
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
    within_similarity: float,
) -> ContractedComponent:
    """Contract same-sample technical duplicates inside one component."""
    duplicated_component = _component_has_duplicate_samples(component, summary_records)
    contracted_groups: list[tuple[str, tuple[str, ...], str]] = []
    grouped_members: set[str] = set()

    # Same-sample contraction is intentionally local to one raw connected
    # component. We never merge records from different components here.
    if duplicated_component:
        for sample, members, mode in _same_sample_groups(
            component,
            summary_records,
            seed_order,
            edges,
            within_similarity,
        ):
            contracted_groups.append((sample, members, mode))
            grouped_members.update(members)

    used_reconciliation = bool(contracted_groups)

    for core in sorted(component, key=seed_order.get):
        if core in grouped_members:
            continue
        sample = summary_records[core].sample
        contracted_groups.append((sample, (core,), "none"))

    # Stable contracted-node ordering is what later preserves deterministic
    # locus naming and output row order.
    contracted_groups.sort(key=lambda item: min(seed_order[core] for core in item[1]))
    node_members: dict[str, tuple[str, ...]] = {}
    node_samples: dict[str, frozenset[str]] = {}
    node_order: dict[str, int] = {}
    node_modes: dict[str, str] = {}

    for local_idx, (sample, members, mode) in enumerate(contracted_groups, start=1):
        node_id = f"contract_{component_id}_{local_idx}"
        ordered_members = tuple(sorted(members, key=seed_order.get))
        node_members[node_id] = ordered_members
        node_samples[node_id] = frozenset({sample})
        node_order[node_id] = min(seed_order[core] for core in ordered_members)
        node_modes[node_id] = mode

    ordered_node_ids = sorted(node_members, key=node_order.get)
    contracted_edges: dict[tuple[str, str], tuple[float, float]] = {}
    for left_idx, left in enumerate(ordered_node_ids):
        left_members = node_members[left]
        for right in ordered_node_ids[left_idx + 1:]:
            best: tuple[float, float] | None = None
            # Lift raw-seed edges up to contracted nodes by keeping the strongest
            # supporting edge between the two same-sample groups.
            for left_core in left_members:
                for right_core in node_members[right]:
                    edge = tuple(sorted((left_core, right_core)))
                    score = edges.get(edge)
                    if score is not None and (best is None or score > best):
                        best = score
            if best is not None:
                contracted_edges[(left, right)] = best

    return ContractedComponent(
        component_id=component_id,
        duplicated_component=duplicated_component,
        same_sample_reconciliation_attempted=duplicated_component,
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
    components: list[tuple[int, set[str]]],
    edges: dict[tuple[str, str], tuple[float, float]],
) -> dict[int, dict[tuple[str, str], tuple[float, float]]]:
    """Return per-component edge dictionaries in one pass over the global edges."""
    node_to_component: dict[str, int] = {}
    for component_id, component in components:
        for node in component:
            node_to_component[node] = component_id

    component_edges = {component_id: {} for component_id, _component in components}
    for edge, score in edges.items():
        component_id = node_to_component.get(edge[0])
        if component_id is None:
            continue
        component_edges[component_id][edge] = score
    return component_edges


def _std_or_nan(values: list[float | int]) -> float:
    """Return pandas-compatible sample std, including NaN for one value."""
    return float(pd.Series(values, dtype="float64").std())


def _round_locus_stat(value: float) -> float:
    """Clamp selected floating-point locus stats to 3 decimals."""
    return round(float(value), 3)


def _part_reconcile_mode(
    ordered_nodes: tuple[str, ...],
    contracted: ContractedComponent,
) -> str:
    """Return the reconcile mode label for one final locus part."""
    if not contracted.used_reconciliation:
        return "none"
    modes = [
        contracted.node_modes[node]
        for node in ordered_nodes
        if len(contracted.node_members[node]) > 1
    ]
    if any(mode == "mixed" for mode in modes):
        return "mixed"
    if any(mode == "joined_only" for mode in modes):
        return "joined_only"
    return "none"


def _part_output_form(
    *,
    mapping_rows: list[dict[str, object]],
    ordered_nodes: tuple[str, ...],
    contracted: ContractedComponent,
) -> str:
    """Return the output form for one final locus part."""
    if all(row["record_type"] == "joined" for row in mapping_rows):
        return "spaced"
    if any(contracted.node_modes[node] in {"joined_only", "mixed"} for node in ordered_nodes):
        return "spaced"
    return "stripped"


def _split_contracted_component(
    contracted: ContractedComponent,
    *,
    max_component_nodes: int,
) -> SplitExecution | None:
    """Return split loci for one contracted graph, or None for oversize fallback."""
    contracted_nodes = set(contracted.node_members)
    if len(contracted_nodes) > max_component_nodes:
        rescued = _presplit_oversize_component(contracted, max_component_nodes)
        if rescued is None:
            return None

        subcomponents: list[set[str]] = []
        used_residue_cleanup = False
        # Oversize rescue first partitions the contracted graph into pieces that
        # fit under the node cap, then runs constrained splitting on each piece.
        for rescued_nodes in rescued:
            split_details = split_constrained_component(
                set(rescued_nodes),
                _subset_component_edges(rescued_nodes, contracted.edges),
                contracted.node_samples,
                contracted.node_order,
                return_details=True,
            )
            subcomponents.extend(split_details.components)
            used_residue_cleanup = used_residue_cleanup or split_details.used_residue_cleanup
        return SplitExecution(
            subcomponents=tuple(subcomponents),
            used_oversize_rescue=True,
            used_residue_cleanup=used_residue_cleanup,
        )

    # The constrained splitter can report whether it needed its residue cleanup
    # fallback; preserve that for audit summaries.
    split_details = split_constrained_component(
        contracted_nodes,
        contracted.edges,
        contracted.node_samples,
        contracted.node_order,
        return_details=True,
    )
    return SplitExecution(
        subcomponents=split_details.components,
        used_oversize_rescue=False,
        used_residue_cleanup=split_details.used_residue_cleanup,
    )


def _iter_ordered_cores(
    ordered_nodes: tuple[str, ...],
    contracted: ContractedComponent,
    seed_order: dict[str, int],
) -> list[tuple[str, str]]:
    """Return `(contract_group, core)` rows in stable global seed order."""
    ordered_cores = [
        (node, core)
        for node in ordered_nodes
        for core in contracted.node_members[node]
    ]
    ordered_cores.sort(key=lambda item: seed_order[item[1]])
    return ordered_cores


def _build_component_part(
    *,
    component_id: int,
    subcomponent_id: int,
    nodes: set[str],
    contracted: ContractedComponent,
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
) -> tuple[ComponentPart, tuple[dict[str, object], ...]]:
    """Expand one contracted split into mapping, stats, and audit rows."""
    _validate_sample_uniqueness(nodes, contracted.node_samples)
    ordered_nodes = tuple(sorted(nodes, key=contracted.node_order.get))
    ordered_cores = _iter_ordered_cores(ordered_nodes, contracted, seed_order)
    mapping_rows: list[dict[str, object]] = []
    n_reads_values: list[int] = []
    length_values: list[int] = []
    merged_values: list[int] = []
    samples = [next(iter(contracted.node_samples[node])) for node in ordered_nodes]
    reconciled_nodes = [node for node in ordered_nodes if len(contracted.node_members[node]) > 1]

    # Mapping rows stay at raw-core granularity even after same-sample
    # contraction. Contracted groups are carried along explicitly so later
    # alignment code can collapse reconciled groups inside one final locus.
    for contract_group, core in ordered_cores:
        info = summary_records[core]
        merged = int(info.record_type == "merged")
        n_reads_values.append(info.n_reads)
        length_values.append(info.length)
        merged_values.append(merged)
        mapping_rows.append({
            "component_id": int(component_id),
            "subcomponent_id": int(subcomponent_id),
            "contract_group": contract_group,
            "sample": info.sample,
            "n_reads": info.n_reads,
            "n_unique": info.n_unique,
            "length": info.length,
            "cluster_length": info.cluster_length,
            "merged": merged,
            "record_type": info.record_type,
            "cluster_id": info.cluster_id,
            "core": core,
        })

    reconcile_mode = _part_reconcile_mode(ordered_nodes, contracted)
    output_form = _part_output_form(
        mapping_rows=mapping_rows,
        ordered_nodes=ordered_nodes,
        contracted=contracted,
    )
    # Locus-level flags are derived after mapping rows exist because output-form
    # decisions depend on both record types and reconciliation mode tags.
    n_reconciled_groups = len(reconciled_nodes)
    n_mixed_records = sum(
        len(contracted.node_members[node])
        for node in reconciled_nodes
        if contracted.node_modes[node] == "mixed"
    )
    audit_rows: list[dict[str, object]] = []
    for row in mapping_rows:
        # Audit rows intentionally mirror mapping rows with a little extra
        # reconciliation context for duplicated components.
        row["reconcile_mode"] = reconcile_mode
        row["reconciled_group"] = (
            row["contract_group"]
            if len(contracted.node_members[row["contract_group"]]) > 1
            else ""
        )
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
        "n_contracted_groups": int(n_reconciled_groups),
        "n_reconciled_groups": int(n_reconciled_groups),
        "n_mixed_records": int(n_mixed_records),
        "n_reads_sum": int(sum(n_reads_values)),
        "n_reads_mean": _round_locus_stat(sum(n_reads_values) / len(n_reads_values)),
        "n_reads_std": _round_locus_stat(_std_or_nan(n_reads_values)),
        "length_mean": _round_locus_stat(sum(length_values) / len(length_values)),
        "length_std": _round_locus_stat(_std_or_nan(length_values)),
        "merged_freq": _round_locus_stat(sum(merged_values) / len(merged_values)),
        "duplicated_component": bool(contracted.duplicated_component),
        "same_sample_reconciliation_attempted": bool(contracted.same_sample_reconciliation_attempted),
        "used_reconciliation": bool(contracted.used_reconciliation),
        "reconcile_mode": reconcile_mode,
        "used_joined_only_reconciliation": bool(reconcile_mode == "joined_only"),
        "used_mixed_reconciliation": bool(reconcile_mode == "mixed"),
        "output_form": output_form,
        "samples": ",".join(samples),
    }
    return (
        ComponentPart(
            component_id=component_id,
            subcomponent_id=subcomponent_id,
            mapping_rows=tuple(mapping_rows),
            stats_row=stats_row,
        ),
        tuple(audit_rows),
    )


def _make_component_audit_summary(
    *,
    component_id: int,
    component: tuple[str, ...],
    contracted: ContractedComponent,
    summary_records: dict[str, SummaryRecord],
    split_execution: SplitExecution,
    n_final_loci: int,
) -> dict[str, object]:
    """Return one compact audit summary row for a processed component."""
    counts = Counter(summary_records[node].sample for node in component)
    component_mode = _part_reconcile_mode(
        tuple(sorted(contracted.node_members, key=contracted.node_order.get)),
        contracted,
    )
    return {
        "component_id": int(component_id),
        "n_input_nodes": int(len(component)),
        "n_contracted_nodes": int(len(contracted.node_members)),
        "n_input_samples": int(len(counts)),
        "n_duplicate_samples": int(sum(1 for value in counts.values() if value > 1)),
        "has_joined": bool(any(summary_records[node].record_type == "joined" for node in component)),
        "has_merged": bool(any(summary_records[node].record_type == "merged" for node in component)),
        "same_sample_reconciliation_attempted": bool(contracted.same_sample_reconciliation_attempted),
        "used_reconciliation": bool(contracted.used_reconciliation),
        "reconcile_mode": component_mode,
        "used_oversize_rescue": bool(split_execution.used_oversize_rescue),
        "used_residue_cleanup": bool(split_execution.used_residue_cleanup),
        "status": "processed",
        "discard_reason": "",
        "discard_limit_nodes": 0,
        "n_final_loci": int(n_final_loci),
    }


def _process_component(
    component_id: int,
    component: tuple[str, ...],
    summary_records: dict[str, SummaryRecord],
    component_edges: dict[tuple[str, str], tuple[float, float]],
    seed_order: dict[str, int],
    within_similarity: float,
    max_component_nodes: int,
) -> tuple[int, ComponentResult]:
    """Contract, split, and expand one connected component."""
    # Each worker performs the full lifecycle for one raw connected component:
    # contract same-sample duplicates, split the contracted graph, then expand
    # final loci back into raw-core rows and compact audit artifacts.
    contracted = _build_contracted_component(
        component_id,
        set(component),
        summary_records,
        component_edges,
        seed_order,
        within_similarity=within_similarity,
    )
    split_execution = _split_contracted_component(
        contracted,
        max_component_nodes=max_component_nodes,
    )
    if split_execution is None:
        # Oversize placeholder loci intentionally bypass graph splitting while
        # still keeping one representative per sample in the pseudoreference.
        return component_id, _make_oversize_component_result(
            component_id=component_id,
            component=set(component),
            contracted=contracted,
            summary_records=summary_records,
            seed_order=seed_order,
            discard_limit_nodes=max_component_nodes,
        )

    parts: list[ComponentPart] = []
    audit_rows: list[dict[str, object]] = []
    for subcomponent_id, nodes in enumerate(split_execution.subcomponents, start=1):
        part, part_audit_rows = _build_component_part(
            component_id=component_id,
            subcomponent_id=subcomponent_id,
            nodes=nodes,
            contracted=contracted,
            summary_records=summary_records,
            seed_order=seed_order,
        )
        parts.append(part)
        audit_rows.extend(part_audit_rows)

    audit_summary = _make_component_audit_summary(
        component_id=component_id,
        component=component,
        contracted=contracted,
        summary_records=summary_records,
        split_execution=split_execution,
        n_final_loci=len(parts),
    )
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


def _process_raw_oversize_component(
    component_id: int,
    component: tuple[str, ...],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    max_component_nodes: int,
) -> tuple[int, ComponentResult]:
    """Return one placeholder locus for a raw component above the node cap."""
    return component_id, _make_raw_oversize_component_result(
        component_id=component_id,
        component=set(component),
        summary_records=summary_records,
        seed_order=seed_order,
        discard_limit_nodes=max_component_nodes,
    )


def _pick_oversize_representatives(
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
) -> tuple[str, ...]:
    """Return one best raw seed per sample from one oversized component."""
    best_by_sample: dict[str, str] = {}

    def _score(seed: str) -> tuple[int, int, int]:
        record = summary_records[seed]
        # Prefer deeper, longer representatives. Stable seed order breaks ties
        # so placeholder-locus selection remains deterministic across runs.
        return (
            int(record.n_reads),
            int(record.cluster_length),
            -int(seed_order[seed]),
        )

    for seed in sorted(component, key=seed_order.get):
        sample = summary_records[seed].sample
        current = best_by_sample.get(sample)
        if current is None or _score(seed) > _score(current):
            best_by_sample[sample] = seed
    return tuple(sorted(best_by_sample.values(), key=seed_order.get))


def _oversize_output_form(
    representatives: tuple[str, ...],
    summary_records: dict[str, SummaryRecord],
) -> str:
    """Return the output form for one oversized unsplit placeholder locus."""
    if representatives and all(summary_records[seed].record_type == "joined" for seed in representatives):
        return "spaced"
    return "stripped"


def _make_placeholder_component_result(
    *,
    component_id: int,
    component: set[str],
    representatives: tuple[str, ...],
    summary_records: dict[str, SummaryRecord],
    discard_limit_nodes: int,
    n_contracted_nodes: int,
    same_sample_reconciliation_attempted: bool,
    used_reconciliation: bool,
    raw_oversize_fast_path: bool,
) -> ComponentResult:
    """Build one placeholder locus from the best representative of each sample."""
    output_form = _oversize_output_form(representatives, summary_records)
    samples = [summary_records[seed].sample for seed in representatives]
    n_reads_values = [summary_records[seed].n_reads for seed in representatives]
    length_values = [summary_records[seed].length for seed in representatives]
    merged_values = [
        int(summary_records[seed].record_type == "merged")
        for seed in representatives
    ]
    mapping_rows: list[dict[str, object]] = []

    for local_idx, seed in enumerate(representatives, start=1):
        record = summary_records[seed]
        merged = int(record.record_type == "merged")
        mapping_rows.append({
            "component_id": int(component_id),
            "subcomponent_id": 1,
            "contract_group": f"contract_{component_id}_{local_idx}",
            "sample": record.sample,
            "n_reads": record.n_reads,
            "n_unique": record.n_unique,
            "length": record.length,
            "cluster_length": record.cluster_length,
            "merged": merged,
            "record_type": record.record_type,
            "cluster_id": record.cluster_id,
            "core": seed,
            "reconcile_mode": "none",
            "reconciled_group": "",
            "output_form": output_form,
        })

    stats_row = {
        "component_id": int(component_id),
        "subcomponent_id": 1,
        "n_samples": len(samples),
        "n_cores": len(representatives),
        "n_contracted_groups": 0,
        "n_reconciled_groups": 0,
        "n_mixed_records": 0,
        "n_reads_sum": int(sum(n_reads_values)),
        "n_reads_mean": _round_locus_stat(sum(n_reads_values) / len(n_reads_values)),
        "n_reads_std": _round_locus_stat(_std_or_nan(n_reads_values)),
        "length_mean": _round_locus_stat(sum(length_values) / len(length_values)),
        "length_std": _round_locus_stat(_std_or_nan(length_values)),
        "merged_freq": _round_locus_stat(sum(merged_values) / len(merged_values)),
        "duplicated_component": bool(_component_has_duplicate_samples(component, summary_records)),
        "same_sample_reconciliation_attempted": bool(same_sample_reconciliation_attempted),
        "used_reconciliation": bool(used_reconciliation),
        "reconcile_mode": "none",
        "used_joined_only_reconciliation": False,
        "used_mixed_reconciliation": False,
        "output_form": output_form,
        "samples": ",".join(samples),
    }

    counts = Counter(summary_records[node].sample for node in component)
    audit_summary = {
        "component_id": int(component_id),
        "n_input_nodes": int(len(component)),
        "n_contracted_nodes": int(n_contracted_nodes),
        "n_input_samples": int(len(counts)),
        "n_duplicate_samples": int(sum(1 for value in counts.values() if value > 1)),
        "has_joined": bool(any(summary_records[node].record_type == "joined" for node in component)),
        "has_merged": bool(any(summary_records[node].record_type == "merged" for node in component)),
        "same_sample_reconciliation_attempted": bool(same_sample_reconciliation_attempted),
        "used_reconciliation": bool(used_reconciliation),
        "reconcile_mode": "none",
        "used_oversize_rescue": False,
        "used_residue_cleanup": False,
        "status": "oversize_unsplit",
        "discard_reason": "",
        "discard_limit_nodes": int(discard_limit_nodes),
        "n_final_loci": 1,
    }
    return ComponentResult(
        component_id=component_id,
        parts=(
            ComponentPart(
                component_id=component_id,
                subcomponent_id=1,
                mapping_rows=tuple(mapping_rows),
                stats_row=stats_row,
            ),
        ),
        audit_summary=audit_summary,
        audit_rows=(),
        audit_fasta=(),
        raw_oversize_fast_path=raw_oversize_fast_path,
    )


def _make_raw_oversize_component_result(
    *,
    component_id: int,
    component: set[str],
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    discard_limit_nodes: int,
) -> ComponentResult:
    """Return one placeholder locus for a raw component above the node cap."""
    representatives = _pick_oversize_representatives(component, summary_records, seed_order)
    return _make_placeholder_component_result(
        component_id=component_id,
        component=component,
        representatives=representatives,
        summary_records=summary_records,
        discard_limit_nodes=discard_limit_nodes,
        n_contracted_nodes=len(representatives),
        same_sample_reconciliation_attempted=False,
        used_reconciliation=False,
        raw_oversize_fast_path=True,
    )


def _make_oversize_component_result(
    *,
    component_id: int,
    component: set[str],
    contracted: ContractedComponent,
    summary_records: dict[str, SummaryRecord],
    seed_order: dict[str, int],
    discard_limit_nodes: int,
) -> ComponentResult:
    """Return one placeholder-locus result for an oversized unsplit component."""
    # Oversize components are retained as one placeholder locus so later read
    # mapping still has a sink for those reads, even though the component is too
    # large to split safely at denovo time.
    representatives = _pick_oversize_representatives(component, summary_records, seed_order)
    return _make_placeholder_component_result(
        component_id=component_id,
        component=component,
        representatives=representatives,
        summary_records=summary_records,
        discard_limit_nodes=discard_limit_nodes,
        n_contracted_nodes=len(contracted.node_members),
        same_sample_reconciliation_attempted=bool(contracted.same_sample_reconciliation_attempted),
        used_reconciliation=bool(contracted.used_reconciliation),
        raw_oversize_fast_path=False,
    )


def _write_component_audits(
    audit_dir: Path,
) -> Path:
    """Prepare compact audit outputs and return the summary TSV path."""
    if audit_dir.exists():
        for path in audit_dir.iterdir():
            if path.is_dir():
                continue
            path.unlink()
    else:
        audit_dir.mkdir(parents=True, exist_ok=True)
    summary_path = audit_dir / "components.summary.tsv"
    with open(summary_path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(
            out,
            delimiter="	",
            fieldnames=AUDIT_SUMMARY_FIELDS,
        )
        writer.writeheader()
    return summary_path


def _write_component_audit_result(
    audit_dir: Path,
    summary_path: Path,
    result: ComponentResult,
) -> None:
    """Append one duplicated-component audit result to disk."""
    if (
        result.audit_summary["status"] == "processed"
        and result.audit_summary["n_duplicate_samples"] <= 0
    ):
        return

    with open(summary_path, "at", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="	", fieldnames=AUDIT_SUMMARY_FIELDS)
        writer.writerow(result.audit_summary)

    if result.audit_summary["status"] != "processed":
        return

    members_path = audit_dir / f"component_{result.component_id}.members.tsv"
    with open(members_path, "wt", encoding="utf-8", newline="") as out:
        if result.audit_rows:
            writer = csv.DictWriter(out, delimiter="	", fieldnames=list(result.audit_rows[0]))
            writer.writeheader()
            writer.writerows(result.audit_rows)

    fasta_path = audit_dir / f"component_{result.component_id}.fa"
    with open(fasta_path, "wt", encoding="utf-8") as out:
        for core, seq in result.audit_fasta:
            out.write(f">{core}\n{seq}\n")


def _strip_component_audit_payload(result: ComponentResult) -> ComponentResult:
    """Drop bulky audit payloads once they have been written to disk."""
    if not result.audit_rows and not result.audit_fasta:
        return result
    return ComponentResult(
        component_id=result.component_id,
        parts=result.parts,
        audit_summary=result.audit_summary,
        audit_rows=(),
        audit_fasta=(),
        raw_oversize_fast_path=result.raw_oversize_fast_path,
    )


def _iter_component_jobs(
    components: list[tuple[int, set[str]]],
    *,
    summary_records: dict[str, SummaryRecord],
    component_edges: dict[int, dict[tuple[str, str], tuple[float, float]]],
    raw_oversize_component_ids: set[int],
    seed_order: dict[str, int],
    within_similarity: float,
    max_component_nodes: int,
):
    """Yield per-component jobs lazily to avoid duplicating full graph state."""
    for component_id, component in components:
        ordered = tuple(sorted(component, key=seed_order.get))
        if component_id in raw_oversize_component_ids:
            yield component_id, (
                _process_raw_oversize_component,
                {
                    "component_id": component_id,
                    "component": ordered,
                    "summary_records": {seed: summary_records[seed] for seed in ordered},
                    "seed_order": {seed: seed_order[seed] for seed in ordered},
                    "max_component_nodes": max_component_nodes,
                },
            )
            continue
        yield component_id, (
            _process_component,
            {
                "component_id": component_id,
                "component": ordered,
                # Ship only the per-component slices each worker needs rather
                # than broadcasting the full summary/edge tables to every job.
                "summary_records": {seed: summary_records[seed] for seed in ordered},
                "component_edges": component_edges[component_id],
                "seed_order": {seed: seed_order[seed] for seed in ordered},
                "within_similarity": within_similarity,
                "max_component_nodes": max_component_nodes,
            },
        )


def _subset_component_edges(
    nodes: set[str] | tuple[str, ...],
    edges: dict[tuple[str, str], tuple[float, float]],
) -> dict[tuple[str, str], tuple[float, float]]:
    """Return one induced edge map for a node subset."""
    node_set = set(nodes)
    return {
        edge: score
        for edge, score in edges.items()
        if edge[0] in node_set and edge[1] in node_set
    }


def _presplit_oversize_component(
    contracted: ContractedComponent,
    max_component_nodes: int,
) -> list[set[str]] | None:
    """Return a valid PID presplit for one oversize contracted graph."""
    if len({next(iter(samples)) for samples in contracted.node_samples.values()}) <= 1:
        return None
    thresholds = sorted({pid for pid, _qcov in contracted.edges.values()})
    best_partition: list[set[str]] | None = None
    best_score: tuple[int, int, int, float] | None = None
    nodes = set(contracted.node_members)

    for threshold in thresholds:
        subset = subset_edges_above_threshold(nodes, contracted.edges, threshold)
        subcomponents = sort_components(
            connected_components(nodes, subset),
            contracted.node_order,
        )
        if any(len(component) > max_component_nodes for component in subcomponents):
            continue
        # Prefer partitions that preserve larger multi-node loci while avoiding
        # unnecessary fragmentation into many small components or singletons.
        score = (
            sum(len(component) for component in subcomponents if len(component) > 1),
            -len(subcomponents),
            -sum(1 for component in subcomponents if len(component) == 1),
            -float(threshold),
        )
        if best_score is None or score > best_score:
            best_score = score
            best_partition = subcomponents
    return best_partition


def make_global_tables(
    outdir: Path,
    *,
    cores: int = 1,
    log_level: str = "INFO",
    within_similarity: float = 0.95,
) -> GraphTableSummary:
    """Build locus mapping/stats tables from the across-sample denovo graph."""
    if cores < 1:
        raise ValueError("cores must be >= 1")
    out_mapping_tsv = outdir.parent / DENOVO_MAPPING_FILENAME
    out_stats_tsv = outdir.parent / DENOVO_STATS_FILENAME
    audit_dir = outdir.parent / "denovo.audit"

    ordered_seeds, summary_records, seed_order, sample_count = _load_summary_records(outdir)
    # return graph of strongest edges for every seed pair [big data object]
    edges = get_edges_dict(outdir)

    # Component ids are assigned in stable seed order and must stay stable even
    # when later processing runs in parallel.
    components = sort_components(
        connected_components(ordered_seeds, edges),
        seed_order,
    )
    max_component_nodes = OVERSIZE_COMPONENT_FACTOR * sample_count

    component_entries = list(enumerate(components, start=1))
    raw_oversize_component_ids = {
        component_id
        for component_id, component in component_entries
        if len(component) > max_component_nodes
    }
    component_edges = _bucket_component_edges(
        [
            (component_id, component)
            for component_id, component in component_entries
            if component_id not in raw_oversize_component_ids
        ],
        edges,
    )
    audit_summary_path = _write_component_audits(audit_dir)
    component_meta = {
        component_id: {
            "raw_nodes": len(component),
            "input_samples": _component_input_sample_count(component, summary_records),
        }
        for component_id, component in component_entries
    }

    # counters for ...
    total_loci = 0
    total_mapping_rows = 0
    next_component_id = 1
    pending_results: dict[int, ComponentResult] = {}
    mapping_writer: csv.DictWriter | None = None
    stats_writer: csv.DictWriter | None = None
    duplicated_components_seen = 0
    same_sample_reconciliation_attempted = 0
    components_reconciled = 0
    joined_only_reconciled_loci = 0
    mixed_reconciled_loci = 0
    mixed_reconciled_groups = 0
    raw_oversize_results = 0
    raw_oversize_nodes = 0
    contracted_oversize_results = 0
    contracted_oversize_nodes = 0
    rescued_oversize_results = 0
    blocked_flush_started_at: float | None = None
    last_blocked_flush_log_at: float | None = None

    def _flush_result(result: ComponentResult) -> int:
        nonlocal mapping_writer, stats_writer, total_loci, total_mapping_rows
        loci_written = 0
        for part in result.parts:
            total_loci += 1
            loci_written += 1
            total_mapping_rows += len(part.mapping_rows)
            locus_name = f"locus_{part.component_id}_{part.subcomponent_id}"
            for row in part.mapping_rows:
                data = {
                    "locus": int(total_loci),
                    "locus_name": locus_name,
                    **row,
                }
                if mapping_writer is None:
                    mapping_writer = csv.DictWriter(
                        mapping_out,
                        delimiter="	",
                        fieldnames=list(data),
                    )
                    mapping_writer.writeheader()
                mapping_writer.writerow(data)
            stats_data = {
                "locus": int(total_loci),
                "locus_name": locus_name,
                **part.stats_row,
            }
            if stats_writer is None:
                stats_writer = csv.DictWriter(
                    stats_out,
                    delimiter="	",
                    fieldnames=list(stats_data),
                )
                stats_writer.writeheader()
            stats_writer.writerow(stats_data)
        return loci_written

    def _flush_pending_results() -> None:
        nonlocal next_component_id
        # Parallel workers can finish out of order, but final locus numbering
        # and output row order must follow original component order. Buffer
        # completed results until the next expected component is available.
        while next_component_id in pending_results:
            _flush_result(pending_results.pop(next_component_id))
            next_component_id += 1

    def _maybe_log_blocked_flush() -> None:
        nonlocal blocked_flush_started_at, last_blocked_flush_log_at
        if not pending_results or next_component_id in pending_results:
            blocked_flush_started_at = None
            last_blocked_flush_log_at = None
            return
        now = time.monotonic()
        if blocked_flush_started_at is None:
            blocked_flush_started_at = now
            return
        if now - blocked_flush_started_at < FLUSH_STALL_LOG_SECONDS:
            return
        if (
            last_blocked_flush_log_at is not None
            and now - last_blocked_flush_log_at < FLUSH_STALL_LOG_SECONDS
        ):
            return
        meta = component_meta[next_component_id]
        logger.warning(
            f"waiting to flush split outputs for component {next_component_id} "
            f"(raw_nodes={meta['raw_nodes']}, input_samples={meta['input_samples']}); "
            f"{len(pending_results)} later component results already buffered"
        )
        last_blocked_flush_log_at = now

    def _record_completed_result(result: ComponentResult) -> None:
        nonlocal duplicated_components_seen, same_sample_reconciliation_attempted
        nonlocal components_reconciled, joined_only_reconciled_loci
        nonlocal mixed_reconciled_loci, mixed_reconciled_groups
        nonlocal raw_oversize_results, raw_oversize_nodes
        nonlocal contracted_oversize_results, contracted_oversize_nodes
        nonlocal rescued_oversize_results
        _write_component_audit_result(audit_dir, audit_summary_path, result)
        duplicated_components_seen += int(int(result.audit_summary["n_duplicate_samples"]) > 0)
        same_sample_reconciliation_attempted += int(
            bool(result.audit_summary["same_sample_reconciliation_attempted"])
        )
        components_reconciled += int(bool(result.audit_summary["used_reconciliation"]))
        if result.audit_summary["status"] == "oversize_unsplit":
            if result.raw_oversize_fast_path:
                raw_oversize_results += 1
                raw_oversize_nodes += int(result.audit_summary["n_input_nodes"])
            else:
                contracted_oversize_results += 1
                contracted_oversize_nodes += int(result.audit_summary["n_input_nodes"])
        if result.audit_summary["used_oversize_rescue"]:
            rescued_oversize_results += 1
        for part in result.parts:
            reconcile_mode = str(part.stats_row.get("reconcile_mode", "none"))
            joined_only_reconciled_loci += int(reconcile_mode == "joined_only")
            mixed_reconciled_loci += int(reconcile_mode == "mixed")
            mixed_reconciled_groups += int(part.stats_row.get("n_reconciled_groups", 0))
        pending_results[result.component_id] = _strip_component_audit_payload(result)
        _flush_pending_results()
        _maybe_log_blocked_flush()

    with open(out_mapping_tsv, "wt", encoding="utf-8", newline="") as mapping_out, open(
        out_stats_tsv,
        "wt",
        encoding="utf-8",
        newline="",
    ) as stats_out:
        if component_entries:
            jobs_iter = _iter_component_jobs(
                component_entries,
                summary_records=summary_records,
                component_edges=component_edges,
                raw_oversize_component_ids=raw_oversize_component_ids,
                seed_order=seed_order,
                within_similarity=within_similarity,
                max_component_nodes=max_component_nodes,
            )
            if cores == 1 or len(component_entries) == 1:
                # Keep the serial path explicit so the output-ordering logic is
                # identical to the parallel path, just without worker dispatch.
                for _component_id, (func, kwargs) in jobs_iter:
                    _key, result = func(**kwargs)
                    _record_completed_result(result)
            else:
                for _key, (_component_id, result) in run_with_pool_iter(
                    jobs_iter,
                    log_level=log_level,
                    max_workers=min(cores, len(component_entries)),
                    msg="Splitting global clusters",
                    njobs=len(component_entries),
                ):
                    _record_completed_result(result)

        if mapping_writer is None:
            pd.DataFrame().to_csv(out_mapping_tsv, sep="	", index=False)
        if stats_writer is None:
            pd.DataFrame().to_csv(out_stats_tsv, sep="	", index=False)

    if raw_oversize_results:
        logger.warning(
            f"retaining {raw_oversize_results} raw oversize clusters as unsplit placeholder loci "
            f"(limit={max_component_nodes} raw nodes; total_oversize_nodes={raw_oversize_nodes})"
        )
    if contracted_oversize_results:
        logger.warning(
            f"retaining {contracted_oversize_results} post-contraction oversize clusters as unsplit placeholder loci "
            f"(limit={max_component_nodes} contracted nodes; total_oversize_nodes={contracted_oversize_nodes})"
        )

    logger.info(
        "built {} loci from {} graph components (rescued oversize: {}, raw placeholders: {}, post-contraction placeholders: {})",
        total_loci,
        len(component_entries),
        rescued_oversize_results,
        raw_oversize_results,
        contracted_oversize_results,
    )
    logger.debug("wrote locus mapping to {}", out_mapping_tsv)
    logger.debug("wrote locus stats to {}", out_stats_tsv)
    logger.debug("wrote denovo component audits to {}", audit_dir)
    return GraphTableSummary(
        loci_written=total_loci,
        consensus_records=total_mapping_rows,
        duplicated_components_seen=duplicated_components_seen,
        same_sample_reconciliation_attempted=same_sample_reconciliation_attempted,
        components_reconciled=components_reconciled,
        joined_only_reconciled_loci=joined_only_reconciled_loci,
        mixed_reconciled_loci=mixed_reconciled_loci,
        mixed_reconciled_groups=mixed_reconciled_groups,
        rescued_oversize_components=rescued_oversize_results,
        raw_oversize_placeholder_components=raw_oversize_results,
        post_contraction_oversize_placeholder_components=contracted_oversize_results,
    )
