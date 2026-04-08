#!/usr/bin/env python

"""Constrained helpers for denovo graph splitting."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


def connected_components(
    nodes: Iterable[str],
    edges: dict[tuple[str, str], tuple[float, float]],
) -> list[set[str]]:
    """Return connected components for one undirected graph."""
    nodes = list(nodes)
    parent = {node: node for node in nodes}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, right in edges:
        union(left, right)

    out: dict[str, set[str]] = {}
    for node in nodes:
        out.setdefault(find(node), set()).add(node)
    return list(out.values())


def has_duplicate_samples(
    component: set[str],
    node_samples: dict[str, frozenset[str]],
) -> bool:
    """Return True when a component contains repeated samples."""
    counts = Counter()
    for node in component:
        for sample in node_samples[node]:
            counts[sample] += 1
    return any(value > 1 for value in counts.values())


def subset_edges_above_threshold(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    threshold: float,
) -> dict[tuple[str, str], tuple[float, float]]:
    """Return only edges above one PID threshold for one node subset."""
    out: dict[tuple[str, str], tuple[float, float]] = {}
    for left, right in combinations(sorted(nodes), 2):
        edge = tuple(sorted((left, right)))
        value = edges.get(edge)
        if value is None:
            continue
        pid, qcov = value
        if pid >= threshold:
            out[edge] = (pid, qcov)
    return out


def sort_components(
    components: Iterable[set[str]],
    node_order: dict[str, int],
) -> list[set[str]]:
    """Return components sorted by the first node seen in canonical order."""
    return sorted(
        (set(component) for component in components),
        key=lambda component: min(node_order[node] for node in component),
    )


@dataclass(frozen=True)
class SplitDetails:
    """Split output plus constrained-specific audit flags."""

    components: tuple[set[str], ...]
    used_residue_cleanup: bool


def _duplicate_burden(
    components: list[set[str]],
    node_samples: dict[str, frozenset[str]],
) -> int:
    """Return summed duplicate excess across one set of components."""
    burden = 0
    for component in components:
        counts = Counter()
        for node in component:
            for sample in node_samples[node]:
                counts[sample] += 1
        burden += sum(count - 1 for count in counts.values() if count > 1)
    return burden


def _partition_makes_progress(
    nodes: set[str],
    clean: list[set[str]],
    residue: list[set[str]],
) -> bool:
    """Return True when a candidate partition shrinks the duplicated residue."""
    return bool(clean) or len(residue) != 1 or residue[0] != nodes


def _best_residue_partition(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    node_samples: dict[str, frozenset[str]],
    node_order: dict[str, int],
) -> tuple[list[set[str]], list[set[str]]] | None:
    """Return the best clean/residue threshold partition for one residue graph."""
    thresholds = sorted({pid for pid, _qcov in edges.values()})
    baseline_score = (_duplicate_burden([nodes], node_samples), 0, 0, 0.0)
    best_partition: tuple[list[set[str]], list[set[str]]] | None = None
    best_score: tuple[int, int, int, float] | None = None

    for threshold in thresholds:
        subset = subset_edges_above_threshold(nodes, edges, threshold)
        subcomponents = sort_components(connected_components(nodes, subset), node_order)
        clean = [component for component in subcomponents if not has_duplicate_samples(component, node_samples)]
        residue = [component for component in subcomponents if has_duplicate_samples(component, node_samples)]
        if not _partition_makes_progress(nodes, clean, residue):
            continue
        score = (
            _duplicate_burden(residue, node_samples),
            -sum(len(component) for component in clean),
            -len(clean),
            -float(threshold),
        )
        if best_score is None or score < best_score:
            best_score = score
            best_partition = (clean, residue)

    if best_partition is None or best_score is None or best_score >= baseline_score:
        return None
    return best_partition


def _split_residue(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    node_samples: dict[str, frozenset[str]],
    node_order: dict[str, int],
) -> tuple[list[set[str]], bool]:
    """Recursively clean one residue graph before falling back to singletons."""
    if len(nodes) <= 1:
        return ([set(nodes)] if nodes else []), False
    if not has_duplicate_samples(nodes, node_samples):
        return [set(nodes)], False

    residue_edges = {
        edge: score
        for edge, score in edges.items()
        if edge[0] in nodes and edge[1] in nodes
    }
    partition = _best_residue_partition(nodes, residue_edges, node_samples, node_order)
    if partition is None:
        singles = [{node} for node in sorted(nodes, key=node_order.get)]
        return singles, False

    clean, residue = partition
    out = list(clean)
    used_cleanup = True
    for component in residue:
        split_components, used_child_cleanup = _split_residue(
            component,
            residue_edges,
            node_samples,
            node_order,
        )
        out.extend(split_components)
        used_cleanup = used_cleanup or used_child_cleanup
    return sort_components(out, node_order), used_cleanup


def split_component(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    node_samples: dict[str, frozenset[str]],
    node_order: dict[str, int],
    *,
    return_details: bool = False,
) -> list[set[str]] | SplitDetails:
    """Build a constrained forest, then recover residue components."""
    parent = {node: node for node in nodes}
    cluster_samples = {node: set(node_samples[node]) for node in nodes}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if cluster_samples[left_root] & cluster_samples[right_root]:
            return
        if node_order[left_root] > node_order[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        cluster_samples[left_root].update(cluster_samples[right_root])
        cluster_samples.pop(right_root, None)

    ordered_edges = sorted(
        edges.items(),
        key=lambda item: (
            -item[1][0],
            -item[1][1],
            min(node_order[item[0][0]], node_order[item[0][1]]),
            max(node_order[item[0][0]], node_order[item[0][1]]),
        ),
    )
    for (left, right), _score in ordered_edges:
        union(left, right)

    forest_components = sort_components(
        ({node for node in nodes if find(node) == root} for root in {find(node) for node in nodes}),
        node_order,
    )

    out: list[set[str]] = []
    residue_nodes: set[str] = set()
    for component in forest_components:
        if len(component) > 1 and not has_duplicate_samples(component, node_samples):
            out.append(component)
        else:
            residue_nodes.update(component)

    used_residue_cleanup = False
    if residue_nodes:
        residue_edges = {
            edge: score
            for edge, score in edges.items()
            if edge[0] in residue_nodes and edge[1] in residue_nodes
        }
        residue_components = sort_components(connected_components(residue_nodes, residue_edges), node_order)
        for component in residue_components:
            if not has_duplicate_samples(component, node_samples):
                out.append(component)
                continue
            split_components, used_cleanup = _split_residue(
                component,
                residue_edges,
                node_samples,
                node_order,
            )
            out.extend(split_components)
            used_residue_cleanup = used_residue_cleanup or used_cleanup

    components = tuple(sort_components(out, node_order))
    if return_details:
        return SplitDetails(
            components=components,
            used_residue_cleanup=used_residue_cleanup,
        )
    return list(components)
