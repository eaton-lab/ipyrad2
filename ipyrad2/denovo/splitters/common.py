#!/usr/bin/env python

"""Common helpers for denovo graph splitters."""

from __future__ import annotations

from collections import Counter
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
