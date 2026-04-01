#!/usr/bin/env python

"""Sample-constrained maximum-spanning-forest denovo splitter."""

from __future__ import annotations

from .common import sort_components


def split_component(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    node_samples: dict[str, frozenset[str]],
    node_order: dict[str, int],
) -> list[set[str]]:
    """Build a sample-constrained forest by greedily accepting strong edges."""
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

    out: dict[str, set[str]] = {}
    for node in nodes:
        out.setdefault(find(node), set()).add(node)
    return sort_components(out.values(), node_order)
