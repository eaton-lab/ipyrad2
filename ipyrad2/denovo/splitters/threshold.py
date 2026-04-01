#!/usr/bin/env python

"""Legacy threshold-sweep denovo graph splitter."""

from __future__ import annotations

from .common import connected_components, has_duplicate_samples, sort_components, subset_edges_above_threshold


def split_component(
    nodes: set[str],
    edges: dict[tuple[str, str], tuple[float, float]],
    node_samples: dict[str, frozenset[str]],
    node_order: dict[str, int],
) -> list[set[str]]:
    """Split one component by ascending PID threshold."""
    if len(nodes) <= 1:
        return [set(nodes)] if nodes else []

    remaining = set(nodes)
    thresholds = sorted({pid for pid, _qcov in edges.values()})
    out: list[set[str]] = []

    for threshold in thresholds:
        subset = subset_edges_above_threshold(remaining, edges, threshold)
        subcomponents = connected_components(remaining, subset)
        to_remove: set[str] = set()
        for component in sort_components(subcomponents, node_order):
            if has_duplicate_samples(component, node_samples):
                continue
            out.append(component)
            to_remove.update(component)
        remaining.difference_update(to_remove)
        if not remaining:
            break

    if remaining:
        for node in sorted(remaining, key=node_order.get):
            out.append({node})

    return sort_components(out, node_order)
