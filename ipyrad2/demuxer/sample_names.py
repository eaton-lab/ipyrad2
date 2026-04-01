#!/usr/bin/env python

"""Helpers for demux sample names and technical replicate labels."""

from __future__ import annotations

from typing import Tuple


TECHNICAL_REPLICATE_TOKEN = "-technical-replicate-"


def split_technical_replicate_name(sample_name: str) -> Tuple[str, int | None]:
    """Return the base sample name and replicate index when the suffix is well-formed."""
    base_name, token, suffix = sample_name.rpartition(TECHNICAL_REPLICATE_TOKEN)
    if not token or not suffix.isdigit():
        return sample_name, None
    return base_name, int(suffix)


def technical_replicate_base_name(sample_name: str) -> str:
    """Return the replicate base name or the original name when not a replicate."""
    base_name, _replicate_idx = split_technical_replicate_name(sample_name)
    return base_name


def is_technical_replicate_name(sample_name: str) -> bool:
    """Return True when a sample name ends with a technical-replicate suffix."""
    _base_name, replicate_idx = split_technical_replicate_name(sample_name)
    return replicate_idx is not None


def final_output_sample_name(sample_name: str, merge_technical_replicates: bool) -> str:
    """Return the final output sample name after optional replicate merging."""
    if merge_technical_replicates:
        return technical_replicate_base_name(sample_name)
    return sample_name
