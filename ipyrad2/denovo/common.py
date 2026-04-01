#!/usr/bin/env python

"""Shared denovo sequence helpers."""

from __future__ import annotations


CLUSTER_JOINED_SPACER_LEN = 24
OUTPUT_JOINED_SPACER_LEN = 50


def infer_record_type(seed: str) -> str:
    """Infer the denovo record type from one seed/core name."""
    kind = seed.rsplit(";", 1)[-1][:1].upper()
    if kind == "M":
        return "merged"
    if kind == "J":
        return "joined"
    return "single"


def strip_joined_spacer(seq: str, spacer_len: int = CLUSTER_JOINED_SPACER_LEN) -> str:
    """Remove exactly one N-spacer from a joined record when present."""
    s = seq.upper()
    token = "N" * spacer_len
    return s.replace(token, "", 1) if token in s else s


def split_joined_sequence(
    seq: str,
    spacer_len: int = CLUSTER_JOINED_SPACER_LEN,
) -> tuple[str, str, str, bool]:
    """Return stripped sequence plus left/right arms for a joined record."""
    s = seq.upper()
    token = "N" * spacer_len
    if token not in s:
        return s, s, "", False
    left, right = s.split(token, 1)
    return left + right, left, right, True
