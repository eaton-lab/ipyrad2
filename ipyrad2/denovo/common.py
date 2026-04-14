#!/usr/bin/env python

"""Shared denovo sequence helpers."""

from __future__ import annotations


CLUSTER_JOINED_SPACER_LEN = 24
OUTPUT_JOINED_SPACER_LEN = 50
DENOVO_MAPPING_FILENAME = "denovo.loci.mapping.tsv"
DENOVO_STATS_FILENAME = "denovo.loci.stats.tsv"
DENOVO_SAMPLE_GRAPH_SUMMARY_FILENAME = "denovo.sample_graph_summary.tsv"


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


def get_arm_boundary(
    seq: str,
    spacer_len: int = CLUSTER_JOINED_SPACER_LEN,
) -> tuple[str, int]:
    """Return stripped cluster sequence plus left-arm boundary in stripped coordinates."""
    cluster_sequence, left_arm, _right_arm, _joined = split_joined_sequence(
        seq,
        spacer_len=spacer_len,
    )
    return cluster_sequence, len(left_arm)


def split_cluster_sequence_at_boundary(
    cluster_sequence: str,
    arm_boundary: int,
) -> tuple[str, str]:
    """Split one stripped cluster sequence into left/right arms."""
    boundary = max(0, min(int(arm_boundary), len(cluster_sequence)))
    return cluster_sequence[:boundary], cluster_sequence[boundary:]


def insert_joined_spacer(
    cluster_sequence: str,
    arm_boundary: int,
    spacer_len: int = CLUSTER_JOINED_SPACER_LEN,
) -> str:
    """Reinsert the joined-read spacer at one stripped arm boundary."""
    left_arm, right_arm = split_cluster_sequence_at_boundary(cluster_sequence, arm_boundary)
    if not right_arm:
        return cluster_sequence
    return left_arm + ("N" * spacer_len) + right_arm
