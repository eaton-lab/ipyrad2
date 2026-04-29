#!/usr/bin/env python

"""Helpers for deterministic external sorting in assemble pipelines."""

from __future__ import annotations

from typing import Sequence


def assemble_sort_with_args(args: Sequence[str]) -> list[str]:
    """Return an external sort command with bytewise-stable collation."""
    return ["env", "LC_ALL=C", "sort", *args]
