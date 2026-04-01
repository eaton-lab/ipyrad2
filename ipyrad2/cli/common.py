#!/usr/bin/env python

"""Shared CLI parser helpers."""

from __future__ import annotations

import argparse
from argparse import ArgumentTypeError
from functools import partial


RAW_HELP_FORMATTER = partial(
    argparse.RawDescriptionHelpFormatter,
    width=140,
    max_help_position=60,
)


def intlike(s: str) -> int:
    """Allow int-like numeric input and coerce it to an int."""
    try:
        return int(round(float(s)))                # 3.5 -> 4; 1e5 -> 100_000
    except ValueError as exc:
        raise ArgumentTypeError(f"{s!r} is not a number") from exc
