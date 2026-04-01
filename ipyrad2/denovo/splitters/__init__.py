#!/usr/bin/env python

"""Denovo graph splitter registry."""

from __future__ import annotations

from ...utils.exceptions import IPyradError
from .constrained import split_component as split_component_constrained
from .threshold import split_component as split_component_threshold


SPLITTERS = {
    "threshold": split_component_threshold,
    "constrained": split_component_constrained,
}


def get_splitter(name: str):
    """Return one registered graph splitter."""
    try:
        return SPLITTERS[name]
    except KeyError as exc:
        choices = ", ".join(sorted(SPLITTERS))
        raise IPyradError(f"Unsupported denovo graph splitter '{name}'. Choose from: {choices}") from exc
