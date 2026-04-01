#!/usr/bin/env python

"""Helpers for compact command logging."""

from __future__ import annotations


MAX_LOGGED_MATCHED_PATH_CHARS = 240

_PATH_LIST_FLAGS_BY_SUBCOMMAND = {
    "demux": {"-d", "--fastqs"},
    "trim": {"-d", "--fastqs"},
    "denovo": {"-d", "--fastqs"},
    "map": {"-d", "--fastqs"},
    "assemble": {"-d", "--rad-bams", "-w", "--wgs-bams"},
}


def _truncate_path_tokens(paths: list[str], max_chars: int) -> list[str]:
    """Truncate a list of matched paths at whole-token boundaries."""
    if not paths:
        return []

    included: list[str] = []
    for path in paths:
        candidate = " ".join(included + [path])
        if len(candidate) <= max_chars or not included:
            included.append(path)
            continue
        break

    omitted = len(paths) - len(included)
    if not omitted:
        return paths

    marker = f"...[truncated; {len(paths)} total matched paths]"
    return included + [marker]


def format_logged_command(
    argv: list[str] | tuple[str, ...],
    max_path_chars: int = MAX_LOGGED_MATCHED_PATH_CHARS,
) -> str:
    """Format a CLI command string, truncating only long matched path lists."""
    tokens = [str(token) for token in argv]
    if not tokens:
        return "ipyrad2"

    path_flags = _PATH_LIST_FLAGS_BY_SUBCOMMAND.get(tokens[0], set())
    formatted: list[str] = []
    idx = 0

    while idx < len(tokens):
        token = tokens[idx]
        formatted.append(token)
        idx += 1

        if token not in path_flags:
            continue

        paths: list[str] = []
        while idx < len(tokens) and not tokens[idx].startswith("-"):
            paths.append(tokens[idx])
            idx += 1

        formatted.extend(_truncate_path_tokens(paths, max_path_chars))

    return f"ipyrad2 {' '.join(formatted)}"
