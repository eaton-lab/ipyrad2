#!/usr/bin/env python

"""Helpers for compact command logging."""

from __future__ import annotations

import shlex


MAX_LOGGED_MATCHED_PATH_CHARS = 240

_PATH_LIST_FLAGS_BY_SUBCOMMAND = {
    "demux": {"-d", "--fastqs"},
    "trim": {"-d", "--fastqs"},
    "denovo": {"-d", "--fastqs"},
    "map": {"-d", "--fastqs"},
    "assemble": {"-d", "--rad-bams", "-w", "--wgs-bams"},
}


def format_stats_command(argv: list[str], max_chars: int = 120) -> str:
    """Return one bounded display line; the original argv belongs in JSON.

    Placeholders and truncation are deliberately labelled as a summary, not
    presented as a runnable command. Control characters are escaped for display.
    """
    if max_chars < 32:
        raise ValueError("A stats command line needs at least 32 characters.")

    def display(token):
        return ''.join(c if c.isprintable() else c.encode('unicode_escape').decode('ascii')
                       for c in shlex.quote(str(token)))

    tokens = [display(t) for t in argv]
    escaped = any(not c.isprintable() for token in argv for c in str(token))
    full = ("Command (summary): " if escaped else "Command: ") + " ".join(tokens)
    if len(full) <= max_chars:
        return full
    subcommand = str(argv[1]) if len(argv) > 1 else ""
    path_flags = _PATH_LIST_FLAGS_BY_SUBCOMMAND.get(subcommand, set())
    unit = "BAMs" if subcommand == "assemble" else "FASTQs"
    groups = []
    idx = 0
    while idx < len(argv):
        token = str(argv[idx])
        idx += 1
        group = [display(token)]
        if token in path_flags:
            start = idx
            while idx < len(argv) and not str(argv[idx]).startswith("-"):
                idx += 1
            if idx > start:
                group.append(f"<{idx - start} {unit}>")
        elif token.startswith("-") and idx < len(argv) and not str(argv[idx]).startswith("-"):
            group.append(display(argv[idx]))
            idx += 1
        groups.append(" ".join(group))
    line = "Command (summary):"
    for group in groups:
        if len(line) + 1 + len(group) + 4 > max_chars:
            return line + " ..."
        line += " " + group
    return line


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
