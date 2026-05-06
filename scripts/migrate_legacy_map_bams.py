#!/usr/bin/env python

"""Migrate legacy ipyrad2 map BAMs and sidecars to current conventions."""

from __future__ import annotations

import argparse
from pathlib import Path

from ipyrad2.mapper.legacy_map_bams import migrate_legacy_map_outputs
from ipyrad2.utils.logger import set_log_level


def _parse_args() -> argparse.Namespace:
    """Return parsed CLI args for the legacy BAM migration script."""
    parser = argparse.ArgumentParser(
        description=(
            "Copy legacy ipyrad2 map BAMs and stats into a new directory while "
            "rewriting BAM RG sample names from {name}.trimmed to {name}."
        )
    )
    parser.add_argument(
        "--indir",
        type=Path,
        required=True,
        help="Directory containing legacy ipyrad2 map BAMs and optional map stats sidecars.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory to write migrated BAMs, fresh CSI indexes, and rewritten stats files.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads to pass to samtools for BAM rewrite and CSI indexing. [default=%(default)s]",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the migration plan without writing any BAMs or stats files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite conflicting targets inside --outdir if they already exist.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Stderr log level for progress messages. [default=%(default)s]",
    )
    return parser.parse_args()


def main() -> int:
    """Run the legacy map BAM migration script."""
    args = _parse_args()
    set_log_level(args.log_level)
    migrate_legacy_map_outputs(
        indir=args.indir,
        outdir=args.outdir,
        threads=args.threads,
        dry_run=args.dry_run,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
