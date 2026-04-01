#!/usr/bin/env python

"""Shared HDF5 writer helpers for assemble output databases."""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

import numpy as np


MIB = 1024 * 1024
GIB = 1024 * MIB
DEFAULT_TOTAL_RAM_BYTES = 8 * GIB
_FAI_COLUMNS = {
    "scaffold": 0,
    "length": 1,
    "sumsize": 2,
    "a": 3,
    "b": 4,
    "line_bases": 3,
    "line_width": 4,
}


def _clamp_int(value: int, minimum: int, maximum: int) -> int:
    """Return `value` bounded to the inclusive [minimum, maximum] range."""
    return max(minimum, min(int(value), maximum))


def detect_total_ram_bytes(default: int = DEFAULT_TOTAL_RAM_BYTES) -> int:
    """Best-effort total system RAM detection using only the standard library."""
    try:
        pages = int(os.sysconf("SC_PHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        if pages > 0 and page_size > 0:
            return pages * page_size
    except (AttributeError, OSError, ValueError):
        pass
    return int(default)


def choose_hdf5_cache_settings(total_ram_bytes: int | None = None) -> dict[str, int | str]:
    """Choose bounded HDF5 raw chunk-cache settings from host RAM size."""
    if total_ram_bytes is None:
        total_ram_bytes = detect_total_ram_bytes()

    rdcc_nbytes = _clamp_int(total_ram_bytes // 128, 128 * MIB, 1 * GIB)
    if rdcc_nbytes <= 256 * MIB:
        rdcc_nslots = 524_287
    elif rdcc_nbytes <= 512 * MIB:
        rdcc_nslots = 1_000_003
    else:
        rdcc_nslots = 2_000_003
    return {
        "libver": "latest",
        "rdcc_nbytes": rdcc_nbytes,
        "rdcc_nslots": rdcc_nslots,
    }


def choose_unsigned_int_dtype(max_value: int) -> np.dtype:
    """Return the narrowest safe unsigned dtype for non-negative integer data."""
    if max_value < 0:
        raise ValueError("max_value must be >= 0")
    if max_value <= np.iinfo(np.uint32).max:
        return np.dtype(np.uint32)
    return np.dtype(np.uint64)


def format_bytes(nbytes: int) -> str:
    """Format byte counts for concise debug logging."""
    return f"{nbytes / MIB:.1f} MiB"


@lru_cache(maxsize=None)
def _read_fai_rows(fai_path: str) -> tuple[tuple[str | int, ...], ...]:
    """Read and cache the `.fai` rows for one reference."""
    rows: list[tuple[str | int, ...]] = []
    with open(fai_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed FASTA index row in {fai_path}: {line}")
            rows.append((
                fields[0],
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
                int(fields[4]),
            ))
    return tuple(rows)


def get_fai_values(reference: Path, key: str) -> np.ndarray:
    """Return one `.fai` column as a numpy array for the indexed reference."""
    try:
        column_idx = _FAI_COLUMNS[key]
    except KeyError as exc:
        keys = ", ".join(sorted(_FAI_COLUMNS))
        raise KeyError(f"Unsupported FASTA index column {key!r}; expected one of: {keys}") from exc

    fai = reference.with_suffix(reference.suffix + ".fai")
    rows = _read_fai_rows(str(fai))
    if column_idx == 0:
        return np.array([row[column_idx] for row in rows], dtype=object)
    return np.array([row[column_idx] for row in rows], dtype=np.int64)
