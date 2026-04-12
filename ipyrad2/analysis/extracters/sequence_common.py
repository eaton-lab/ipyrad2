#!/usr/bin/env python

"""Shared helpers for sequence-HDF5 sample and population selection."""

from __future__ import annotations

from pathlib import Path
import itertools

import h5py
import numpy as np
import pandas as pd
from loguru import logger

from ...utils.exceptions import IPyradError
from ...utils.pops import expand_imap_patterns, parse_imap, parse_minmap, parse_pops_file


REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"
SEQUENCE_CHUNK_SITES = 5000


def _decode_h5_names(values) -> list[str]:
    """Decode HDF5 string arrays into a list of Python strings."""
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def plan_sequence_chunk_spans(
    spans: list[tuple[int, int]] | tuple[tuple[int, int], ...],
    *,
    target_sites: int = SEQUENCE_CHUNK_SITES,
) -> list[tuple[tuple[int, int], ...]]:
    """Pack ordered phy spans into moderate-size chunks for HDF5 reads."""
    chunks: list[tuple[tuple[int, int], ...]] = []
    current: list[tuple[int, int]] = []
    current_sites = 0

    for start, end in spans:
        width = end - start
        if width <= 0:
            continue
        if current and current_sites + width > target_sites:
            chunks.append(tuple(current))
            current = []
            current_sites = 0
        if current and current[-1][1] == start:
            current[-1] = (current[-1][0], end)
        else:
            current.append((start, end))
        current_sites += width

    if current:
        chunks.append(tuple(current))
    return chunks


def load_sequence_chunk_from_phy(
    phy: h5py.Dataset,
    sidxs: list[int],
    spans: tuple[tuple[int, int], ...],
) -> np.ndarray:
    """Load one `(samples, sites)` sequence chunk from one or more phy spans."""
    total_sites = sum(end - start for start, end in spans)
    block = np.empty((len(sidxs), total_sites), dtype=np.uint8)
    offset = 0
    for start, end in spans:
        width = end - start
        block[:, offset : offset + width] = phy[sidxs, start:end]
        offset += width
    block[block == 45] = 78
    return block


def normalize_sequence_population_inputs(imap, minmap):
    """Normalize imap/minmap inputs from dicts or files."""
    if isinstance(imap, str):
        imap = Path(imap)
    if isinstance(minmap, str):
        minmap = Path(minmap)

    if imap is not None and not isinstance(imap, (Path, dict)):
        raise IPyradError("imap must be one of Path, str, or Dict")
    if minmap is not None and not isinstance(minmap, (Path, dict)):
        raise IPyradError("minmap must be one of Path, str, or Dict")

    if isinstance(minmap, Path):
        minmap = parse_minmap(minmap)

    if isinstance(imap, Path):
        parsed_minmap = None
        try:
            imap, parsed_minmap = parse_pops_file(imap)
        except IPyradError:
            logger.info(
                "imap file doesn't include minmap info, parsing standard imap file format."
            )
            imap = parse_imap(imap)
        if minmap is None:
            minmap = parsed_minmap

    return imap, minmap


def load_sequence_scaffold_table(data: Path | str) -> pd.DataFrame:
    """Return scaffold names and lengths from a sequence-capable HDF5."""
    with h5py.File(data, "r") as io5:
        if io5.attrs["version"] < 2.0:
            raise IPyradError("hdf5 database version must be >= 2.0")
        scaff_names = _decode_h5_names(io5.attrs["scaffold_names"])
        scaff_lengths = [int(value) for value in io5.attrs["scaffold_lengths"]]
    return pd.DataFrame(
        columns=["scaffold_name", "scaffold_length"],
        data={"scaffold_name": scaff_names, "scaffold_length": scaff_lengths},
    )


def load_sequence_sample_names(data: Path | str) -> list[str]:
    """Return the ordered sample names stored in a sequence HDF5."""
    with h5py.File(data, "r") as io5:
        return _decode_h5_names(io5.attrs["names"])


def resolve_sequence_sample_subset(
    data: Path | str,
    *,
    exclude=None,
    include_reference: bool = False,
    imap=None,
) -> tuple[list[str], list[int], set[str], dict[str, list[str]]]:
    """Return selected sample names, their HDF5 row indexes, and final excludes."""
    exclude_set = set(exclude if exclude else [])
    all_names = load_sequence_sample_names(data)
    dbnames = set(all_names)

    if imap:
        imap, _unmatched = expand_imap_patterns(
            imap,
            all_names,
            mapping_name="IMAP",
            available_name="the HDF5 database",
        )
        imapset = set(itertools.chain(*imap.values()))
        badnames = imapset.difference(dbnames)
        if badnames:
            badlist = ", ".join(sorted(badnames))
            raise IPyradError(
                f"IMAP samples are not in the HDF5 database: {badlist}"
            )
        if (
            include_reference
            and REFERENCE_SAMPLE_NAME not in exclude_set
            and REFERENCE_SAMPLE_NAME not in imapset
        ):
            raise IPyradError(
                "assembly_reference_sequence was requested with -R, "
                "but it must also be assigned to an IMAP group."
            )
    else:
        imap = {}
        imapset = set()

    if (
        REFERENCE_SAMPLE_NAME in dbnames
        and REFERENCE_SAMPLE_NAME not in exclude_set
        and not include_reference
        and REFERENCE_SAMPLE_NAME not in imapset
    ):
        exclude_set.add(REFERENCE_SAMPLE_NAME)

    if imap:
        exclude_set.update(set(all_names).difference(imapset))
        logger.debug(
            "dropping samples that are either not in the imap dict, or are in the exclude list: {}",
            exclude_set,
        )

    sidxs = [idx for idx, name in enumerate(all_names) if name not in exclude_set]
    snames = [name for idx, name in enumerate(all_names) if idx in sidxs]
    return snames, sidxs, exclude_set, imap


def build_sequence_imap_minmap(
    snames: list[str],
    *,
    min_sample_coverage: int | float,
    imap=None,
    minmap=None,
) -> tuple[dict[str, list[str]], dict[str, int | float]]:
    """Return the final IMAP/minmap after sample selection."""
    if not imap:
        return {"all": list(snames)}, {"all": int(min_sample_coverage)}

    if minmap is None:
        logger.info("No minmap provided. Assuming minimum one sample per population.")
        minmap = {pop: 1 for pop in imap.keys()}
    if set(minmap) != set(imap):
        raise IPyradError("imap and minmap keys must match.")

    included_names = set(snames)
    filtered_imap = {
        key: [name for name in names if name in included_names]
        for key, names in imap.items()
    }
    return filtered_imap, {key: minmap[key] for key in filtered_imap}


def sync_sequence_imap_after_sample_drop(
    snames: list[str],
    *,
    user_imap: bool,
    imap: dict[str, list[str]],
    minmap: dict[str, int | float],
    min_sample_coverage: int | float,
) -> tuple[dict[str, list[str]], dict[str, int | float]]:
    """Rebuild IMAP after dropping samples by max missingness."""
    if not user_imap:
        return {"all": list(snames)}, {"all": int(min_sample_coverage)}

    selected = set(snames)
    new_imap = {
        group: [name for name in names if name in selected]
        for group, names in imap.items()
    }
    empty = sorted(group for group, names in new_imap.items() if not names)
    if empty:
        raise IPyradError(
            "IMAP group(s) became empty after max_sample_missing filtering: "
            + ", ".join(empty)
        )
    return new_imap, minmap.copy()
