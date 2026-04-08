#!/usr/bin/env python

"""Write the sequence-alignment datasets for the final assemble HDF5 output.

`phy` stores the assembled alignment matrix, while `phymap` maps alignment
windows back to their genome coordinates using 0-based half-open indexing.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import h5py
import numpy as np
from loguru import logger
from .hdf5_utils import choose_hdf5_cache_settings
from .hdf5_utils import choose_unsigned_int_dtype
from .hdf5_utils import format_bytes
from .hdf5_utils import get_fai_values
from .hdf5_utils import get_retained_fai_rows
from .loci import get_retained_loci_manifest_path, iter_parse_loci, resolve_locus_for_output


@dataclass
class SeqsHdf5Writer:
    """Open HDF5 handles plus running counters for streamed sequence output."""

    io5: h5py.File
    phy: h5py.Dataset
    phymap: h5py.Dataset
    nsites: int = 0
    nloci: int = 0


def choose_chunk_cols(
    nsamples: int,
    itemsize: int,
    target_mb: int = 16,
    typical_window: int = 100_000,  # consider increasing
    min_cols: int = 8_192,
    max_cols: int = 262_144,
) -> int:
    """Choose an alignment-column chunk size that stays efficient to read."""
    # size-driven cols
    by_size = (target_mb * 1024 * 1024) // (nsamples * itemsize)
    cols = int(max(min_cols, min(by_size, max_cols)))
    return max(4096, (int((cols + typical_window) / 2) // 4096) * 4096)


def _ensure_chunk_capacity(
    chunkarr: np.ndarray,
    nsamples: int,
    required_cols: int,
) -> np.ndarray:
    """Grow the working alignment buffer when a locus exceeds the initial chunk."""
    if required_cols <= chunkarr.shape[1]:
        return chunkarr
    new_cols = max(chunkarr.shape[1] * 2, required_cols)
    grown = np.full((nsamples, new_cols), np.uint8(ord("N")), dtype=np.uint8)
    grown[:, :chunkarr.shape[1]] = chunkarr
    return grown


def open_seqs_hdf5_writer(
    *,
    name: str,
    outdir: Path,
    snames: list[str],
    reference: Path,
    loci_bed: Path,
) -> SeqsHdf5Writer:
    """Create a resizable HDF5 writer for streamed final-locus output."""
    seqs_database = outdir / f"{name}.hdf5"

    snames = sorted(snames)
    snames = ["assembly_reference_sequence"] + snames
    nsamples = len(snames)
    if loci_bed.exists() and loci_bed.stat().st_size:
        retained_rows = get_retained_fai_rows(reference, loci_bed)
    else:
        retained_rows = tuple()
    if not retained_rows and not (loci_bed.exists() and loci_bed.stat().st_size):
        scaffold_names = [str(value) for value in get_fai_values(reference, "scaffold")]
        scaffold_lengths = [int(value) for value in get_fai_values(reference, "length")]
    else:
        scaffold_names = [str(row[0]) for row in retained_rows]
        scaffold_lengths = [int(row[1]) for row in retained_rows]
    if loci_bed.exists() and loci_bed.stat().st_size and not scaffold_names:
        raise ValueError(f"No retained scaffolds found in final BED: {loci_bed}")

    chunk_size = choose_chunk_cols(
        nsamples,
        np.dtype(np.uint8).itemsize,
        target_mb=256,
    )
    phy_chunk_cols = max(1, chunk_size)
    phymap_chunk_rows = 4096
    map_dtype = choose_unsigned_int_dtype(
        max(
            len(scaffold_names) - 1,
            sum(scaffold_lengths),
            max(scaffold_lengths, default=0),
        )
    )

    kwargs = choose_hdf5_cache_settings()
    logger.debug(
        (
            "seqs writer config: phy_chunk_cols={}, phymap_chunk_rows={}, "
            "cache={}, cache_slots={}, map_dtype={}"
        ),
        phy_chunk_cols,
        phymap_chunk_rows,
        format_bytes(int(kwargs["rdcc_nbytes"])),
        kwargs["rdcc_nslots"],
        map_dtype.name,
    )
    io5 = h5py.File(seqs_database, "w", **kwargs)
    io5.attrs["version"] = 2.0
    io5.attrs["names"] = snames
    io5.attrs["reference"] = str(reference)
    io5.attrs["scaffold_lengths"] = scaffold_lengths
    io5.attrs["scaffold_names"] = scaffold_names

    phy = io5.create_dataset(
        "phy",
        shape=(nsamples, 0),
        maxshape=(nsamples, None),
        dtype=np.uint8,
        chunks=(nsamples, phy_chunk_cols),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
    )
    phy.attrs["nsites"] = 0

    phymap = io5.create_dataset(
        "phymap",
        shape=(0, 5),
        maxshape=(None, 5),
        dtype=map_dtype,
        chunks=(phymap_chunk_rows, 5),
    )
    phymap.attrs["indexing"] = [0, 0, 0, 0, 0]
    phymap.attrs["columns"] = ["scaff", "phy0", "phy1", "pos0", "pos1"]
    return SeqsHdf5Writer(io5=io5, phy=phy, phymap=phymap)


def append_seqs_hdf5_chunk(
    writer: SeqsHdf5Writer,
    chunkarr: np.ndarray,
    chunkmap: np.ndarray,
) -> None:
    """Append one ordered chunk of alignment columns and phymap rows."""
    if chunkarr.size:
        ncols = int(chunkarr.shape[1])
        writer.phy.resize((writer.phy.shape[0], writer.nsites + ncols))
        writer.phy[:, writer.nsites:writer.nsites + ncols] = chunkarr
        writer.nsites += ncols
        writer.phy.attrs["nsites"] = writer.nsites
    if chunkmap.size:
        nrows = int(chunkmap.shape[0])
        writer.phymap.resize((writer.nloci + nrows, writer.phymap.shape[1]))
        writer.phymap[writer.nloci:writer.nloci + nrows, :] = chunkmap
        writer.nloci += nrows


def finalize_seqs_hdf5_writer(
    writer: SeqsHdf5Writer,
    *,
    expected_nsites: int | None = None,
    expected_nloci: int | None = None,
    retained_scaffold_names: list[str] | None = None,
    retained_scaffold_lengths: list[int] | None = None,
) -> None:
    """Validate counts and compact scaffold metadata on a still-open HDF5 writer."""
    if expected_nsites is not None and writer.nsites != expected_nsites:
        raise ValueError(
            f"seqs writer filled {writer.nsites} sites but loci summary reported {expected_nsites}"
        )
    if expected_nloci is not None and writer.nloci != expected_nloci:
        raise ValueError(
            f"seqs writer filled {writer.nloci} phymap rows but loci summary reported {expected_nloci}"
        )
    if retained_scaffold_names is None:
        return
    if retained_scaffold_lengths is None:
        raise ValueError("retained_scaffold_lengths must be provided with retained_scaffold_names")

    old_names = [str(value) for value in writer.io5.attrs["scaffold_names"]]
    remap = {old_names.index(name): new_idx for new_idx, name in enumerate(retained_scaffold_names)}
    if writer.nloci:
        phymap = writer.phymap[:]
        phymap[:, 0] = np.array([remap[int(idx)] for idx in phymap[:, 0]], dtype=writer.phymap.dtype)
        writer.phymap[:] = phymap
    writer.io5.attrs["scaffold_names"] = retained_scaffold_names
    writer.io5.attrs["scaffold_lengths"] = retained_scaffold_lengths


def close_seqs_hdf5_writer(writer: SeqsHdf5Writer) -> None:
    """Close the streamed sequence HDF5 writer without additional validation."""
    writer.io5.close()


def write_seqs_hdf5(
    name: str,
    outdir: Path,
    tmpdir: Path,
    snames: list[str],
    reference: Path,
    loci_bed: Path,
    nsites_after_filtering: int,
    nloci_after_filtering: int,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float = 0.10,
) -> None:
    """Write the final alignment datasets into the combined assemble HDF5."""
    database = tmpdir / f"{name}.database.fa"
    snames = sorted(snames)
    snames = ["assembly_reference_sequence"] + snames
    retained_rows = get_retained_fai_rows(reference, loci_bed)
    scaffold_names = [str(row[0]) for row in retained_rows]
    scaff2idx = {name: idx for idx, name in enumerate(scaffold_names)}
    if nloci_after_filtering and not retained_rows:
        raise ValueError(f"No retained scaffolds found in final BED: {loci_bed}")
    map_dtype = choose_unsigned_int_dtype(
        max(
            len(scaffold_names) - 1,
            nsites_after_filtering,
            nloci_after_filtering,
            max((int(row[1]) for row in retained_rows), default=0),
        )
    )
    nsamples = len(snames)
    chunk_size = choose_chunk_cols(nsamples, np.dtype(np.uint8).itemsize, target_mb=256)
    params = (min_locus_sample_coverage, min_locus_trim_sample_coverage, min_locus_length, max_locus_hetero_frequency, max_locus_variant_frequency)
    retained_manifest = _load_retained_loci_manifest(name, tmpdir)
    iter_chunks = iter_super_matrix_chunks_from_database(
        snames,
        database,
        scaff2idx,
        chunk_size,
        map_dtype,
        retained_manifest,
        *params,
        max_sample_hetero_frequency,
    )
    writer = open_seqs_hdf5_writer(
        name=name,
        outdir=outdir,
        snames=snames[1:],
        reference=reference,
        loci_bed=loci_bed,
    )
    try:
        for chunkarr, chunkmap in iter_chunks:
            assert chunkarr.shape[0] == nsamples
            if chunkmap.size:
                chunkmap = chunkmap.astype(writer.phymap.dtype, copy=False)
            append_seqs_hdf5_chunk(writer, chunkarr, chunkmap)
        finalize_seqs_hdf5_writer(
            writer,
            expected_nsites=nsites_after_filtering,
            expected_nloci=nloci_after_filtering,
            retained_scaffold_names=[str(row[0]) for row in retained_rows],
            retained_scaffold_lengths=[int(row[1]) for row in retained_rows],
        )
    finally:
        close_seqs_hdf5_writer(writer)
    logger.debug("wrote seqs database to {}", outdir / f"{name}.hdf5")


def iter_super_matrix_chunks_from_database(
    snames: list[str],
    database: Path,
    scaff2idx: dict[str, int],
    chunk_size: int,
    map_dtype: np.dtype,
    retained_manifest: dict[str, tuple[str, set[str]]],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float,
):
    """Parse and filter loci and yield in chunks."""

    # NOTE: snames is already sorted with ref optionally on top.
    # the global order of names in the supermatrix
    sidxs = {j: i for (i, j) in enumerate(snames)}
    nsamples = len(sidxs)

    # create initial chunkarr and chunkmap
    chunkarr = np.full((nsamples, chunk_size * 2), np.uint8(ord('N')), dtype=np.uint8)
    chunkmap = []

    # iterate to fill chunkarr
    phypos = 0
    cursor = 0
    for oheader, ldict in iter_parse_loci(database):
        retained = retained_manifest.get(oheader)
        if retained is None:
            continue

        # trim and filter the locus
        args = (
            oheader,
            ldict,
            min_locus_sample_coverage,
            min_locus_trim_sample_coverage,
            min_locus_length,
            max_locus_hetero_frequency,
            max_locus_variant_frequency,
            max_sample_hetero_frequency,
        )
        header, tnames, tseqs, snparr, filters, _ = resolve_locus_for_output(
            *args,
            forced_masked_samples=retained[1],
        )
        if not sum(filters.values()):
            expected_header = retained[0]
            if header != expected_header:
                raise ValueError(
                    f"Retained locus header mismatch for {oheader}: expected {expected_header}, observed {header}"
                )

            # parse new trimmed locus header
            scaff, positions = header.strip().split(":")
            pos0, pos1 = [int(i) for i in positions.split("-")]
            try:
                scaff_idx = scaff2idx[scaff]
            except KeyError as exc:
                raise ValueError(
                    f"Retained locus scaffold {scaff!r} is absent from the final BED-derived scaffold table."
                ) from exc

            # enter data into chunkarr
            chunkarr = _ensure_chunk_capacity(chunkarr, nsamples, cursor + tseqs.shape[1])
            for loc_sidx, loc_name in enumerate(tnames):
                global_sidx = sidxs[loc_name]
                chunkarr[global_sidx, cursor:cursor + tseqs.shape[1]] = tseqs[loc_sidx]
            chunkmap.append((scaff_idx, phypos, phypos + tseqs.shape[1], pos0, pos1))
            cursor += tseqs.shape[1]
            phypos += tseqs.shape[1]

            # if the chunk is full, yield it and refresh
            if cursor >= chunk_size:
                # trim extra, fill empty to N, convert map to array
                chunkarr = chunkarr[:, :cursor]
                chunkmap = np.array(chunkmap, dtype=map_dtype)
                yield chunkarr, chunkmap

                # reset
                chunkarr = np.full((nsamples, chunk_size * 2), np.uint8(ord('N')), dtype=np.uint8)
                chunkmap = []
                cursor = 0  # reset cursor but not phypos

    # trim extra and convert the final map chunk to an array
    if cursor:
        chunkarr = chunkarr[:, :cursor]
        chunkmap = np.array(chunkmap, dtype=map_dtype)
        yield chunkarr, chunkmap


def _load_retained_loci_manifest(name: str, tmpdir: Path) -> dict[str, tuple[str, set[str]]]:
    """Load retained final loci plus per-locus masked sample names."""
    path = get_retained_loci_manifest_path(name, tmpdir)
    retained: dict[str, tuple[str, set[str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            masked = {
                item
                for item in row["masked_samples"].split(",")
                if item
            }
            retained[row["raw_header"]] = (row["final_header"], masked)
    return retained
