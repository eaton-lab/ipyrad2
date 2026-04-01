#!/usr/bin/env python

"""Write the sequence-alignment datasets for the final assemble HDF5 output.

`phy` stores the assembled alignment matrix, while `phymap` maps alignment
windows back to their genome coordinates using 0-based half-open indexing.
"""

from __future__ import annotations

from pathlib import Path
import h5py
import numpy as np
from loguru import logger
from .hdf5_utils import choose_hdf5_cache_settings
from .hdf5_utils import choose_unsigned_int_dtype
from .hdf5_utils import format_bytes
from .hdf5_utils import get_fai_values
from .loci import iter_parse_loci, filter_trim_locus


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


def write_seqs_hdf5(
    name: str,
    outdir: Path,
    tmpdir: Path,
    snames: list[str],
    reference: Path,
    nsites_after_filtering: int,
    nloci_after_filtering: int,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
) -> None:
    """Write the final alignment datasets into the combined assemble HDF5."""
    # paths
    database = tmpdir / f"{name}.database.fa"
    seqs_database = outdir / f"{name}.hdf5"

    # get global sorted names with reference sequence on top
    snames = sorted(snames)
    snames = ["assembly_reference_sequence"] + snames
    nsamples = len(snames)
    scaffold_names = [str(i) for i in get_fai_values(reference, "scaffold")]
    scaffold_lengths = [int(i) for i in get_fai_values(reference, "length")]

    # get optimal chunk size
    chunk_size = choose_chunk_cols(
        nsamples,
        np.dtype(np.uint8).itemsize,
        target_mb=256,
    )
    phy_chunk_cols = min(max(1, nsites_after_filtering), chunk_size)
    phymap_chunk_rows = min(max(1, nloci_after_filtering), 4096)
    map_dtype = choose_unsigned_int_dtype(
        max(
            len(scaffold_names) - 1,
            nsites_after_filtering,
            nloci_after_filtering,
            max(scaffold_lengths, default=0),
        )
    )

    # get the data generator
    params = (min_locus_sample_coverage, min_locus_trim_sample_coverage, min_locus_length, max_locus_hetero_frequency, max_locus_variant_frequency)
    iter_chunks = iter_super_matrix_chunks_from_database(
        snames,
        database,
        reference,
        chunk_size,
        map_dtype,
        *params,
    )

    kwargs = choose_hdf5_cache_settings()
    logger.debug(
        (
            "seqs writer config: phy_chunk_cols={}, phymap_chunk_rows={}, "
            "cache={}, cache_slots={}, map_dtype={}, nsites={}, nloci={}"
        ),
        phy_chunk_cols,
        phymap_chunk_rows,
        format_bytes(int(kwargs["rdcc_nbytes"])),
        kwargs["rdcc_nslots"],
        map_dtype.name,
        nsites_after_filtering,
        nloci_after_filtering,
    )
    with h5py.File(seqs_database, "w", **kwargs) as io5:

        # database metadata.
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = snames
        io5.attrs["reference"] = str(reference)
        io5.attrs["scaffold_lengths"] = scaffold_lengths
        io5.attrs["scaffold_names"] = scaffold_names

        phy = io5.create_dataset(
            "phy",
            shape=(nsamples, nsites_after_filtering),
            dtype=np.uint8,
            chunks=(nsamples, phy_chunk_cols),
            compression="gzip",          # write-once → favor space
            compression_opts=4,          # 2–6 is a good speed/ratio range
            shuffle=True,                # improves gzip on byte-y data
        )
        phy.attrs["nsites"] = int(nsites_after_filtering)

        phymap = io5.create_dataset(
            "phymap",
            shape=(nloci_after_filtering, 5),
            dtype=map_dtype,
            chunks=(phymap_chunk_rows, 5),
        )
        phymap.attrs["indexing"] = [0, 0, 0, 0, 0]
        phymap.attrs["columns"] = ["scaff", "phy0", "phy1", "pos0", "pos1"]

        # Fill the preallocated datasets in streaming batches. The final locus
        # pass already told us exactly how many retained sites and loci exist,
        # so we do not need repeated HDF5 growth during writing.
        col_total = 0
        phymap_total = 0
        for chunkarr, chunkmap in iter_chunks:
            assert chunkarr.shape[0] == nsamples
            m = chunkarr.shape[1]
            phy[:, col_total:col_total + m] = chunkarr
            col_total += m

            if chunkmap.size:
                nrows = chunkmap.shape[0]
                phymap[phymap_total:phymap_total + nrows, :] = chunkmap
                phymap_total += nrows

    if col_total != nsites_after_filtering:
        raise ValueError(
            f"seqs writer filled {col_total} sites but loci summary reported {nsites_after_filtering}"
        )
    if phymap_total != nloci_after_filtering:
        raise ValueError(
            f"seqs writer filled {phymap_total} phymap rows but loci summary reported {nloci_after_filtering}"
        )
    logger.debug(f"wrote seqs database to {seqs_database}")


def iter_super_matrix_chunks_from_database(
    snames: list[str],
    database: Path,
    reference: Path,
    chunk_size: int,
    map_dtype: np.dtype,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
):
    """Parse and filter loci and yield in chunks."""

    # NOTE: snames is already sorted with ref optionally on top.
    # the global order of names in the supermatrix
    sidxs = {j: i for (i, j) in enumerate(snames)}
    nsamples = len(sidxs)

    # Precompute scaffold indexes once so locus parsing does not repeatedly
    # scan the scaffold-name list for every retained locus.
    scaff_names = list(get_fai_values(reference, "scaffold"))
    scaff2idx = {str(name): idx for idx, name in enumerate(scaff_names)}

    # create initial chunkarr and chunkmap
    chunkarr = np.full((nsamples, chunk_size * 2), np.uint8(ord('N')), dtype=np.uint8)
    chunkmap = []

    # iterate to fill chunkarr
    phypos = 0
    cursor = 0
    for oheader, ldict in iter_parse_loci(database):
        # trim and filter the locus
        args = (
            oheader,
            ldict,
            min_locus_sample_coverage,
            min_locus_trim_sample_coverage,
            min_locus_length,
            max_locus_hetero_frequency,
            max_locus_variant_frequency,
        )
        header, tnames, tseqs, snparr, filters, _ = filter_trim_locus(*args)
        if not sum(filters.values()):

            # parse new trimmed locus header
            scaff, positions = header.strip().split(":")
            pos0, pos1 = [int(i) for i in positions.split("-")]
            scaff_idx = scaff2idx[scaff]

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
