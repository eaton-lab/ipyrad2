#!/usr/bin/env python

"""Store sequence data in h5 for fast extraction.

Q and A to remember why we do it this way
-----------------------------------------
Q: Why store the full sequences including PE inserts?
A: Even though it is wasteful of space, we want the H5 database
   to match the coordinates of the loci assembled given the parameters
   set during assembly. This ensures there is a 1-to-1 match between
   the length of a locus and the length of the window in the phymap,
   and allows extracting loci just like in the loci file.
   The window can still be filtered later using tools like wex.

Q: Why parse and filter from the database.fa, rather than loading
   the .loci file?
A: It is pretty fast, and this way we can fill the h5 databases in
   parallel while also writing the .loci file.

Q: How does indexing work in the phymap?
A: It stores the mapping of phy coordinates to genome positions using
   Python (and bedtools) 0-based half-open indexing. Tools such as wex
   which read this phymap will convert the coordinates when given
   samtools region inputs. E.g., Chr1:100-250 will extract the phy
   corresponding to 99-250 in phymap.

PHYMAP format
-------------
scaff    phy0    phy1   pos0    pos1
    0       0     100   1000    1100
    0     100     200   3200    3300
    ...

"""

from typing import List
from pathlib import Path
import h5py
import pandas as pd
import numpy as np
from loguru import logger
from .loci import iter_parse_loci, filter_trim_locus


def get_fai_values(reference: Path, key: str) -> np.ndarray:
    """Returns the fai table from the reference as an array."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    columns = ['scaffold', 'length', 'sumsize', 'a', 'b']
    table = pd.read_csv(fai, names=columns, sep="\t")
    return table[key].values


def choose_chunk_cols(
    nsamples: int,
    itemsize: int,
    target_mb: int = 16,
    typical_window: int = 100_000,  # consider increasing
    min_cols: int = 8_192,
    max_cols: int = 512_000
) -> int:
    """Get chunksize to make reading all rows and thousands to
    tens/hundreds of thousands of columns efficient. We will set
    chunk columns to this size so most reads touch 1–3 chunks.
    """
    # size-driven cols
    by_size = (target_mb * 1024 * 1024) // (nsamples * itemsize)
    return by_size
    logger.warning(by_size)
    cols = int(max(min_cols, min(by_size, max_cols)))
    logger.warning(cols)
    # bias toward typical window (round to nearest 4096)
    cols = max(4096, (int((cols + typical_window) / 2) // 4096) * 4096)
    logger.warning(cols)
    return cols


def write_seqs_hdf5(
    name: str,
    outdir: Path,
    snames: List[str],
    reference: Path,
    exclude_reference: bool,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
):
    """Write seqs h5 database from loci file.
    """
    # paths
    database = outdir / "tmpdir" / f"{name}.database.fa"
    seqs_database = outdir / f"{name}.hdf5"

    # get global sorted names with reference sequence on top
    snames = sorted(snames)
    snames = ["assembly_reference_sequence"] + snames
    nsamples = len(snames)

    # get optimal chunk size
    chunk_size = choose_chunk_cols(
        nsamples,
        np.dtype(np.uint8).itemsize,
        target_mb=256,
        # typical_window=100_000,
    )

    # get the data generator
    params = (min_locus_sample_coverage, min_locus_trim_sample_coverage, min_locus_length, max_locus_hetero_frequency, max_locus_variant_frequency)
    iter_chunks = iter_super_matrix_chunks_from_database(snames, database, reference, chunk_size, *params)

    # open H5: 512 MB raw data chunk cache, many hash slots reduces collisions
    maps = []
    kwargs = dict(libver="latest", rdcc_nbytes=512*1024*1024, rdcc_nslots=2_000_003)
    with h5py.File(seqs_database, "w", **kwargs) as io5:

        # database metadata.
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = snames
        io5.attrs["reference"] = str(reference)
        io5.attrs["scaffold_lengths"] = [int(i) for i in get_fai_values(reference, "length")]
        io5.attrs["scaffold_names"] = [str(i) for i in get_fai_values(reference, "scaffold")]

        phy = io5.create_dataset(
            "phy",
            shape=(nsamples, 0),
            maxshape=(nsamples, None),
            dtype=np.uint8,
            chunks=(nsamples, chunk_size),
            compression="gzip",          # write-once → favor space
            compression_opts=4,          # 2–6 is a good speed/ratio range
            shuffle=True,                # improves gzip on byte-y data
        )

        # Append column batches
        col_total = 0
        for chunkarr, chunkmap in iter_chunks:
            assert chunkarr.shape[0] == nsamples
            m = chunkarr.shape[1]
            maps.append(chunkmap)

            # Optional pre-grow to amortize resizes (e.g., grow in big steps)
            need = col_total + m
            cap  = phy.shape[1]
            if need > cap:
                # grow to next multiple of chunk_cols to avoid frequent resizes
                new_cap = ((need + chunk_size - 1) // chunk_size) * chunk_size
                phy.resize(new_cap, axis=1)
                logger.debug(f"new cap = {new_cap}")
            phy[:, col_total:col_total + m] = chunkarr
            col_total += m

        # Optionally shrink tail if we overgrew
        phy.attrs["nsites"] = col_total
        phy.resize(col_total, axis=1)

        # store the phymap -------------------------------------------
        phymap = io5.create_dataset(
            name="phymap",
            data=np.concatenate(maps, axis=0),
        )
        phymap.attrs["indexing"] = [0, 0, 0, 0, 0]
        phymap.attrs["columns"] = ["scaff", "phy0", "phy1", "pos0", "pos1"]
    logger.debug(f"wrote seqs database to {seqs_database}")


def iter_super_matrix_chunks_from_database(
    snames: List[str],
    database: Path,
    reference: Path,
    chunk_size: int,
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

    # scaffnames
    scaff_names = list(get_fai_values(reference, "scaffold"))

    # create initial chunkarr and chunkmap
    chunkarr = np.full((nsamples, chunk_size * 2), np.uint8(ord('N')), dtype=np.uint8)
    chunkmap = []

    # iterate to fill chunkarr
    phypos = 0
    cursor = 0
    for oheader, ldict in iter_parse_loci(database):
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
            scaff_idx = scaff_names.index(scaff)

            # enter data into chunkarr
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
                # chunkarr[chunkarr == 0] = 78
                chunkmap = np.array(chunkmap, dtype=np.uint64)
                yield chunkarr, chunkmap

                # reset
                chunkarr = np.full((nsamples, chunk_size * 2), np.uint8(ord('N')), dtype=np.uint8)
                # chunkarr = np.zeros((nsamples, CHUNKSIZE * 2), dtype=np.uint8)
                chunkmap = []
                cursor = 0  # reset cursor but not phypos

    # trim extra, fill empty to N, convert map to array
    if cursor:
        chunkarr = chunkarr[:, :cursor]
        # chunkarr[chunkarr == 0] = 78
        chunkmap = np.array(chunkmap, dtype=np.uint64)
        yield chunkarr, chunkmap



if __name__ == "__main__":

    pass

    # loci_file = data.stepdir / f"{data.name}.loci.txt"
    DIR = Path("/home/deren/Documents/ipyrad-tests/OUT/")
    REF = Path("/home/deren/Documents/tools/ipyrad/tests/ipsimdata/pairddrad_example_genome.fa")

    SNAMES = [
        "1A_0_R", "1B_0_R", "1C_0_R", "1D_0_R",
        "2E_0_R", "2F_0_R", "2G_0_R", "2H_0_R",
        "3I_0_R", "3J_0_R", "3K_0_R", "3L_0_R",
    ]

    write_seqs_hdf5(
        name="assembly",
        outdir=DIR,
        snames=SNAMES,
        reference=REF,
        exclude_reference=False,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=4,
        min_locus_length=35,
        max_locus_hetero_frequency=0.3,
        max_locus_variant_frequency=1.0,
    )
