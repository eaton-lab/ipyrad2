#!/usr/bin/env python

"""
"""

from typing import Dict, Iterator, Tuple, List
from pathlib import Path
import h5py
import pandas as pd
import numpy as np

CHUNKSIZE = 50_000


def get_fai_values(reference: Path, key: str):
    """Returns the fai table from the reference as an array."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    columns = ['scaffold', 'length', 'sumsize', 'a', 'b']
    table = pd.read_csv(fai, names=columns, sep="\t")
    return table[key].values


def write_seqs_hdf5(
    name: str,
    outdir: Path,
    snames: List[str],
    reference: Path,
    exclude_reference: bool,
    nloci: int,
    nsites: int,
):
    """Write seqs h5 database from loci file.
    """
    # get global sorted names
    snames = sorted(snames)
    if not exclude_reference:
        snames = ["assembly_reference_sequence"] + snames
    nsamples = len(snames)

    # open h5 and write
    out_handle = outdir / f"{name}.seqs.hdf5"
    with h5py.File(out_handle, 'w') as io5:

        # store meta-data
        io5.attrs["names"] = snames
        io5.attrs["version"] = 2.0
        io5.attrs["nsites"] = int(nsites)
        io5.attrs["reference"] = str(reference)
        # Store phymap meta-data as list of strings not bytes
        io5.attrs["scaffold_lengths"] = [int(i) for i in get_fai_values(reference, "length")]
        io5.attrs["scaffold_names"] = [str(i) for i in get_fai_values(reference, "scaffold")]
        # store ordered names and column labels
        io5.attrs["indexing"] = [0, 0, 0, 0, 0]
        io5.attrs["columns"] = ["scaff", "phy0", "phy1", "pos0", "pos1"]

        # create the seq and seq positions arr datasets.
        _ = io5.create_dataset(
            name="phy",
            shape=(nsamples, nsites),
            dtype=np.uint8,
        )
        _ = io5.create_dataset(
            name="phymap",
            shape=(nloci, 5),
            dtype=np.uint64,
        )

        # start filling the supermatrix
        start_pos = 0
        start_arr = 0
        start_map = 0
        loci_file = outdir / f"{name}.loci.txt"
        for chunkarr, chunkmap in iter_supermatrix_chunks(loci_file, reference, snames):
            # advance map_chunk
            chunkmap[:, [1, 2]] += start_pos

            # get end positions
            end_arr = start_arr + chunkarr.shape[1]
            end_map = start_map + chunkmap.shape[0]

            # insert to dataset array
            io5["phy"][:, start_arr:end_arr] = chunkarr
            io5["phymap"][start_map:end_map] = chunkmap

            # increment start positions
            start_arr += chunkarr.shape[1]
            start_map += chunkmap.shape[0]
            start_pos = int(chunkarr[-1, 2])


def iter_supermatrix_chunks(loci_file: Path, reference: Path, snames: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Generate (array, array) of CHUNKSIZE from parsed loci file data.
    """
    # the global order of names in the supermatrix
    sidxs = {j: i for (i, j) in enumerate(snames)}
    nsamples = len(sidxs)

    # scaffnames
    scaff_names = list(get_fai_values(reference, "scaffold"))

    # create initial chunkarr and chunkmap
    chunkarr = np.zeros((nsamples, CHUNKSIZE * 2), dtype=np.uint8)
    chunkmap = []

    # iterate to fill chunkarr
    phypos = 0
    cursor = 0
    for locus_dict, header in iter_loci(loci_file):

        _ = locus_dict.pop("snpstring").decode()
        snames = list(locus_dict.keys())
        seqs = [list(i) for i in locus_dict.values()]
        arr = np.array(seqs, dtype=np.uint8)

        # parse header
        fldx, scaff, pos0, pos1 = header
        scaff_idx = scaff_names.index(scaff)

        # enter data into chunkarr
        for loc_sidx, sname in enumerate(snames):
            global_sidx = sidxs[sname]
            chunkarr[global_sidx, cursor:cursor + arr.shape[1]] = arr[loc_sidx]
        chunkmap.append((scaff_idx, phypos, phypos + arr.shape[1], pos0, pos1))
        cursor += arr.shape[1]
        phypos += arr.shape[1]

        # if the chunk is full, yield it and refresh
        if cursor >= CHUNKSIZE:
            # trim extra, fill empty to N, convert map to array
            chunkarr = chunkarr[:, :cursor]
            chunkarr[chunkarr == 0] = 78
            chunkmap = np.array(chunkmap, dtype=np.uint64)
            yield chunkarr, chunkmap

            # reset
            chunkarr = np.zeros((nsamples, CHUNKSIZE * 2), dtype=np.uint8)
            chunkmap = []
            cursor = 0  # reset cursor but not phypos

    # trim extra, fill empty to N, convert map to array
    if cursor:
        chunkarr = chunkarr[:, :cursor]
        chunkarr[chunkarr == 0] = 78
        chunkmap = np.array(chunkmap, dtype=np.uint64)
        yield chunkarr, chunkmap
        # missing_cells += np.sum((chunkarr == 45) | (chunkarr == 78))


def iter_loci(loci_file: Path) -> Iterator[Tuple[Dict[str, str], Tuple[int, str, int, int]]]:
    """Yields loci from each ordered .loci file until empty.
    """
    with open(loci_file, 'r', encoding="utf-8") as indata:

        # get pad length from first line then reset
        snp_pad = indata.readline().rfind(" ") + 1
        indata.seek(0)

        # iterate over loci and yield one at a time.
        names_to_seqs = {}
        for line in indata:
            # within a locus, fill the data.
            if not line.startswith("//"):
                name, seq = line.split()
                names_to_seqs[name] = bytes(seq, "utf-8")
            else:
                names_to_seqs['snpstring'] = bytes(line[snp_pad:].split("|")[0], "utf-8")

                # parse reference position from snpstring
                line_chunks = line.split("|")
                chrom_int, chrom_name, pos = line_chunks[1].split(":")
                chrom_int = int(chrom_int)
                pos0, pos1 = [int(i) for i in pos.split("-")]

                # end of locus, yield the dict.
                yield names_to_seqs, (chrom_int, chrom_name, pos0, pos1)
                names_to_seqs = {}



if __name__ == "__main__":

    # loci_file = data.stepdir / f"{data.name}.loci.txt"
    loci_file = Path("/tmp/OUT_klmnop/") / "assembly.loci.txt"

    write_seqs_hdf5(
        name="TEST",
        outdir=Path("/tmp/OUT_klmnop"),
        snames=...,
        reference=...,
        exclude_reference=False,
        nloci=4759,
        nsites=911965,
    )
