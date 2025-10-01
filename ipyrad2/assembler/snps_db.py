#!/usr/bin/env python

"""Write a SNPs database file.

genos
-----
(nsnps, nsamples, 2)


"""

from __future__ import annotations
from typing import Dict, Iterator, Tuple, Union, List, Optional
from pathlib import Path
import gzip
import h5py
import pandas as pd
import numpy as np
from loguru import logger


def get_fai_values(reference: Path, key: str) -> np.ndarray:
    """Returns the fai table from the reference as an array."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    columns = ['scaffold', 'length', 'sumsize', 'a', 'b']
    table = pd.read_csv(fai, names=columns, sep="\t")
    return table[key].values


def write_snps_hdf5_2(
    name: str,
    outdir: Path,
    snames: List[str],
    reference: Path,
    exclude_reference: bool,
    vcf_path: Path,
    loci_bed: Path,
):
    """Write seqs h5 database from loci file.

    Parameters
    -----------
    vcf_path: Path
        vcf can be gzipped or not.
    """
    # paths
    database = outdir / f"{name}.hdf5"

    # get global sorted names
    snames = sorted(snames)
    if not exclude_reference:
        snames = ["assembly_reference_sequence"] + snames
    nsamples = len(snames)

    # open H5: 256 MB raw data chunk cache, many hash slots reduces collisions
    maps = []
    kwargs = dict(rdcc_nbytes=256*1024*1024, rdcc_nslots=1_000_003)
    with h5py.File(database, 'w', **kwargs) as io5:

        # database metadata.
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = snames
        io5.attrs["reference"] = str(reference)
        io5.attrs["scaffold_lengths"] = [int(i) for i in get_fai_values(reference, "length")]
        io5.attrs["scaffold_names"] = [str(i) for i in get_fai_values(reference, "scaffold")]

        # 2D array of the variant sites
        _ = io5.create_dataset(
            "snps",
            shape=(nsamples, 0),
            maxshape=(nsamples, None),
            dtype=np.uint8,
            chunks=(nsamples, CHUNKSIZE),
            compression="lzf",      # or compression="gzip", compression_opts=4
            shuffle=True
        )

        # 2D array of scaff and loc positions 0-indexed
        snpsmap = io5.create_dataset(
            name="snpsmap",
            shape=(nsnps, 5),
            dtype=np.uint64,
        )
        snpsmap.attrs["indexing"] = [0, 0, 0, 0, 0]
        snpsmap.attrs["columns"] = ["loc", "loc_idx", "loc_pos", "scaff", "pos"]

        # stores genotype calls (ALTS) as (0,0, 0,1, 0,2, etc) where
        # for reference data 0 is the REF allele, else it is the most
        # common allele in denovo.
        _ = io5.create_dataset(
            name="genos",
            shape=(nsnps, nsamples, 2),
            dtype=np.uint8,
        )

        # fill the arrays


    logger.debug(f"wrote snps database to {database}")



def load_bed_index_nonoverlap(
    bed_path: Union[str, Path],
    scaffold_order: Optional[list[str]] = None,
) -> Tuple[Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]], Dict[str, int]]:
    """
    Load a non-overlapping BED (optionally .gz) into per-chrom arrays.

    Returns:
      - bed_index: chrom -> (starts[int64], ends[int64], orig_bed_row_index[int64])
      - scaff2idx: chrom -> 0-based scaffold index (order from scaffold_order or BED encounter order)
    """
    bed_path = str(bed_path)
    compression = "gzip" if bed_path.endswith(".gz") else None

    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        compression=compression,
        dtype={0: str, 1: "Int64", 2: "Int64"},
        engine="c",
        na_filter=False,
    )
    if df.shape[1] < 3:
        raise ValueError("BED must have at least 3 columns (chrom, start, end).")

    df = df[[0, 1, 2]].rename(columns={0: "chrom", 1: "start", 2: "end"})
    df = df[df["start"].notna() & df["end"].notna() & (df["start"] >= 0) & (df["end"] > df["start"])]

    # Keep original row index so callers can refer back to the input BED
    df["bed_idx"] = df.index.to_numpy()

    bed_index: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    chrom_order: list[str] = []

    for chrom, sub in df.groupby("chrom", sort=False):  # preserve encounter order
        chrom_order.append(chrom)
        sub = sub.sort_values("start", kind="mergesort")  # stable
        starts = sub["start"].to_numpy(dtype=np.int64, copy=True)
        ends   = sub["end"].to_numpy(dtype=np.int64, copy=True)
        idxs   = sub["bed_idx"].to_numpy(dtype=np.int64, copy=True)

        # Non-overlap sanity check
        if starts.size and np.any(starts[1:] < ends[:-1]):
            raise ValueError(f"BED intervals overlap on {chrom}, but non-overlap was assumed.")

        bed_index[chrom] = (starts, ends, idxs)

    if scaffold_order is not None:
        scaff2idx = {name: i for i, name in enumerate(scaffold_order)}
    else:
        scaff2idx = {name: i for i, name in enumerate(chrom_order)}

    return bed_index, scaff2idx



def iter_vcf_filtered_snps_with_bed(
    vcf_path: Union[str, Path],
    bed_path: Union[str, Path],
    scaffold_order: Optional[list[str]] = None,
) -> Iterator[Tuple[int, int, int, int, int, str, str, str, str]]:
    """
    Iterate a VCF(.gz) and yield ONLY variants with FILTER==PASS and NOT having INDEL in INFO, as:
      (bed_idx, var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name, REF, ALT, QUAL)

    - VCF POS is 1-based; convert to pos0 (0-based).
    - BED intervals are 0-based, half-open [start, end), and non-overlapping.
    - If any such variant is not contained in the BED, raises ValueError immediately.
    """
    bed_index, scaff2idx = load_bed_index_nonoverlap(bed_path, scaffold_order)

    # per-interval counters (for var_idx_in_bed)
    counters: Dict[str, np.ndarray] = {
        chrom: np.zeros(len(bed_index[chrom][0]), dtype=np.int64)
        for chrom in bed_index
    }

    vcf_path = str(vcf_path)
    _open = gzip.open if vcf_path.endswith(".gz") else open
    with _open(vcf_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue  # need at least CHROM POS ID REF ALT QUAL FILTER INFO
            chrom, spos, _id, REF, ALT, QUAL, FILTER, INFO = parts[:8]

            # keep only FILTER==PASS
            if FILTER != "PASS":
                continue
            # exclude INDELs: presence of literal "INDEL" flag in INFO
            # (INFO is ';' separated; flags appear as bare tokens)
            if "INDEL" in INFO.split(";"):
                continue

            try:
                pos = int(spos)
            except ValueError:
                continue
            pos0 = pos - 1
            scaff_name = chrom
            scaff_idx = scaff2idx.get(scaff_name, -1)

            tup = bed_index.get(chrom)
            if tup is None:
                raise ValueError(f"Variant {chrom}:{pos} not covered by BED (chrom absent).")

            starts, ends, bedix = tup
            i = np.searchsorted(starts, pos0, side="right") - 1
            if i >= 0 and pos0 < ends[i]:
                var_idx_in_bed = int(counters[chrom][i])
                counters[chrom][i] += 1
                offset_in_bed = int(pos0 - starts[i])
                yield int(bedix[i]), var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name, REF, ALT, QUAL
            else:
                # not contained in any interval on this chrom
                # Find nearest interval for better error message
                left_end = ends[i] if i >= 0 else None
                right_start = starts[i+1] if (i+1) < len(starts) else None
                raise ValueError(
                    f"Variant {chrom}:{pos} (pos0={pos0}) not covered by BED. "
                    f"Nearest interval ends at {left_end} and next starts at {right_start}."
                )


if __name__ == "__main__":

    VCF = "/home/deren/Documents/ipyrad-tests/Ama-out/test.vcf.gz"
    BED = "/home/deren/Documents/ipyrad-tests/Ama-out/test.bed"
    bdict = load_bed_index_nonoverlap(BED)
    # print(bdict["A_tuberculatus_Chr01"])

    # ii = iter_vcf_with_bed(VCF, bdict, False)


    ii = iter_vcf_filtered_snps_with_bed(VCF, BED)
    (bed_idx, var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name) = next(ii)
    print((bed_idx, var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name))


if __name__ == "__main___":
    pass
