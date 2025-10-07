#!/usr/bin/env python

"""Write a SNPs database file.

snpsmap
-------
description: The map of SNP positions on RAD loci and genome scaffolds.
Loci, scaffolds, and positions are all stored 0-indexed.
dtype: np.uint8
shape: (nsnps, 5)
attrs["columns"]: ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
attrs["indexing"]: [0, 0, 0, 0, 0]

genos
-----
description: The ordered diploid genotype call of the ordered samples at
every SNP position. If call is "1/1" and
dtype: np.uint8
shape: (nsamples, nsnps, 3)
example:
>>> genos[0, :, 2] = [84, 84, 84, 87, 78]  # get one column
>>> genos[:, 0] != genos[:, 1]             # get hetero sites

reference
---------
description: The ordered REF allele at every SNP position.
dtype: np.uint8
shape: (nsnps,)
"""

from __future__ import annotations
from typing import Dict, Iterator, Tuple, Union, List, Optional
from pathlib import Path
import gzip
import h5py
import pandas as pd
import numpy as np
from loguru import logger


# IUPAC ambiguity code for heterozygous SNPs (unordered pairs)
_IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def get_fai_values(reference: Path, key: str) -> np.ndarray:
    """Returns the fai table from the reference as an array."""
    fai = reference.with_suffix(reference.suffix + ".fai")
    columns = ['scaffold', 'length', 'sumsize', 'a', 'b']
    table = pd.read_csv(fai, names=columns, sep="\t")
    return table[key].values


# UPDATE THIS: CHATGPT MADE. CHUNKS SHOULD JUST USE TARGET MB. MAYBE 512MB
def _choose_chunk_snps(
    nsamples: int,
    target_mb: int = 16,
    typical_window: int = 100_000,
    min_snps: int = 8_192,
    max_snps: int = 262_144,
) -> int:
    """#SNPs per chunk so each genos chunk ≈ target_mb and aligns with your window."""
    bytes_per_snp = max(1, nsamples * 3)             # uint8 * 3 planes
    by_size = int((target_mb * 1024 * 1024) // bytes_per_snp)
    snps = min(max(by_size, min_snps), max_snps)
    # bias toward your window (don’t exceed it much)
    snps = min(snps, typical_window)
    # round to nice boundary
    return max(4096, (snps // 4096) * 4096)


def write_snps_hdf5(
    name: str,
    outdir: Path,
    snames: List[str],
    reference: Path,
    scaffold_order: Optional[List[str]] = None,
):
    """Stream VCF→HDF5 with read-optimized chunking."""
    # paths
    database = outdir / f"{name}.hdf5"
    vcf_path = outdir / f"{name}.vcf.gz"
    loci_bed = outdir / f"{name}.bed"

    # sorted names
    snames = sorted(snames)
    nsamples = len(snames)

    # pick chunk size along SNP axis (tuned for read-many)
    chunk_snps = _choose_chunk_snps(nsamples, target_mb=16, typical_window=100_000)

    # HDF5 file/open options
    kwargs = dict(libver="latest", rdcc_nbytes=512*1024*1024, rdcc_nslots=2_000_003)
    with h5py.File(database, "a", **kwargs) as io5:
        # ---- metadata ----
        scaff_names = [str(i) for i in get_fai_values(reference, "scaffold")]
        scaff_lens  = [int(i) for i in get_fai_values(reference, "length")]

        io5.attrs["version"] = 2.0
        io5.attrs["names"] = snames
        io5.attrs["reference"] = str(reference)
        io5.attrs["scaffold_names"] = scaff_names
        io5.attrs["scaffold_lengths"] = scaff_lens

        # ---- datasets (extendable along SNP axis) ----
        # SNP map: (n_snps, 5)
        snpsmap = io5.create_dataset(
            "snpsmap",
            shape=(0, 5),
            maxshape=(None, 5),
            dtype=np.uint64,
            chunks=(chunk_snps, 5),
            compression="gzip", compression_opts=4, shuffle=True
        )
        snpsmap.attrs["columns"] = ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
        snpsmap.attrs["indexing"] = [0, 0, 0, 0, 0]

        # Genotypes: (nsamples, n_snps, 3)  ← note fixed sample axis
        genos = io5.create_dataset(
            "genos",
            shape=(nsamples, 0, 3),
            maxshape=(nsamples, None, 3),
            dtype=np.uint8,
            chunks=(nsamples, chunk_snps, 3),
            compression="gzip", compression_opts=4, shuffle=True
        )

        # Reference ord per SNP: (n_snps,)
        reference_ord = io5.create_dataset(
            "reference",
            shape=(0,),
            maxshape=(None,),
            dtype=np.uint8,
            chunks=(chunk_snps,),
            compression="gzip", compression_opts=4, shuffle=True
        )

        # --- Streaming buffers (one chunk per flush) ---
        buf_map = np.empty((chunk_snps, 5), dtype=np.uint64)
        buf_gen = np.empty((nsamples, chunk_snps, 3), dtype=np.uint8)
        buf_ref = np.empty((chunk_snps,), dtype=np.uint8)
        fill = 0
        total = 0
        N_ord = np.uint8(ord("N"))

        def flush(n: int):
            nonlocal total
            if n == 0:
                return
            new_total = total + n
            snpsmap.resize((new_total, 5))
            genos.resize((nsamples, new_total, 3))
            reference_ord.resize((new_total,))
            # write slices
            snpsmap[total:new_total, :] = buf_map[:n, :]
            genos[:, total:new_total, :] = buf_gen[:, :n, :]
            reference_ord[total:new_total] = buf_ref[:n]
            total = new_total

        # drain generator
        it = iter_vcf_filtered_snps_with_bed(vcf_path, loci_bed, snames, scaffold_order)
        for bed_idx, var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name, REF, ALT, QUAL, GT in it:
            # map row
            buf_map[fill, 0] = np.uint64(bed_idx)
            buf_map[fill, 1] = np.uint64(var_idx_in_bed)
            buf_map[fill, 2] = np.uint64(offset_in_bed)
            buf_map[fill, 3] = np.uint64(scaff_idx)
            buf_map[fill, 4] = np.uint64(pos0)

            # genotypes come as (nsamples,3) → place into current SNP column
            buf_gen[:, fill, :] = GT

            # reference ord
            buf_ref[fill] = np.uint8(ord(REF)) if REF in ("A","C","G","T") else N_ord

            fill += 1
            if fill == chunk_snps:
                flush(fill)
                fill = 0
        flush(fill)
        io5.attrs["nsnps"] = int(total)
    logger.debug(f"wrote snps dataset to {database} (nsnps={total:,})")


def load_bed_index_nonoverlap(
    bed_path: Union[str, Path],
    scaffold_order: Optional[List[str]] = None,
) -> Tuple[Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]], Dict[str, int]]:
    """Load a non-overlapping BED into per-chrom arrays and a scaffold index map."""
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
    df["bed_idx"] = df.index.to_numpy()

    bed_index: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    chrom_order: List[str] = []

    for chrom, sub in df.groupby("chrom", sort=False):  # preserve encounter order
        chrom_order.append(chrom)
        sub = sub.sort_values("start", kind="mergesort")
        starts = sub["start"].to_numpy(dtype=np.int64, copy=True)
        ends   = sub["end"].to_numpy(dtype=np.int64, copy=True)
        idxs   = sub["bed_idx"].to_numpy(dtype=np.int64, copy=True)
        if starts.size and np.any(starts[1:] < ends[:-1]):
            raise ValueError(f"BED intervals overlap on {chrom}, but non-overlap was assumed.")
        bed_index[chrom] = (starts, ends, idxs)

    scaff2idx = {name: i for i, name in enumerate(scaffold_order)} if scaffold_order \
                else {name: i for i, name in enumerate(chrom_order)}
    return bed_index, scaff2idx


def _parse_gt_to_uint8(gt_field: str) -> np.ndarray:
    """
    Parse a VCF GT subfield into two uint8 allele indexes.
    - REF=0, first ALT=1, etc.
    - Missing allele -> 255
    """
    out = np.array([255, 255], dtype=np.uint8)
    if not gt_field or gt_field == ".":
        return out
    sep = '/' if '/' in gt_field else ('|' if '|' in gt_field else None)
    if sep is None:
        # single allele (haploid-like)
        try:
            v = int(gt_field)
            out[0] = v if v >= 0 else 255
        except ValueError:
            pass
        return out
    a, b = gt_field.split(sep, 1)
    for i, tok in enumerate((a, b)):
        if tok == '.' or tok == '':
            out[i] = 255
        else:
            try:
                v = int(tok)
                out[i] = v if v >= 0 else 255
            except ValueError:
                out[i] = 255
    return out


def _call_char_ord(alleles: List[str], a0: int, a1: int) -> np.uint8:
    """
    Convert two allele indexes to a single-byte ASCII code:
      - 'N' if missing/invalid or not single-base A/C/G/T
      - IUPAC code for heterozygous SNPs
      - base letter for homozygous SNPs
    """
    N = np.uint8(ord('N'))

    # both alleles must be present
    if a0 == 255 or a1 == 255:
        return N
    # indexes in range?
    if a0 < 0 or a1 < 0 or a0 >= len(alleles) or a1 >= len(alleles):
        return N

    b0, b1 = alleles[a0], alleles[a1]

    # Only accept single-base A/C/G/T
    valid = {"A", "C", "G", "T"}
    if b0 not in valid or b1 not in valid:
        return N

    if b0 == b1:
        return np.uint8(ord(b0))
    # heterozygous: map to IUPAC (unordered)
    code = _IUPAC.get(frozenset((b0, b1)))
    return np.uint8(ord(code)) if code else N


def iter_vcf_filtered_snps_with_bed(
    vcf_path: Union[str, Path],
    bed_path: Union[str, Path],
    snames: List[str],
    scaffold_order: Optional[List[str]] = None,
) -> Iterator[Tuple[int, int, int, int, int, str, str, str, str, np.ndarray]]:
    """
    Iterate a VCF(.gz) and yield records with FILTER==PASS and NOT having 'INDEL' in INFO.
    For each record, also return genotypes for `snames` as a (nsamples, 3) uint8 array:
      [:,0:2] = allele indexes (0=REF, 1=first ALT, ...; 255=missing)
      [:,2]   = ASCII ord of a single-letter call:
                'N' for missing/undeterminable,
                IUPAC code for heterozygous SNPs,
                base letter for homozygous SNPs.

    Yields:
      (bed_idx, var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name, REF, ALT, QUAL, GT)

    Raises:
      - ValueError if a kept variant is not inside a BED interval
      - ValueError if any `snames` are absent from the VCF header
      - ValueError if a record lacks a GT field in FORMAT
    """
    bed_index, scaff2idx = load_bed_index_nonoverlap(bed_path, scaffold_order)

    # per-interval counters (for var_idx_in_bed)
    counters: Dict[str, np.ndarray] = {
        chrom: np.zeros(len(bed_index[chrom][0]), dtype=np.int64)
        for chrom in bed_index
    }

    vcf_path = str(vcf_path)
    _open = gzip.open if vcf_path.endswith(".gz") else open

    sample_cols: Optional[List[int]] = None

    with _open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                name_to_col = {nm: 9 + i for i, nm in enumerate(vcf_samples)}
                missing = [s for s in snames if s not in name_to_col]
                if missing:
                    raise ValueError(f"Samples not found in VCF header: {missing}")
                sample_cols = [name_to_col[s] for s in snames]
                continue

            if sample_cols is None:
                raise ValueError("VCF header (#CHROM) not found before records.")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, spos, _id, REF, ALT, QUAL, FILTER, INFO = parts[:8]

            # gates
            if FILTER != "PASS":
                continue
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
            if not (i >= 0 and pos0 < ends[i]):
                left_end = ends[i] if i >= 0 else None
                right_start = starts[i+1] if (i+1) < len(starts) else None
                raise ValueError(
                    f"Variant {chrom}:{pos} (pos0={pos0}) not covered by BED. "
                    f"Nearest interval ends at {left_end} and next starts at {right_start}."
                )

            if len(parts) < 9:
                raise ValueError(f"Record {chrom}:{pos} missing FORMAT/sample columns.")
            FORMAT = parts[8]
            fmt_keys = FORMAT.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                raise ValueError(f"Record {chrom}:{pos} lacks GT in FORMAT: {FORMAT}")

            # Build GT array (nsamples, 3)
            ns = len(snames)
            GT = np.empty((ns, 3), dtype=np.uint8)

            # Prepare allele string list: [REF, ALT1, ALT2, ...]
            alt_list = ALT.split(",") if ALT else []
            alleles_str = [REF] + alt_list

            for k, col in enumerate(sample_cols):
                # default missing
                a = np.array([255, 255], dtype=np.uint8)

                if col < len(parts):
                    sample_field = parts[col]
                    if sample_field and sample_field != ".":
                        sub = sample_field.split(":")
                        if gt_idx < len(sub):
                            a = _parse_gt_to_uint8(sub[gt_idx])

                GT[k, 0:2] = a

                # third column: ord of call char
                GT[k, 2] = _call_char_ord(alleles_str, int(a[0]), int(a[1]))

            var_idx_in_bed = int(counters[chrom][i])
            counters[chrom][i] += 1
            offset_in_bed = int(pos0 - starts[i])

            yield int(bedix[i]), var_idx_in_bed, offset_in_bed, scaff_idx, pos0, scaff_name, REF, ALT, QUAL, GT


if __name__ == "__main__":

    from ipyrad2.utils.logger import set_log_level
    set_log_level("DEBUG")
    REF = Path("/home/deren/Documents/ipyrad-tests/examples/Atub-genome/AmaTu_v01_no00_renamed.fa")
    VCF = Path("/home/deren/Documents/ipyrad-tests/Ama-out/assembly.vcf.gz")
    BED = Path("/home/deren/Documents/ipyrad-tests/Ama-out/assembly.bed")
    # bdict = load_bed_index_nonoverlap(BED)
    # print(bdict["A_tuberculatus_Chr01"])
    # ii = iter_vcf_with_bed(VCF, bdict, False)
    snames = [
        "SLH_AL_0072-contemp",
        "SLH_AL_0077-contemp",
        "SLH_AL_0078-contemp",
        "SLH_AL_0079-contemp",
        "SLH_AL_0080-contemp",
        "SLH_AL_0084-contemp",
        "SLH_AL_0086-contemp",
    ]

    # ii = iter_vcf_filtered_snps_with_bed(VCF, BED, snames)

    # for i in range(10):
    #     data = next(ii)
    #     print(data)

    write_snps_hdf5("TEST3", Path("/tmp"), snames, REF, VCF, BED, )
