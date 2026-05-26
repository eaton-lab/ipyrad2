#!/usr/bin/env python

"""Write the SNP datasets into the final assemble HDF5 output.

snpsmap
-------
description: The map of SNP positions on RAD loci and genome scaffolds.
Loci, scaffolds, and positions are all stored 0-indexed.
dtype: adaptive uint32/uint64
shape: (nsnps, 5)
attrs["columns"]: ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
attrs["indexing"]: [0, 0, 0, 0, 0]

genos
-----
description: The ordered diploid genotype call of the ordered samples at
every SNP position.
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

sample_dp
---------
description: The ordered per-sample FORMAT/DP value at every SNP position.
dtype: np.uint32
shape: (nsamples, nsnps)

site_qual
---------
description: The ordered VCF QUAL value at every SNP position.
dtype: np.float32
shape: (nsnps,)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Dict, List, Tuple, Union

import h5py
from loguru import logger
import numpy as np
import pandas as pd

from ..utils.parallel import run_pipeline
from ..utils.parallel import run_with_pool
from .hdf5_utils import choose_hdf5_cache_settings
from .hdf5_utils import choose_unsigned_int_dtype
from .hdf5_utils import format_bytes
from .hdf5_utils import get_fai_values


BIN = Path(sys.prefix) / "bin"
BIN_BCF = str(BIN / "bcftools")
_N_ORD = np.uint8(ord("N"))
_VALID_BASE_BYTES = np.array(
    [ord("A"), ord("C"), ord("G"), ord("T")],
    dtype=np.uint8,
)
_IUPAC_BYTE_LOOKUP = np.full((256, 256), _N_ORD, dtype=np.uint8)
for (base0, base1), code in {
    ("A", "G"): "R",
    ("C", "T"): "Y",
    ("G", "C"): "S",
    ("A", "T"): "W",
    ("G", "T"): "K",
    ("A", "C"): "M",
}.items():
    b0 = ord(base0)
    b1 = ord(base1)
    _IUPAC_BYTE_LOOKUP[b0, b1] = ord(code)
    _IUPAC_BYTE_LOOKUP[b1, b0] = ord(code)


@dataclass(frozen=True)
class _SnpChunkPlan:
    """Metadata for one chunk-local SNP query job."""

    chunk_idx: int
    chunk_bed: Path
    bed_offset: int


def _choose_chunk_snps(
    nsamples: int,
    target_mb: int = 16,
    min_snps: int = 4_096,
    max_snps: int = 262_144,
) -> int:
    """Choose the SNP-axis chunk size for the read-optimized HDF5 datasets."""
    bytes_per_snp = max(1, nsamples * 3)
    by_size = int((target_mb * 1024 * 1024) // bytes_per_snp)
    return max(min_snps, min(by_size, max_snps))


def _choose_worker_count(cores: int, threads: int) -> int:
    """Bound SNP-writer workers to the assemble job budget."""
    return max(1, int(cores) // max(1, int(threads)))


def _choose_chunk_count(
    nloci: int,
    worker_count: int,
    requested_chunks: int | None = None,
) -> int:
    """Choose a bounded number of loci chunks for SNP conversion."""
    if requested_chunks is not None:
        return max(1, min(int(requested_chunks), nloci))
    return max(1, min(nloci, max(8, worker_count * 8)))


def load_bed_index_nonoverlap(
    bed_path: Union[str, Path],
) -> Tuple[Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]], Dict[str, int]]:
    """Load a non-overlapping BED into per-chrom arrays and a scaffold index map."""
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        dtype={0: str, 1: "Int64", 2: "Int64"},
        engine="c",
        na_filter=False,
    )
    df = df[[0, 1, 2]].rename(columns={0: "chrom", 1: "start", 2: "end"})
    df["bed_idx"] = df.index.to_numpy()

    bed_index: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    chrom_order: List[str] = []

    for chrom, sub in df.groupby("chrom", sort=False):
        chrom_order.append(chrom)
        sub = sub.sort_values("start", kind="mergesort")
        starts = sub["start"].to_numpy(dtype=np.int64, copy=True)
        ends = sub["end"].to_numpy(dtype=np.int64, copy=True)
        idxs = sub["bed_idx"].to_numpy(dtype=np.int64, copy=True)
        if starts.size and np.any(starts[1:] < ends[:-1]):
            raise ValueError(
                f"BED intervals overlap on {chrom}, but non-overlap was assumed."
            )
        bed_index[chrom] = (starts, ends, idxs)

    scaff2idx = {name: i for i, name in enumerate(chrom_order)}
    return bed_index, scaff2idx


def _write_chunk_beds_with_offsets(
    *,
    loci_bed: Path,
    workdir: Path,
    nchunks: int,
) -> list[_SnpChunkPlan]:
    """Split one BED into chunk BEDs while tracking global BED row offsets."""
    lines = [line for line in loci_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"No loci found in {loci_bed}.")

    nchunks = max(1, min(int(nchunks), len(lines)))
    q, r = divmod(len(lines), nchunks)
    plans: list[_SnpChunkPlan] = []
    offset = 0
    cursor = 0
    for chunk_idx in range(nchunks):
        size = q + (1 if chunk_idx < r else 0)
        chunk_lines = lines[cursor : cursor + size]
        chunk_bed = workdir / f"snp-writer-chunk-{chunk_idx}.bed"
        chunk_bed.write_text("\n".join(chunk_lines) + "\n", encoding="utf-8")
        plans.append(
            _SnpChunkPlan(
                chunk_idx=chunk_idx,
                chunk_bed=chunk_bed,
                bed_offset=offset,
            )
        )
        offset += size
        cursor += size
    return plans


def _return_gt_allele(gt_field: str, index: int) -> np.uint8:
    """Return one allele index from a bare GT token, else 255 for missing."""
    try:
        value = gt_field[index]
    except (IndexError, TypeError):
        return np.uint8(255)
    if value == ".":
        return np.uint8(255)
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return np.uint8(255)
    return np.uint8(parsed if parsed >= 0 else 255)


_v_return_gt_allele = np.vectorize(_return_gt_allele, otypes=[np.uint8])


def _return_dp_value(sample_field: str) -> np.uint32:
    """Return one non-negative DP value from a `GT:DP` query token, else 0."""
    if not isinstance(sample_field, str):
        return np.uint32(0)
    parts = sample_field.split(":", 1)
    if len(parts) != 2:
        return np.uint32(0)
    dp_field = parts[1].strip()
    if dp_field in {"", "."}:
        return np.uint32(0)
    try:
        parsed = int(dp_field)
    except ValueError:
        return np.uint32(0)
    return np.uint32(max(0, parsed))


_v_return_dp_value = np.vectorize(_return_dp_value, otypes=[np.uint32])


def _safe_parse_qual(value) -> np.float32:
    """Return one float QUAL value or NaN when unavailable."""
    text = str(value).strip()
    if text in {"", "."}:
        return np.float32(np.nan)
    try:
        return np.float32(float(text))
    except ValueError:
        return np.float32(np.nan)


def _chunk_query_to_gt_arrays(
    chunkdf: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Convert one bcftools GT-only query chunk into genotype and REF arrays."""
    nsnps = int(chunkdf.shape[0])
    nsamples = max(0, int(chunkdf.shape[1]) - 5)

    ref = np.frombuffer(
        "".join(chunkdf.iloc[:, 3].astype(str)).encode("ascii"),
        dtype=np.uint8,
    )
    alts = chunkdf.iloc[:, 4].astype(bytes).to_numpy(copy=False)
    alt_bytes = np.char.replace(alts, b",", b"")
    alts1 = np.zeros(nsps := nsnps, dtype=np.uint8)
    alts2 = np.zeros(nsps, dtype=np.uint8)
    alts3 = np.zeros(nsps, dtype=np.uint8)
    lengths = np.fromiter((len(value) for value in alt_bytes), count=nsps, dtype=np.int64)
    if np.any(lengths >= 1):
        alts1[lengths >= 1] = [value[0] for value in alt_bytes[lengths >= 1]]
    if np.any(lengths >= 2):
        alts2[lengths >= 2] = [value[1] for value in alt_bytes[lengths >= 2]]
    if np.any(lengths >= 3):
        alts3[lengths >= 3] = [value[2] for value in alt_bytes[lengths >= 3]]

    gt_fields = chunkdf.iloc[:, 5:]
    g0 = _v_return_gt_allele(gt_fields, 0)
    g1 = _v_return_gt_allele(gt_fields, 2)
    sample_dp = _v_return_dp_value(gt_fields)
    site_qual = np.array(
        [_safe_parse_qual(value) for value in chunkdf.iloc[:, 2]],
        dtype=np.float32,
    )

    alleles = np.column_stack((ref, alts1, alts2, alts3)).astype(np.uint8, copy=False)
    row_index = np.arange(nsnps, dtype=np.int64)[:, None]
    a0 = alleles[row_index, np.minimum(g0, 3)]
    a1 = alleles[row_index, np.minimum(g1, 3)]

    snps = np.full((nsnps, nsamples), _N_ORD, dtype=np.uint8)
    homo = (g0 == g1) & (g0 <= 3)
    homo_bases = a0
    homo_valid = homo & np.isin(homo_bases, _VALID_BASE_BYTES)
    snps[homo_valid] = homo_bases[homo_valid]

    het = (g0 != g1) & (g0 <= 3) & (g1 <= 3)
    het_codes = _IUPAC_BYTE_LOOKUP[a0, a1]
    het_valid = het & (het_codes != _N_ORD)
    snps[het_valid] = het_codes[het_valid]

    genos = np.empty((nsamples, nsnps, 3), dtype=np.uint8)
    genos[:, :, 0] = g0.T
    genos[:, :, 1] = g1.T
    genos[:, :, 2] = snps.T
    return genos, ref, snps, sample_dp.T, site_qual


def _build_chunk_snpsmap(
    chunkdf: pd.DataFrame,
    *,
    bed_index: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    scaff2idx: Dict[str, int],
    bed_offset: int,
    counters: Dict[str, np.ndarray],
    map_dtype: np.dtype,
) -> np.ndarray:
    """Build `snpsmap` rows for one GT query chunk."""
    chunkdf = chunkdf.reset_index(drop=True)
    snpsmap = np.empty((chunkdf.shape[0], 5), dtype=map_dtype)

    for chrom, sub in chunkdf.groupby(0, sort=False):
        rows = sub.index.to_numpy(dtype=np.int64, copy=False)
        pos0 = sub.iloc[:, 1].to_numpy(dtype=np.int64, copy=False) - 1
        tup = bed_index.get(chrom)
        if tup is None:
            raise ValueError(f"Variant chunk contains chromosome absent from BED: {chrom}")
        starts, ends, bedix = tup
        interval_idx = np.searchsorted(starts, pos0, side="right") - 1
        valid = (interval_idx >= 0) & (pos0 < ends[interval_idx])
        if not np.all(valid):
            bad_pos0 = int(pos0[np.flatnonzero(~valid)[0]])
            raise ValueError(
                f"Variant {chrom}:{bad_pos0 + 1} (pos0={bad_pos0}) not covered by chunk BED."
            )

        snpsmap[rows, 0] = bed_offset + bedix[interval_idx]
        snpsmap[rows, 2] = pos0 - starts[interval_idx]
        snpsmap[rows, 3] = scaff2idx[chrom]
        snpsmap[rows, 4] = pos0

        # Query output is sorted by coordinate, so variants inside the same BED
        # interval arrive contiguously and can be counted in stable runs.
        breaks = np.flatnonzero(np.diff(interval_idx)) + 1
        for run in np.split(np.arange(interval_idx.size), breaks):
            local_idx = int(interval_idx[run[0]])
            start = int(counters[chrom][local_idx])
            count = int(run.size)
            snpsmap[rows[run], 1] = np.arange(start, start + count, dtype=map_dtype)
            counters[chrom][local_idx] += count

    return snpsmap


def _write_snp_chunk_worker(
    *,
    chunk_idx: int,
    chunk_bed: Path,
    bed_offset: int,
    vcf_path: Path,
    snames: list[str],
    workdir: Path,
    scaff2idx: Dict[str, int],
    map_dtype: np.dtype,
    vcf_chunk_rows: int,
) -> Dict[str, object]:
    """Materialize one BED-local SNP chunk to temporary `.npy` arrays."""
    workdir.mkdir(parents=True, exist_ok=True)
    query_tsv = workdir / f"snp-query-{chunk_idx}.tsv"
    out_map = workdir / f"snp-chunk-{chunk_idx}.snpsmap.npy"
    out_gen = workdir / f"snp-chunk-{chunk_idx}.genos.npy"
    out_ref = workdir / f"snp-chunk-{chunk_idx}.reference.npy"
    out_dp = workdir / f"snp-chunk-{chunk_idx}.sample_dp.npy"
    out_qual = workdir / f"snp-chunk-{chunk_idx}.site_qual.npy"

    try:
        sample_arg = ",".join(snames)
        query_cmd = [
            BIN_BCF,
            "query",
            "-R",
            str(chunk_bed),
            "-s",
            sample_arg,
            "-i",
            'FILTER="PASS" && TYPE="snp"',
            "-f",
            "%CHROM\t%POS\t%QUAL\t%REF\t%ALT[\t%GT:%DP]\n",
            str(vcf_path),
        ]
        run_pipeline([query_cmd], outfile=query_tsv)

        if not query_tsv.exists() or query_tsv.stat().st_size == 0:
            np.save(out_map, np.empty((0, 5), dtype=map_dtype), allow_pickle=False)
            np.save(
                out_gen,
                np.empty((len(snames), 0, 3), dtype=np.uint8),
                allow_pickle=False,
            )
            np.save(out_ref, np.empty((0,), dtype=np.uint8), allow_pickle=False)
            np.save(
                out_dp,
                np.empty((len(snames), 0), dtype=np.uint32),
                allow_pickle=False,
            )
            np.save(out_qual, np.empty((0,), dtype=np.float32), allow_pickle=False)
            return {
                "chunk_idx": chunk_idx,
                "nsnps": 0,
                "snpsmap": str(out_map),
                "genos": str(out_gen),
                "reference": str(out_ref),
                "sample_dp": str(out_dp),
                "site_qual": str(out_qual),
            }

        bed_index, _local_scaff2idx = load_bed_index_nonoverlap(chunk_bed)
        counters: Dict[str, np.ndarray] = {
            chrom: np.zeros(len(starts), dtype=np.int64)
            for chrom, (starts, _ends, _idxs) in bed_index.items()
        }
        maps: list[np.ndarray] = []
        genos: list[np.ndarray] = []
        refs: list[np.ndarray] = []
        depths: list[np.ndarray] = []
        quals: list[np.ndarray] = []
        total = 0

        for chunkdf in pd.read_csv(
            query_tsv,
            sep="\t",
            header=None,
            dtype=str,
            na_filter=False,
            engine="c",
            chunksize=int(vcf_chunk_rows),
        ):
            chunkdf = chunkdf.reset_index(drop=True)
            chunk_map = _build_chunk_snpsmap(
                chunkdf,
                bed_index=bed_index,
                scaff2idx=scaff2idx,
                bed_offset=bed_offset,
                counters=counters,
                map_dtype=map_dtype,
            )
            chunk_genos, chunk_ref, _chunk_snps, chunk_dp, chunk_qual = _chunk_query_to_gt_arrays(
                chunkdf
            )
            maps.append(chunk_map)
            genos.append(chunk_genos)
            refs.append(chunk_ref)
            depths.append(chunk_dp)
            quals.append(chunk_qual)
            total += int(chunkdf.shape[0])

        map_arr = (
            np.concatenate(maps, axis=0)
            if maps
            else np.empty((0, 5), dtype=map_dtype)
        )
        gen_arr = (
            np.concatenate(genos, axis=1)
            if genos
            else np.empty((len(snames), 0, 3), dtype=np.uint8)
        )
        ref_arr = (
            np.concatenate(refs, axis=0)
            if refs
            else np.empty((0,), dtype=np.uint8)
        )
        dp_arr = (
            np.concatenate(depths, axis=1)
            if depths
            else np.empty((len(snames), 0), dtype=np.uint32)
        )
        qual_arr = (
            np.concatenate(quals, axis=0)
            if quals
            else np.empty((0,), dtype=np.float32)
        )
        np.save(out_map, map_arr, allow_pickle=False)
        np.save(out_gen, gen_arr, allow_pickle=False)
        np.save(out_ref, ref_arr, allow_pickle=False)
        np.save(out_dp, dp_arr, allow_pickle=False)
        np.save(out_qual, qual_arr, allow_pickle=False)
        return {
            "chunk_idx": chunk_idx,
            "nsnps": int(total),
            "snpsmap": str(out_map),
            "genos": str(out_gen),
            "reference": str(out_ref),
            "sample_dp": str(out_dp),
            "site_qual": str(out_qual),
        }
    finally:
        query_tsv.unlink(missing_ok=True)


def _append_snp_chunk(
    *,
    result: Dict[str, object],
    snpsmap_ds: h5py.Dataset,
    genos_ds: h5py.Dataset,
    reference_ds: h5py.Dataset,
    sample_dp_ds: h5py.Dataset,
    site_qual_ds: h5py.Dataset,
    total: int,
) -> int:
    """Append one temporary chunk's arrays into the final HDF5 datasets."""
    map_path = Path(str(result["snpsmap"]))
    gen_path = Path(str(result["genos"]))
    ref_path = Path(str(result["reference"]))
    dp_path = Path(str(result["sample_dp"]))
    qual_path = Path(str(result["site_qual"]))
    try:
        map_arr = np.load(map_path, allow_pickle=False)
        if map_arr.shape[0] == 0:
            return total
        gen_arr = np.load(gen_path, allow_pickle=False)
        ref_arr = np.load(ref_path, allow_pickle=False)
        dp_arr = np.load(dp_path, allow_pickle=False)
        qual_arr = np.load(qual_path, allow_pickle=False)
        new_total = total + int(map_arr.shape[0])
        nsamples = int(genos_ds.shape[0])
        snpsmap_ds.resize((new_total, 5))
        genos_ds.resize((nsamples, new_total, 3))
        reference_ds.resize((new_total,))
        sample_dp_ds.resize((nsamples, new_total))
        site_qual_ds.resize((new_total,))
        snpsmap_ds[total:new_total, :] = map_arr
        genos_ds[:, total:new_total, :] = gen_arr
        reference_ds[total:new_total] = ref_arr
        sample_dp_ds[:, total:new_total] = dp_arr
        site_qual_ds[total:new_total] = qual_arr
        return new_total
    finally:
        map_path.unlink(missing_ok=True)
        gen_path.unlink(missing_ok=True)
        ref_path.unlink(missing_ok=True)
        dp_path.unlink(missing_ok=True)
        qual_path.unlink(missing_ok=True)


def write_snps_hdf5(
    name: str,
    outdir: Path,
    snames: List[str],
    reference: Path,
    *,
    tmpdir: Path | None = None,
    cores: int = 1,
    threads: int = 1,
    log_level: str = "INFO",
    chunk_count: int | None = None,
    vcf_chunk_rows: int = 100_000,
) -> int:
    """Write the final assemble SNP HDF5 by chunked GT-only VCF conversion."""
    database = outdir / f"{name}.hdf5"
    vcf_path = outdir / f"{name}.vcf.gz"
    loci_bed = outdir / f"{name}.bed"
    workdir = Path(tmpdir) if tmpdir is not None else outdir
    snp_tmpdir = workdir / "snp_writer_tmp"
    snp_tmpdir.mkdir(parents=True, exist_ok=True)

    snames = sorted(snames)
    nsamples = len(snames)
    chunk_snps = _choose_chunk_snps(nsamples, target_mb=128, max_snps=131_072)
    bed_index, scaff2idx = load_bed_index_nonoverlap(loci_bed)
    nloci = sum(len(starts) for starts, _ends, _idxs in bed_index.values())
    max_locus_length = max(
        (
            int((ends - starts).max())
            for starts, ends, _idxs in bed_index.values()
            if len(starts)
        ),
        default=0,
    )
    map_dtype = choose_unsigned_int_dtype(
        max(
            nloci,
            max_locus_length,
            len(scaff2idx) - 1,
            max((int(i) for i in get_fai_values(reference, "length")), default=0),
        )
    )
    worker_count = _choose_worker_count(cores, threads)
    nchunks = _choose_chunk_count(nloci, worker_count, requested_chunks=chunk_count)
    plans = _write_chunk_beds_with_offsets(
        loci_bed=loci_bed,
        workdir=snp_tmpdir,
        nchunks=nchunks,
    )

    kwargs = choose_hdf5_cache_settings()
    string_dtype = h5py.string_dtype(encoding="utf-8")
    logger.debug(
        "snps writer config: chunk_snps={}, cache={}, cache_slots={}, map_dtype={}, nloci={}, workers={}, chunks={}",
        chunk_snps,
        format_bytes(int(kwargs["rdcc_nbytes"])),
        kwargs["rdcc_nslots"],
        map_dtype.name,
        nloci,
        worker_count,
        len(plans),
    )

    jobs = {
        plan.chunk_idx: (
            _write_snp_chunk_worker,
            dict(
                chunk_idx=plan.chunk_idx,
                chunk_bed=plan.chunk_bed,
                bed_offset=plan.bed_offset,
                vcf_path=vcf_path,
                snames=snames,
                workdir=snp_tmpdir,
                scaff2idx=scaff2idx,
                map_dtype=map_dtype,
                vcf_chunk_rows=vcf_chunk_rows,
            ),
        )
        for plan in plans
    }
    chunk_results = run_with_pool(
        jobs,
        log_level,
        max_workers=worker_count,
        msg="Building SNP database chunks",
    )

    try:
        with h5py.File(database, "a", **kwargs) as io5:
            snpsmap_ds = io5.create_dataset(
                "snpsmap",
                shape=(0, 5),
                maxshape=(None, 5),
                dtype=map_dtype,
                chunks=(chunk_snps, 5),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            snpsmap_ds.attrs["columns"] = ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
            snpsmap_ds.attrs["indexing"] = [0, 0, 0, 0, 0]

            genos_ds = io5.create_dataset(
                "genos",
                shape=(nsamples, 0, 3),
                maxshape=(nsamples, None, 3),
                dtype=np.uint8,
                chunks=(nsamples, chunk_snps, 3),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            genos_ds.attrs["names"] = np.array(snames, dtype=string_dtype)

            reference_ds = io5.create_dataset(
                "reference",
                shape=(0,),
                maxshape=(None,),
                dtype=np.uint8,
                chunks=(chunk_snps,),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            sample_dp_ds = io5.create_dataset(
                "sample_dp",
                shape=(nsamples, 0),
                maxshape=(nsamples, None),
                dtype=np.uint32,
                chunks=(nsamples, chunk_snps),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            site_qual_ds = io5.create_dataset(
                "site_qual",
                shape=(0,),
                maxshape=(None,),
                dtype=np.float32,
                chunks=(chunk_snps,),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )

            total = 0
            for chunk_idx in sorted(chunk_results):
                total = _append_snp_chunk(
                    result=chunk_results[chunk_idx],
                    snpsmap_ds=snpsmap_ds,
                    genos_ds=genos_ds,
                    reference_ds=reference_ds,
                    sample_dp_ds=sample_dp_ds,
                    site_qual_ds=site_qual_ds,
                    total=total,
                )
            io5.attrs["nsnps"] = int(total)
    finally:
        for plan in plans:
            plan.chunk_bed.unlink(missing_ok=True)

    if total == 0:
        logger.debug("wrote empty SNP dataset to {}", database)
    else:
        logger.debug(
            "wrote snps dataset to {} (nsnps={:,}, chunks={})",
            database,
            total,
            len(plans),
        )
    return int(total)
