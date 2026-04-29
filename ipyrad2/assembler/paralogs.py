#!/usr/bin/env python

"""Filter for paralogs both within and across samples.

Within each sample we:
    1. Call variants in mpileup/call
    2. Get beds of indel-affected regions
    3. Extract table of snps not affected by indels
    4. Analyze table to mark potential paralogs within samples
    5. Analyze tables across samples to mark potential paralogs across samples


# Metrics to compute per locus (region)
- site_count
- mean_dp, max_dp
- ax_abs_dp_z
- n_sites_dp_z_gt_T
- n_sites_3allele (third allele supported)
- max_third_frac
- n_sites_maf_gt_m (e.g. maf ≥ 0.2)
- scl_frac
"""

###########

import sys
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

import numpy as np
import pandas as pd

from ipyrad2.assembler.beds import sort_bed_by_reference_order
from ipyrad2.assembler.sort_utils import assemble_sort_with_args
from ipyrad2.utils.parallel import run_pipeline

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")


# -------------------------
# Variant calling + masking
# -------------------------


def call_vcf_from_bam(
    bam: Path,
    regions_bed: Path,
    reference_fasta: Path,
    out_vcf_gz: Path,
    *,
    min_map_q: int,
    min_base_q: int,
    max_depth: int = 250,
    threads_mpileup: int = 1,
    threads_call: int = 1,
    threads_norm: int = 1,
) -> Path:
    """Call SNPs/indels in regions_bed for one sample BAM; writes bgzipped VCF + index."""
    mpileup = [
        BIN_BCF,
        "mpileup",
        "-f",
        str(reference_fasta),
        "-q",
        str(min_map_q),
        "-Q",
        str(min_base_q),
        "-d",
        str(max_depth),
        "-a",
        "FMT/DP,FMT/AD",
        "-R",
        str(regions_bed),
        "--threads",
        str(threads_mpileup),
        "-Ob",
        str(bam),
    ]
    call = [
        BIN_BCF,
        "call",
        "-m",
        "-v",
        "-a",
        "GQ",
        "--ploidy",
        "2",
        "--threads",
        str(threads_call),
        "-Ob",
    ]
    norm = [
        BIN_BCF,
        "norm",
        "-f",
        str(reference_fasta),
        "-Oz",
        "-o",
        str(out_vcf_gz),
        "--threads",
        str(threads_norm),
        "--write-index",
    ]
    run_pipeline([mpileup, call, norm])
    return out_vcf_gz


def make_indel_mask_bed(
    vcf_gz: Path,
    out_bed: Path,
    *,
    pad_bp: int = 10,
) -> Path:
    """Write BED of indel intervals padded by +/- pad_bp."""
    view = [BIN_BCF, "view", "-v", "indels", str(vcf_gz)]
    query = [BIN_BCF, "query", "-f", "%CHROM\t%POS0\t%END\n"]
    pad = [
        "awk",
        "-v",
        f"X={pad_bp}",
        'BEGIN{OFS="\\t"}{s=$2-X; if(s<0)s=0; e=$3+X; print $1,s,e}',
    ]
    sort_ = assemble_sort_with_args(["-k1,1", "-k2,2n"])
    run_pipeline([view, query, pad, sort_], out_bed)
    return out_bed


def extract_snps_table_tsv(vcf_gz: Path, out_tsv: Path) -> Path:
    """Write TSV of SNPs: chrom, pos0, pos, DP, GQ, AD, GT (sorted)."""
    view = [BIN_BCF, "view", "-v", "snps", str(vcf_gz)]
    query = [
        BIN_BCF,
        "query",
        "-f",
        "%CHROM\t%POS0\t%POS\t[%DP]\t[%GQ]\t[%AD]\t[%GT]\n",
    ]
    sort_ = assemble_sort_with_args(["-k1,1", "-k2,2n"])
    run_pipeline([view, query, sort_], out_tsv)
    return out_tsv


def mask_snps_table_with_bed(snps_tsv: Path, mask_bed: Path, out_tsv: Path) -> Path:
    """Remove SNP rows that overlap mask_bed (bedtools intersect -v)."""
    cmd = [BIN_BED, "intersect", "-v", "-a", str(snps_tsv), "-b", str(mask_bed)]
    run_pipeline([cmd], out_tsv)
    return out_tsv


# -------------------------
# Softclipping by region
# -------------------------


def build_softclip_bam(bam: Path, out_bam: Path, softclip_len_threshold: int) -> None:
    """Write BAM containing reads with sclen > softclip_len_threshold (header preserved)."""
    cmd = [BIN_SAM, "view", "-hb", "-e", f"sclen > {softclip_len_threshold}", str(bam)]
    run_pipeline([cmd], out_bam)


def bedtools_coverage_counts(a_bed: Path, b_bam: Path, out_tsv: Path) -> None:
    """Run: bedtools coverage -a a_bed -b b_bam -counts; writes TSV."""
    cmd = [BIN_BED, "coverage", "-a", str(a_bed), "-b", str(b_bam), "-counts"]
    run_pipeline([cmd], out_tsv)


def read_bedtools_coverage_counts(path: Path, count_col: str) -> pd.DataFrame:
    """Read bedtools coverage -counts output for a 3-col BED; keep chrom,start,end,count."""
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: "string"})
    last = df.shape[1] - 1
    df = df[[0, 1, 2, last]].copy()
    df.columns = ["chrom", "start", "end", count_col]
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    df[count_col] = df[count_col].astype("int64")
    return df


def all_reads_by_region(
    bam: Path,
    regions_bed: Path,
    out_tsv: Path,
) -> pd.DataFrame:
    """Compute per-region read counts (all reads) using bedtools coverage -counts."""
    with NamedTemporaryFile(suffix=".all.cov.tsv", delete=False) as tmp:
        all_cov = Path(tmp.name)

    try:
        bedtools_coverage_counts(regions_bed, bam, all_cov)
        all_df = read_bedtools_coverage_counts(all_cov, "all_reads")
        all_df["rid"] = (
            all_df["chrom"].astype(str)
            + ":"
            + all_df["start"].astype(str)
            + "-"
            + all_df["end"].astype(str)
        )
        all_df.to_csv(out_tsv, sep="\t", index=False)
        return all_df
    finally:
        try:
            all_cov.unlink(missing_ok=True)
        except Exception:
            pass


def softclip_fraction_by_region(
    bam: Path,
    regions_bed: Path,
    out_tsv: Path,
    *,
    softclip_len_threshold: int,
) -> pd.DataFrame:
    """
    Compute per-region fraction of reads with sclen > softclip_len_threshold.

    Returns DataFrame with columns:
      chrom, start, end, all_reads, scl_reads, scl_frac, rid
    """
    with NamedTemporaryFile(suffix=".softclip.bam", delete=False) as tmp:
        clipped_bam = Path(tmp.name)

    with NamedTemporaryFile(suffix=".all.cov.tsv", delete=False) as tmp:
        all_cov = Path(tmp.name)

    with NamedTemporaryFile(suffix=".scl.cov.tsv", delete=False) as tmp:
        scl_cov = Path(tmp.name)

    try:
        build_softclip_bam(bam, clipped_bam, softclip_len_threshold)
        bedtools_coverage_counts(regions_bed, bam, all_cov)
        bedtools_coverage_counts(regions_bed, clipped_bam, scl_cov)

        all_df = read_bedtools_coverage_counts(all_cov, "all_reads")
        scl_df = read_bedtools_coverage_counts(scl_cov, "scl_reads")

        merged = all_df.merge(scl_df, on=["chrom", "start", "end"], how="left")
        merged["scl_reads"] = merged["scl_reads"].fillna(0).astype("int64")

        denom = merged["all_reads"].to_numpy(dtype=float)
        numer = merged["scl_reads"].to_numpy(dtype=float)
        merged["scl_frac"] = np.divide(
            numer, denom, out=np.zeros_like(numer), where=denom > 0
        )

        merged["rid"] = (
            merged["chrom"].astype(str)
            + ":"
            + merged["start"].astype(str)
            + "-"
            + merged["end"].astype(str)
        )

        merged.to_csv(out_tsv, sep="\t", index=False)
        return merged

    finally:
        try:
            clipped_bam.unlink(missing_ok=True)
            all_cov.unlink(missing_ok=True)
            scl_cov.unlink(missing_ok=True)
        except Exception:
            pass


# -------------------------
# SNP table -> locus metrics
# -------------------------


def read_snps_table(path: str | Path) -> pd.DataFrame:
    """Read SNP TSV: chrom, start, end, DP, GQ, AD, GT."""
    dtypes = {
        "chrom": "string",
        "start": "int64",
        "end": "int64",
        "DP": "int64",
        "GQ": "int64",
        "AD": "string",
        "GT": "string",
    }
    path = Path(path)
    if (not path.exists()) or path.stat().st_size == 0:
        return pd.DataFrame(
            {key: pd.Series(dtype=value) for key, value in dtypes.items()}
        )
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "DP", "GQ", "AD", "GT"],
            dtype=dtypes,
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(
            {key: pd.Series(dtype=value) for key, value in dtypes.items()}
        )
    return df.sort_values(["chrom", "start"]).reset_index(drop=True)


def read_regions_bed(path: str | Path) -> pd.DataFrame:
    """Read 3-col BED, preserving canonical input order, and add rid keys."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": "string", "start": "int64", "end": "int64"},
    )
    df = df.reset_index(drop=True)
    df["rid"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + "-"
        + df["end"].astype(str)
    )
    return df


def compute_dp_z(dp: np.ndarray) -> np.ndarray:
    """Compute z-scores of log1p(DP) across all rows."""
    x = np.log1p(dp.astype(float))
    mu = float(x.mean()) if x.size else 0.0
    sd = float(x.std(ddof=1)) if x.size > 1 else 1.0
    sd = max(sd, 1e-12)
    return (x - mu) / sd


def allele_depth_metrics(
    ad_str: str, dp: int, min_allele_depth: int
) -> tuple[int, float, float]:
    """Return (alleles_ge_min_ad, maf, third_frac) from AD and DP."""
    if dp <= 0 or ad_str is None:
        return 0, 0.0, 0.0

    ad = np.fromstring(ad_str, sep=",", dtype=int)
    if ad.size == 0:
        return 0, 0.0, 0.0

    alleles_ge_min_ad = int((ad >= min_allele_depth).sum())

    ad_sorted = np.sort(ad)[::-1]
    a2 = int(ad_sorted[1]) if ad_sorted.size >= 2 else 0
    a3 = int(ad_sorted[2]) if ad_sorted.size >= 3 else 0

    maf = a2 / dp
    third_frac = a3 / dp
    return alleles_ge_min_ad, maf, third_frac


def add_site_metrics(df_snps: pd.DataFrame, min_allele_depth: int) -> pd.DataFrame:
    """Add dp_z, alleles_ge_min_ad, maf, third_frac."""
    out = df_snps.copy()
    out["dp_z"] = compute_dp_z(out["DP"].to_numpy())

    # `DataFrame.apply(..., result_type="expand")` returns an empty frame that
    # mirrors the input columns when there are no rows, so handle the no-SNP
    # case explicitly before expanding per-row allele metrics.
    if out.empty:
        out["alleles_ge_min_ad"] = pd.Series(dtype="int64")
        out["maf"] = pd.Series(dtype="float64")
        out["third_frac"] = pd.Series(dtype="float64")
        return out

    vals = out.apply(
        lambda r: allele_depth_metrics(str(r["AD"]), int(r["DP"]), min_allele_depth),
        axis=1,
        result_type="expand",
    )
    vals.columns = ["alleles_ge_min_ad", "maf", "third_frac"]
    return pd.concat([out, vals], axis=1)


def assign_regions_nonoverlap(
    df_sites: pd.DataFrame, df_regions: pd.DataFrame
) -> pd.DataFrame:
    """Assign each SNP row to a non-overlapping region rid; drops rows outside any region."""
    sites = df_sites.copy()
    sites["rid"] = pd.NA

    for chrom, sidx in sites.groupby("chrom", sort=False).groups.items():
        r = df_regions[df_regions["chrom"] == chrom]
        if r.empty:
            continue

        r_start = r["start"].to_numpy()
        r_end = r["end"].to_numpy()
        r_id = r["rid"].to_numpy()

        pos = sites.loc[sidx, "start"].to_numpy()
        j = np.searchsorted(r_start, pos, side="right") - 1
        j_clip = np.clip(j, 0, len(r_end) - 1)
        ok = (j >= 0) & (pos < r_end[j_clip])

        rid = np.full(pos.shape, None, dtype=object)
        rid[ok] = r_id[j[ok]]
        sites.loc[sidx, "rid"] = rid

    return sites.dropna(subset=["rid"]).reset_index(drop=True)


def summarize_loci(
    df_sites: pd.DataFrame,
    *,
    max_abs_dp_z_max: float,
    third_frac_cut: float,
    min_3allele_sites: int,
    maf_threshold: float,
    max_sites_above_maf: int,
) -> pd.DataFrame:
    """
    Aggregate site metrics per region.

    Parameters
    ----------
    max_abs_dp_z_max:
        Locus fails if max(abs(dp_z)) exceeds this.
    third_frac_cut:
        Site has strong 3rd-allele evidence if third_frac >= third_frac_cut.
    min_3allele_sites:
        Locus fails if it has >= this many sites with >=3 alleles supported.
    maf_threshold:
        Site is heterozygous with maf above threshold.
    n_sites_above_maf
        Number of heterozygous sites, used as a cutoff for outlier amounts
    """
    df = df_sites.copy()
    df["is_3allele_site"] = df["alleles_ge_min_ad"] >= 3
    df["is_third_strong"] = df["third_frac"] >= third_frac_cut
    df["is_high_maf"] = df["maf"] >= maf_threshold

    locus = (
        df.groupby("rid", as_index=False)
        .agg(
            n_snps=("rid", "size"),
            max_snp_site_depth=("DP", "max"),
            max_snp_depth_z=(
                "dp_z",
                lambda x: float(np.abs(np.asarray(x)).max()) if len(x) else 0.0,
            ),
            n_sites_3allele=("is_3allele_site", "sum"),
            n_sites_3allele_strong=("is_third_strong", "sum"),
            max_frac_3allele=("third_frac", "max"),
            n_sites_above_maf=("is_high_maf", "sum"),
            max_maf=("maf", "max"),
        )
        .sort_values("rid")
        .reset_index(drop=True)
    )

    # Individual failure flags are kept in the table so the caller can see
    # which signal marked a region as paralog-like.
    locus["fail_dp_z"] = locus["max_snp_depth_z"] > max_abs_dp_z_max
    locus["fail_3allele"] = locus["n_sites_3allele_strong"] >= min_3allele_sites
    locus["fail_maf"] = locus["n_sites_above_maf"] > max_sites_above_maf
    return locus


def add_softclip_and_pass(
    df_locus: pd.DataFrame,
    *,
    df_softclip: pd.DataFrame | None,
    softclip_frac_max: float | None,
) -> pd.DataFrame:
    """
    Merge softclip metrics into locus table and add pass column.

    Parameters
    ----------
    df_softclip:
        Output of softclip_fraction_by_region() with columns rid, all_reads, scl_reads, scl_frac.
    softclip_frac_max:
        If set, locus fails if scl_frac > softclip_frac_max.
    """
    out = df_locus.copy()

    if df_softclip is not None:
        out = out.merge(
            df_softclip[["rid", "all_reads", "scl_reads", "scl_frac"]],
            on="rid",
            how="left",
        )
        out = out.fillna({"all_reads": 0, "scl_reads": 0, "scl_frac": 0.0})
    else:
        out["all_reads"] = 0
        out["scl_reads"] = 0
        out["scl_frac"] = 0.0

    if softclip_frac_max is None:
        out["fail_softclip"] = False
    else:
        out["fail_softclip"] = out["scl_frac"] > softclip_frac_max

    fail_cols = [c for c in out.columns if c.startswith("fail_")]
    out["pass"] = ~out[fail_cols].any(axis=1)
    out["paralog_like"] = ~out["pass"]
    return out


def write_outputs(
    df_sites: pd.DataFrame,
    df_locus: pd.DataFrame,
    df_regions: pd.DataFrame,
    prefix: str,
) -> None:
    """Write site/locus TSVs and good/paralog_like 3-col BEDs."""
    prefix = str(prefix)

    df_sites.to_csv(f"{prefix}.site_metrics.tsv", sep="\t", index=False)
    df_locus.to_csv(f"{prefix}.locus_metrics.tsv", sep="\t", index=False)

    bad = set(df_locus.loc[df_locus["paralog_like"], "rid"])
    bad_bed = df_regions[df_regions["rid"].isin(bad)][["chrom", "start", "end"]]
    good_bed = df_regions[~df_regions["rid"].isin(bad)][["chrom", "start", "end"]]

    bad_bed.to_csv(f"{prefix}.paralog_like.bed", sep="\t", header=False, index=False)
    good_bed.to_csv(f"{prefix}.good.bed", sep="\t", header=False, index=False)


# -------------------------
# Orchestration
# -------------------------


def get_sample_paralog_tables(
    bam: Path,
    regions_bed: Path,
    reference_fasta: Path,
    tmpdir: Path,
    prefix: str,
    min_map_q: int,
    min_base_q: int,
    indel_pad_bp: int = 10,
    min_allele_depth: int = 2,
    max_abs_dp_z_max: float = 5.0,
    third_frac_cut: float = 0.10,
    min_3allele_sites: int = 2,
    maf_threshold: float = 0.20,
    max_sites_above_maf: int = 8,
    softclip_len_threshold: int | None = None,
    softclip_frac_max: float | None = None,
    callable_regions_bed: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Within-sample workflow returning (site_table, locus_table).

    If softclip_len_threshold is not None, computes scl_frac by region and merges.
    If softclip_frac_max is not None, marks loci failing scl_frac > softclip_frac_max.
    The returned locus table includes every region in regions_bed, even when a
    locus had data but no retained SNPs after masking.
    """
    tmpdir.mkdir(parents=True, exist_ok=True)

    vcf = tmpdir / f"{prefix}.vcf.gz"
    snps = tmpdir / f"{prefix}.snps.tsv"
    mask = tmpdir / f"{prefix}.indel_mask.bed"
    snps_masked = tmpdir / f"{prefix}.snps.masked.tsv"
    variant_regions_bed = (
        Path(callable_regions_bed) if callable_regions_bed is not None else regions_bed
    )

    if variant_regions_bed.exists() and variant_regions_bed.stat().st_size > 0:
        call_vcf_from_bam(
            bam=bam,
            regions_bed=variant_regions_bed,
            reference_fasta=reference_fasta,
            out_vcf_gz=vcf,
            min_map_q=min_map_q,
            min_base_q=min_base_q,
        )
        make_indel_mask_bed(vcf, mask, pad_bp=indel_pad_bp)
        extract_snps_table_tsv(vcf, snps)
        mask_snps_table_with_bed(snps, mask, snps_masked)
    else:
        mask.write_text("", encoding="utf-8")
        snps.write_text("", encoding="utf-8")
        snps_masked.write_text("", encoding="utf-8")

    df_regions = read_regions_bed(regions_bed)
    df_snps = read_snps_table(snps_masked)

    df_sites = add_site_metrics(df_snps, min_allele_depth=min_allele_depth)
    df_sites = assign_regions_nonoverlap(df_sites, df_regions)

    df_locus_snps = summarize_loci(
        df_sites,
        max_abs_dp_z_max=max_abs_dp_z_max,
        third_frac_cut=third_frac_cut,
        min_3allele_sites=min_3allele_sites,
        maf_threshold=maf_threshold,
        max_sites_above_maf=max_sites_above_maf,
    )

    # Start from the full regions BED so loci with real coverage but no retained
    # SNPs still appear in the locus metrics table and can contribute has_data.
    df_locus = df_regions[["rid"]].merge(df_locus_snps, on="rid", how="left")

    fill0_int = {
        "n_snps": 0,
        "max_snp_site_depth": 0,
        "n_sites_3allele": 0,
        "n_sites_3allele_strong": 0,
        "n_sites_above_maf": 0,
    }
    fill0_float = {
        "max_snp_depth_z": 0.0,
        "max_frac_3allele": 0.0,
        "max_maf": 0.0,
    }
    for key, value in fill0_int.items():
        if key in df_locus.columns:
            df_locus[key] = df_locus[key].fillna(value).astype("int64")
        else:
            df_locus[key] = value
    for key, value in fill0_float.items():
        if key in df_locus.columns:
            df_locus[key] = df_locus[key].fillna(value).astype(float)
        else:
            df_locus[key] = value
    for flag in ["fail_dp_z", "fail_3allele", "fail_maf"]:
        if flag in df_locus.columns:
            df_locus[flag] = (
                pd.array(df_locus[flag], dtype="boolean").fillna(False).astype(bool)
            )
        else:
            df_locus[flag] = False

    df_softclip = None
    if softclip_len_threshold is not None:
        scl_out = tmpdir / f"{prefix}.softclip_by_region.tsv"
        df_softclip = softclip_fraction_by_region(
            bam=bam,
            regions_bed=regions_bed,
            out_tsv=scl_out,
            softclip_len_threshold=softclip_len_threshold,
        )

    df_locus = add_softclip_and_pass(
        df_locus,
        df_softclip=df_softclip,
        softclip_frac_max=softclip_frac_max,
    )

    # Measure all reads per region separately from the SNP table so the
    # across-sample reducer can distinguish "no SNPs" from "no data".
    cov_out = tmpdir / f"{prefix}.all_reads_by_region.tsv"
    df_all = all_reads_by_region(bam=bam, regions_bed=regions_bed, out_tsv=cov_out)
    df_locus = df_locus.merge(
        df_all[["rid", "all_reads"]],
        on="rid",
        how="left",
        suffixes=("", "_cov"),
    )
    if "all_reads_cov" in df_locus.columns:
        df_locus["all_reads"] = (
            df_locus["all_reads_cov"]
            .fillna(df_locus["all_reads"])
            .fillna(0)
            .astype("int64")
        )
        df_locus = df_locus.drop(columns=["all_reads_cov"])
    else:
        df_locus["all_reads"] = df_locus["all_reads"].fillna(0).astype("int64")
    df_locus["has_data"] = df_locus["all_reads"] > 0

    # Write (pass is in locus metrics)
    write_outputs(df_sites, df_locus, df_regions, prefix=str(tmpdir / prefix))
    return df_sites, df_locus


def _read_bed_to_rids(path: Path) -> set[str]:
    """Read a 3-col BED and return rid set."""
    if not path.exists():
        return set()
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            dtype={"chrom": "string", "start": "int64", "end": "int64"},
        )
    except pd.errors.EmptyDataError:
        return set()
    if df.empty:
        return set()
    rid = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + "-"
        + df["end"].astype(str)
    )
    return set(rid.tolist())


def _read_has_data_from_locus_metrics(path: Path) -> set[str]:
    """Return rid set where has_data is True from one locus_metrics TSV."""
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", dtype={"rid": "string"})
    if "has_data" in df.columns:
        return set(df.loc[df["has_data"].astype(bool), "rid"].astype(str).tolist())
    if "all_reads" in df.columns:
        return set(df.loc[df["all_reads"] > 0, "rid"].astype(str).tolist())
    return set()


def aggregate_across_samples(
    *,
    regions_bed: Path,
    sample_prefixes: Iterable[str],
    in_dir: Path,
    out_prefix: Path,
    fail_frac_max: float,
    min_data_samples: int = 1,
) -> pd.DataFrame:
    """Aggregate per-sample paralog outcomes across shared loci."""
    df_regions = read_regions_bed(regions_bed)
    rid_index = pd.Index(df_regions["rid"].astype(str), name="rid")

    out = pd.DataFrame(index=rid_index)
    out["n_data"] = 0
    out["n_good"] = 0
    out["n_fail"] = 0

    prefixes = list(sample_prefixes)
    if not prefixes:
        raise ValueError("No sample prefixes provided.")

    for prefix in prefixes:
        prefix = str(prefix)
        good_bed = in_dir / f"{prefix}.good.bed"
        bad_bed = in_dir / f"{prefix}.paralog_like.bed"
        locus_tsv = in_dir / f"{prefix}.locus_metrics.tsv"

        # Count each sample only when the locus actually had read data. This is
        # the key hardening needed for no-SNP loci that still had coverage.
        rids_has_data = _read_has_data_from_locus_metrics(locus_tsv)
        rids_good = _read_bed_to_rids(good_bed) & rids_has_data
        rids_bad = _read_bed_to_rids(bad_bed) & rids_has_data

        out.loc[list(rids_has_data), "n_data"] += 1
        out.loc[list(rids_good), "n_good"] += 1
        out.loc[list(rids_bad), "n_fail"] += 1

    denom = out["n_data"].to_numpy(dtype=float)
    out["fail_frac_among_data"] = np.divide(
        out["n_fail"].to_numpy(dtype=float),
        denom,
        out=np.zeros_like(denom),
        where=denom > 0,
    )
    out["good_frac_among_data"] = np.divide(
        out["n_good"].to_numpy(dtype=float),
        denom,
        out=np.zeros_like(denom),
        where=denom > 0,
    )
    out["drop_global"] = (out["n_data"] >= min_data_samples) & (
        out["fail_frac_among_data"] > fail_frac_max
    )
    out["keep_global"] = (out["n_data"] >= min_data_samples) & (~out["drop_global"])

    metrics = df_regions.merge(out.reset_index(), on="rid", how="left")
    metrics.to_csv(f"{out_prefix}.shared_metrics.tsv", sep="\t", index=False)

    keep_rids = set(out.index[out["keep_global"]].astype(str))
    keep_bed = df_regions[df_regions["rid"].astype(str).isin(keep_rids)][
        ["chrom", "start", "end"]
    ]
    keep_bed.to_csv(
        f"{out_prefix}.shared_good.final.bed", sep="\t", header=False, index=False
    )

    strict = (
        (out["n_data"] >= min_data_samples)
        & (out["n_fail"] == 0)
        & (out["n_good"] == out["n_data"])
    )
    strict_rids = set(out.index[strict].astype(str))
    strict_bed = df_regions[df_regions["rid"].astype(str).isin(strict_rids)][
        ["chrom", "start", "end"]
    ]
    strict_bed.to_csv(
        f"{out_prefix}.shared_good.strict_all_samples.bed",
        sep="\t",
        header=False,
        index=False,
    )
    return metrics


def write_per_sample_final_good(
    *,
    sample_prefixes: Iterable[str],
    in_dir: Path,
    shared_good_bed: Path,
    out_dir: Path,
    out_suffix: str = ".final.good.bed",
) -> dict[str, Path]:
    """Write and return per-sample retained BEDs after shared paralog filtering.

    Each returned BED is `sample.good.bed ∩ shared_good.final.bed`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_info = out_dir.parent / "REF_info.txt"
    written: dict[str, Path] = {}
    for prefix in sample_prefixes:
        prefix = str(prefix)
        good_bed = in_dir / f"{prefix}.good.bed"
        out_bed = out_dir / f"{prefix}{out_suffix}"
        tmp_out = out_dir / f"{prefix}{out_suffix}.tmp"
        cmd = [
            BIN_BED,
            "intersect",
            "-sorted",
            "-g",
            str(ref_info),
            "-a",
            str(good_bed),
            "-b",
            str(shared_good_bed),
        ]
        run_pipeline([cmd], tmp_out)
        sort_bed_by_reference_order(tmp_out, out_bed, ref_info)
        tmp_out.unlink(missing_ok=True)
        written[prefix] = out_bed
    return written


if __name__ == "__main__":
    for sname in [
        "SLH_AL_0072-conte",
        "SLH_AL_0064-conte",
        "SLH_AL_0077-conte",
        "SLH_AL_0080-conte",
    ]:
        SNAME = sname

        REF = "/home/deren/Documents/ipyrad-tests/examples/Atub-genome/AmaTu_v01_no00_renamed.fa"
        BED = f"/home/deren/Documents/ipyrad-tests/AMA_2026_OUT/TEST2_tmpdir/beds/{sname}.filtered.bed"
        # BED = "/home/deren/Documents/ipyrad-tests/AMA_2026_OUT/TEST2_tmpdir/loci.bed"
        BAM = f"/home/deren/Documents/ipyrad-tests/AMA_2026_MAP/{sname}.filtered.bam"
        TMP = Path("/tmp")
        # SNPS = TMP / f"{SNAME}.snps.masked.tsv"

        sites, locus = get_sample_paralog_tables(
            bam=BAM,
            regions_bed=BED,
            reference_fasta=REF,
            tmpdir=TMP,
            prefix=SNAME,
            min_map_q=40,
            min_base_q=30,
            indel_pad_bp=10,
            # minor allele must occur this many times to be relevant
            min_allele_depth=2,
            # depth outlier gate (single threshold on max |dp_z| per locus)
            max_abs_dp_z_max=5.0,
            # multi-allelic / third-allele evidence
            third_frac_cut=0.10,
            min_3allele_sites=2,
            # MAF evidence (optional but recommended)
            maf_threshold=0.20,
            max_sites_above_maf=8,
            # softclip evidence
            softclip_len_threshold=20,
            softclip_frac_max=0.25,
        )

        print(locus.head())
