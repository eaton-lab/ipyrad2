#!/usr/bin/env python

"""Convert within-sample cluster results to summary tables."""

from pathlib import Path
import pandas as pd
from loguru import logger
from typing import Dict

from .common import (
    CLUSTER_JOINED_SPACER_LEN,
    get_arm_boundary,
    infer_record_type,
)


def get_header_to_seq_dict(consensus_fa: Path) -> Dict[str, str]:
    """Return dict mapping core name to its stripped cluster consensus sequence."""
    seqs = {}
    with open(consensus_fa, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line[0] == ">":
                prefix, mjid, _, _ = line.strip().rsplit(";", 3)
                core = f"{prefix.split('=', 1)[-1]};{mjid}"
            else:
                seqs[core] = line.strip().upper()
    return seqs


def get_header_to_metadata_dict(*fasta_paths: Path) -> Dict[str, tuple[str, int]]:
    """Return record type and arm boundary keyed by raw within-sample input label."""
    out: Dict[str, tuple[str, int]] = {}
    for fasta_path in fasta_paths:
        if not fasta_path.exists():
            continue
        header: str | None = None
        with open(fasta_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    header = line[1:].strip()
                    continue
                if header is None:
                    continue
                record_type = infer_record_type(header)
                _cluster_sequence, arm_boundary = get_arm_boundary(line)
                out[header] = (record_type, int(arm_boundary))
                header = None
    return out


def parse_uc_data(uc_path: Path) -> pd.DataFrame:
    """Parse UC tsv cluster info."""
    rows = []
    with open(uc_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            cat, cidx, _, _, _, _, _, _, label, _ = line.rstrip().split("	")
            if cat not in ("S", "H"):
                continue
            cidx = int(cidx)
            name, mjid, size = label.rsplit(";", 2)
            size = int(size.split("=", 1)[1])
            core = f"{name};{mjid}"
            rows.append((cat, cidx, core, size))
    return pd.DataFrame(rows, columns=["type", "cluster", "core", "size"])


def build_sample_summary(
    sname: str,
    outdir: Path,
    *,
    seed_to_meta: Dict[str, tuple[str, int]] | None = None,
    joined_spacer: int = CLUSTER_JOINED_SPACER_LEN,
) -> pd.DataFrame:
    """Summarize within-sample clusters from vsearch."""
    columns = [
        "sample",
        "cluster_id",
        "seed",
        "length",
        "cluster_length",
        "n_unique",
        "n_reads",
        "record_type",
        "cluster_sequence",
        "arm_boundary",
    ]

    consensus_path = outdir / f"{sname}.consensus.fa"
    uc_path = outdir / f"{sname}.clusters.tsv"

    seed_to_seq = get_header_to_seq_dict(consensus_path)
    if seed_to_meta is None:
        joined_path = outdir / f"{sname}.joined.fa"
        merged_path = outdir / f"{sname}.merged.fa"
        seed_to_meta = get_header_to_metadata_dict(joined_path, merged_path)

    uc = parse_uc_data(uc_path)
    if uc.empty:
        logger.warning(f"No S/H rows parsed in {uc_path}; summary will be empty.")
        df_empty = pd.DataFrame(columns=columns)
        df_empty.to_csv(outdir / f"{sname}.summary.tsv", sep="	", index=False)
        return df_empty

    out_rows = []
    grp = uc.groupby("cluster")
    for i, j in grp:
        nu = len(j['core'].unique())
        nr = sum(j['size'])
        seed = j.loc[j['type'] == "S", 'core'].item()
        cluster_sequence = seed_to_seq[seed]
        record_type, arm_boundary = seed_to_meta.get(
            seed,
            (infer_record_type(seed), len(cluster_sequence)),
        )
        arm_boundary = max(0, min(int(arm_boundary), len(cluster_sequence)))
        has_right_arm = record_type == "joined" and arm_boundary < len(cluster_sequence)
        length = len(cluster_sequence) + (joined_spacer if has_right_arm else 0)
        out_rows.append({
            "sample": sname,
            "cluster_id": str(i),
            "seed": seed,
            "length": length,
            "cluster_length": len(cluster_sequence),
            "n_unique": int(nu),
            "n_reads": int(nr),
            "record_type": record_type,
            "cluster_sequence": cluster_sequence,
            "arm_boundary": int(arm_boundary),
        })
    df = pd.DataFrame(out_rows, columns=columns)
    out_path = outdir / f"{sname}.summary.tsv"
    df.to_csv(out_path, sep="	", index=False)
    logger.debug(f"Wrote {len(df)} clusters → {out_path}")
    return df


def concat_summaries(outdir: Path) -> pd.DataFrame:
    """Concatenate per-sample summary TSVs into one DataFrame."""
    out_tsv = outdir / "concat.summary.tsv"
    tsvs = sorted(outdir.glob("*.summary.tsv"))
    dfs = []
    for p in tsvs:
        df = pd.read_csv(p, sep="	", dtype={
            "sample": "string",
            "cluster_id": "string",
            "length": "Int64",
            "cluster_length": "Int64",
            "seed": "string",
            "n_unique": "Int64",
            "n_reads": "Int64",
            "record_type": "string",
            "cluster_sequence": "string",
            "arm_boundary": "Int64",
        })
        dfs.append(df)
    all_df = pd.concat(dfs, ignore_index=True)
    if out_tsv:
        all_df.to_csv(out_tsv, sep="	", index=False)
        logger.debug("wrote concatenated summaries to {}", out_tsv)
    return all_df
