#!/usr/bin/env python

"""Convert within-sample cluster results to summary tables.
"""

from pathlib import Path
import pandas as pd
from loguru import logger
from typing import Dict


def get_header_to_seq_dict(consensus_fa: Path) -> Dict[str, str]:
    """Return dict mapping core name to (sequence, mj)."""
    seqs = {}
    with open(consensus_fa, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                # >centroid=1A_0;J16321;size=25;seqs=1
                prefix, mjid, _, _ = line.strip().rsplit(";", 3)
                core = f"{prefix.split("=", 1)[-1]};{mjid}"
            else:
                seqs[core] = line.strip().upper()
    return seqs


def parse_uc_data(uc_path: Path) -> pd.DataFrame:
    """Parse UC tsv cluster info."""
    rows = []
    with open(uc_path, "rt") as fh:
        for line in fh:
            # S       4       207     *       *       *       *       *       1A_0;J16321;size=25     *
            cat, cidx, _, _, _, _, _, _, label, _ = line.rstrip().split("\t")
            if cat not in ("S", "H"):    # ignore C/L/D/U/…
                continue
            cidx = int(cidx)
            # e.g., 1A_0;J;14481;size=23
            name, mjid, size = label.rsplit(";", 2)
            size = int(size.split("=", 1)[1])
            core = f"{name};{mjid}"
            rows.append((cat, cidx, core, size, 1 if mjid[0] == "M" else 0))
    return pd.DataFrame(rows, columns=["type", "cluster", "core", "size", "merged"])


def build_sample_summary(sname: str, outdir: Path, joined_spacer: int = 24) -> pd.DataFrame:
    """Summarize within-sample clusters from vsearch.
    """
    columns = ["sample", "cluster_id", "seed", "length", "n_unique", "n_reads", "merged", "consensus"]

    # paths
    consensus_path = outdir / f"{sname}.consensus.fa"
    uc_path = outdir / f"{sname}.clusters.tsv"

    # get {header: sequence}
    seed_to_seq = get_header_to_seq_dict(consensus_path)

    # get df w/ rows ["type", "cluster", "core", "size", "merged"]
    # example: ['S', 10, '1A_0;M200', 22, 0]
    uc = parse_uc_data(uc_path)
    if uc.empty:
        logger.warning(f"No S/H rows parsed in {uc_path}; summary will be empty.")
        df_empty = pd.DataFrame(columns=columns)
        df_empty.to_csv(outdir / f"{sname}.summary.tsv", sep="\t", index=False)
        return df_empty

    # per-cluster totals
    out_rows = []
    grp = uc.groupby("cluster")
    for i, j in grp:
        # logger.warning(f"entering {sname} cluster {i}")
        nu = len(j['core'].unique())
        nr = sum(j['size'])
        seed = j.loc[j['type'] == "S", 'core'].item()
        consensus = seed_to_seq[seed]
        out_rows.append({
            "sample": sname,
            "cluster_id": str(i),
            "seed": seed,
            "length": len(consensus),
            "n_unique": int(nu),
            "n_reads": int(nr),
            "merged": seed.rsplit(";", 1)[1][0] == "M",
            "consensus": consensus,
        })
    df = pd.DataFrame(out_rows, columns=columns)
    out_path = outdir / f"{sname}.summary.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    logger.info(f"Wrote {len(df)} clusters → {out_path}")
    return df


def concat_summaries(outdir: Path) -> pd.DataFrame:
    """Concatenate per-sample summary TSVs into one DataFrame.
    """
    out_tsv = outdir / "concat.summary.tsv"
    tsvs = sorted(outdir.glob("*.summary.tsv"))
    dfs = []
    for p in tsvs:
        df = pd.read_csv(p, sep="\t", dtype={
            "sample": "string",
            "cluster_id": "string",
            "length": "Int64",
            "seed": "string",
            "n_unique": "Int64",
            "n_reads": "Int64",
            "merged": "boolean",
            "consensus": "string",
        })
        dfs.append(df)
    all_df = pd.concat(dfs, ignore_index=True)
    if out_tsv:
        all_df.to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"Wrote concatenated summaries → {out_tsv}")
    return all_df


if __name__ == "__main__":

    pd.set_option('display.max_columns', None)
    DIR = Path("/home/deren/Documents/ipyrad-tests/WMERGE_DENOVO/")
    # df = build_sample_summary("1A_0", DIR)
    # print(df)

    uc = parse_uc_data(DIR / "1A_0.clusters.tsv")
    for i, j in uc.groupby("cluster"):
        print(f"{i}\n{j}")