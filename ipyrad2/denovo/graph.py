

import pandas as pd
from typing import Dict, List, Tuple, Iterable, Set, Iterator
from pathlib import Path
from collections import Counter
from itertools import combinations
from loguru import logger

from ..utils.progress import ProgressBar


def get_edges_dict(outdir) -> Dict[Tuple[str,str], Tuple[float,float]]:
    """Return dict of undirected edges mapped to their highest (id, qcov).

    Returns: Dict
        e.g., {('1A_0;J141', '1B_0;J200'): (0.98, 1.0), ...}
    """
    uc_path = outdir / "global_hits.uc.tsv"
    edges = {}
    with open(uc_path, "rt") as fh:
        for line in fh:
            record = line.rstrip().split("\t")
            # "query+target+id+qstrand+qcov+ql+tl",
            query, target, pid, qstrand, qcov, qlen, tlen = record
            if query == target:
                continue
            # "centroid=1A_0;J141;size=28;seqs=1"
            prefix, mjid, _, _ = query.rsplit(";", 3)
            qcore = f"{prefix.split('=', 1)[-1]};{mjid}"

            prefix, mjid, _, _ = target.rsplit(";", 3)
            tcore = f"{prefix.split('=', 1)[-1]};{mjid}"

            pid = float(pid) / 100.0
            qcov = float(qcov) / 100.0

            # order pair lexicographically)
            edge = tuple(sorted((qcore, tcore)))

            # store if score pid is better than previous if exists
            if edge not in edges:
                edges[edge] = (pid, qcov)
            else:
                if (pid >= edges.get(edge)[0]) & (qcov >= edges[edge][1]):
                    edges[edge] = (pid, qcov)
    return edges


def get_summary_df(outdir: Path) -> pd.DataFrame:
    """Loading a big table with all consensus sequences across all samples.

    Returns: DataFrame
        e.g., [sample, cluster_id,       seed, length, n_unique, n_reads, merged, consensus]
              [  1A_0,          0,  1A_0;J141,    207,        1,      28,  False,  TTGAA...]
              [ ...]
    """
    tsv = outdir / "concat.summary.tsv"
    df = pd.read_csv(tsv, sep="\t")
    return df


def get_clusters_in_graph(nodes: Iterable[str], edges: Dict[Tuple[str,str], Tuple[float,float]]) -> List[Set[str]]:
    """Return connected nodes in the network using a union–find (disjoint-set) structure:

    Returns: List[Set[str]]
        e.g., [{'2G_0;J3802', '1A_0;J3861', '2E_0;J3880', ...}, {...}]
    """
    parent = {n:n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for edge in edges:
        union(*edge)
    comps = {}
    for n in nodes:
        r = find(n)
        comps.setdefault(r, set()).add(n)
    return list(comps.values())


def _duplicates_present(comp: Set[str]) -> bool:
    """Return True if any sample has >1 consensus in a cluster."""
    c = Counter(s.split(";", 1)[0] for s in comp)
    return max(c.values()) > 1


def _get_subset_edges_above_threshold(
                                    nodes: Set[str],
                                    edict: Dict[Tuple[str, str], Tuple[float, float]],
                                    threshold: float) -> Dict[Tuple[str, str],Tuple[float,float]]:
    """Return dict mapping {edge: (pid, qcov} only for edges with pid >= t and for nodes in nodes."""
    out = {}
    for u, v in combinations(nodes, 2):
        a, b = sorted((u, v))
        # Get returns None if the edge doesn't exist
        val = edict.get((a, b))
        if val is None:
            continue
        pid, qcov = val
        if pid >= threshold:
            out[(a, b)] = (pid, qcov)
    return out


def iter_non_duplicated_subcomponent_by_ascending_pid(component: Set[str], edict: Dict[Tuple[str, str], Tuple[float, float]]) -> Iterator[Set[str]]:
    """
    For the given component (set of seeds), find the *lowest* pairwise identity
    threshold (pid) such that if we keep only edges with pid >= t, all resulting
    subcomponents have ≤ max_per_sample seeds per sample.

    Returns:
        (subcomponents, chosen_threshold)
    """
    if len(component) <= 1:
        if component:
            yield set(component)
        return

    # get pids for all component edges
    pids = set()
    for u, v in combinations(component, 2):
        a, b = sorted((u, v))
        val = edict.get((a, b))
        if val is not None:
            pids.add(val[0])

    # Consider thresholds in ascending order (least strict → more strict),
    # plus an extra sentinel above the max to force full disconnection if needed.
    # We want the minimal t that satisfies the per-sample rule.
    remaining = set(component)  # make a copy to pop from.
    thresholds = sorted(pids)
    while 1:
        # try the next lowest threshold
        try:
            thresh = thresholds.pop(0)
        except IndexError:
            break

        # get subset of edges >= threshold: {(u,v):(pid,qcov),...}
        e = _get_subset_edges_above_threshold(remaining, edict, thresh)

        # get components connected by e: [{a1, b1, c1}, {a2, b2, c2, a3, b3, c3}]
        subcomps = get_clusters_in_graph(remaining, e)

        # yield components that don't contain duplicates
        to_remove = set()
        for sub in subcomps:
            if not _duplicates_present(sub):
                yield sub
                to_remove.update(sub)
            else:
                logger.debug(f"duplicates in {sub}")
        remaining.difference_update(to_remove)

        # end loop if all nodes have been assigned to subgroups
        if not remaining:
            break

    # only way to prevent duplicates is to treat each remaining consensus
    # as its own locus:
    for node in remaining:
        yield {node}



def make_global_tables(outdir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parse global_hits.tsv to contruct graph; split graph on increasing
    identify thresholds to split paralogs; build dataframes with stats
    for each locus from consensus summary table.

    Returns:
      mapping_df: rows (locus_id, core, sample, n_reads, length, gc, frac_N)
      depth_wide: rows loci, columns depth.{sample}, alleles.{sample}
    """
    # write to true outdir (parent of this dir)
    out_mapping_tsv = outdir.parent / "loci.mapping.tsv"
    out_stats_tsv = outdir.parent / "loci.stats.tsv"

    # load per-sample summary stats
    df = get_summary_df(outdir)
    df = df.set_index("seed")

    # edges → components
    edict = get_edges_dict(outdir)
    comps = get_clusters_in_graph(df.index, edict)

    # refinement to split paralogs. Note that our goal here is to create
    # a denovo reference, so we want to include paralogs so ensure that
    # reads will later map to their respective separate copies.
    prog = ProgressBar(len(comps), 0, "Splitting paralogs")
    prog.finished = 0
    prog.update()

    refined_parts = []
    for part in comps:
        for sub in iter_non_duplicated_subcomponent_by_ascending_pid(part, edict):
            refined_parts.append(sub)

        prog.finished += 1
        prog.update()
    print("")

    # n_component_splits = len(refined_parts) - len(comps)
    logger.info(f"split {len(comps)} clusters into {len(refined_parts)} non-duplicated subclusters")

    # mapping rows ----------------------------------------------
    logger.info(f"writing mapping table to {out_mapping_tsv}")

    mapping_rows = []
    for k, part in enumerate(refined_parts, start=1):
        # iterate over each consensus in the locus
        for core in sorted(part):
            # store consensus id info
            info = df.loc[core]
            row = {
                "locus": int(k),
                "sample": str(core.rsplit(";", 1)[0]),
                "n_reads": int(info["n_reads"]),
                "n_unique": int(info["n_unique"]),
                "length": int(info["length"]),
                "merged": int(info["merged"]),
                "cluster_id": int(info["cluster_id"]),
                "core": str(core),
            }
            mapping_rows.append(row)
    mapping_df = pd.DataFrame(mapping_rows)
    mapping_df.to_csv(out_mapping_tsv, sep="\t", float_format="%12.6f", index=False)

    # depth wide table
    logger.info(f"writing locus stats to {out_stats_tsv}")

    recs = []
    for locus, sub in mapping_df.groupby("locus"):
        rec = {
            "locus": locus,
            "n_samples": int(sub.shape[0]),
            "n_reads_sum": int(sub["n_reads"].sum()),
            "n_reads_mean": float(sub["n_reads"].mean()),
            "n_reads_std": float(sub["n_reads"].std()),
            "length_mean": float(sub["length"].mean()),
            "length_std": float(sub["length"].std()),
            "merged_freq": float(sub["merged"].mean()),
            "samples": ",".join([str(i) for i in sub["sample"]])
        }
        recs.append(rec)
    stats = pd.DataFrame(recs).sort_values("locus").reset_index(drop=True)
    stats.to_csv(out_stats_tsv, sep="\t", float_format="%12.6f", index=False)
    return mapping_df, stats



if __name__ == "__main__":

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max.rows', None)

    outdir = Path("/home/deren/Documents/ipyrad-tests/WMERGE_DENOVO/tmpdir")

    # load per-sample summary stats
    df = get_summary_df(outdir)
    df = df.set_index("seed")
    # print(df.sort_values(by="n_reads", ascending=True).head(20))
    print(df.sort_values(by=["cluster_id", "sample"], ascending=True).iloc[:, :6].head())
    make_global_tables(outdir)
