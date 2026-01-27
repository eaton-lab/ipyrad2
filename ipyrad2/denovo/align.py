#!/usr/bin/env python

"""Align consensus sequences in a component using MAFFT and write the
consensus to denovo_reference.fa.

"""

from __future__ import annotations
from typing import Dict, List, Tuple, Iterator, Callable, Any
import sys
from pathlib import Path
import pandas as pd
from loguru import logger
from ipyrad2.utils.seqs import IUPAC
from ipyrad2.utils.parallel import run_pipeline, run_with_pool_iter
from ..utils.progress import ProgressBar



BIN = Path(sys.prefix) / "bin"
BIN_MAFFT = str(BIN / "mafft")


# ----------------- small utilities -----------------

def fasta_text(record: List[Tuple[str,str]]) -> str:
    """Render a FASTA string from [(name, seq), ...] with no line wrapping."""
    return "".join(f">{h}\n{s}\n" for h, s in record)


def strip_joined_spacer(seq: str, spacer_len: int = 24) -> str:
    """Remove exactly one N^spacer_len joiner (if present) from a joined read."""
    s = seq.upper()
    token = "N" * spacer_len
    return s.replace(token, "", 1) if token in s else s


def consensus_from_aligned(aligned: List[Tuple[str,str]], min_prop: float = 0.5) -> str:
    """Build a simple gap-aware consensus from an alignment.

    Rules
    -----
    - Ignore gaps and non-ACGT at each column.
    - If the most frequent base has frequency >= `min_prop`, emit that base.
    - Otherwise:
        * If exactly two bases are present, emit the IUPAC code for that pair.
        * Else emit 'N'.
    - All-gap columns are dropped.

    Parameters
    ----------
    aligned
        List of (name, aligned_sequence) with equal lengths.
    min_prop
        Majority threshold per column (0–1).

    Returns
    -------
    str
        Consensus sequence (ungapped).
    """
    if not aligned:
        return ""
    L = len(aligned[0][1])
    if any(len(s) != L for _, s in aligned):
        raise ValueError("Aligned sequences differ in length.")
    out = []
    for j in range(L):
        bases = [s[j].upper() for _, s in aligned if s[j] != "-" and s[j].upper() in "ACGT"]
        if not bases:
            continue
        counts = {b: bases.count(b) for b in "ACGT"}
        tot = sum(counts.values())
        base, n = max(counts.items(), key=lambda kv: kv[1])
        if n / max(1, tot) >= min_prop:
            out.append(base)
        else:
            alleles = [b for b, c in counts.items() if c > 0]
            out.append(IUPAC.get(frozenset(alleles), "N") if len(alleles) == 2 else "N")
    return "".join(out)


# ----------------- MAFFT runner -----------------


def mafft_align_one(record: List[Tuple[str, str]], threads: int = 1) -> List[Tuple[str, str]]:
    """Run MAFFT on records via stdin and return aligned [(name, aligned_seq), ...]."""
    if not record:
        raise RuntimeError("record is empty.")
    # run MAFFT alignment
    argv = [
        BIN_MAFFT,
        "--quiet",
        "--auto",
        "--thread", str(threads),
        "--adjustdirection",
        "-",
    ]
    rc, out, err = run_pipeline([argv], outfile=None, stdin_text=fasta_text(record))
    if rc != 0 or not out:
        raise RuntimeError(f"mafft failed rc={rc} stderr={err.decode('utf-8', 'replace')}")

    # parse aligned FASTA
    aln = []
    name, buf = None, []
    for ln in out.decode("utf-8", "replace").splitlines():
        if ln.startswith(">"):
            if name is not None:
                aln.append((name, "".join(buf)))
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name is not None:
        aln.append((name, "".join(buf)))
    return aln


# ----------------- Worker & pool driver -----------------


def worker_build_consensus(
    locus_id: int,
    record: List[Tuple[str, str]],
    min_prop: float = 0.5,
    threads: int = 1,
) -> Tuple[int, str]:
    """MAFFT worker: compute consensus for one locus.

    Returns
    -------
    (locus_id, consensus)
    """
    if not record:
        raise RuntimeError(f"no record for locus {locus_id}")
    if len(record) == 1:
        return locus_id, record[0][1]
    aln = mafft_align_one(record, threads=threads)
    cons = consensus_from_aligned(aln, min_prop=min_prop)
    return locus_id, cons


# ---------- 1) Stream locus -> [(seed, seq), ...] from DataFrames ----------


def iter_locus(
    mapping_tsv: Path,
    summary_tsv: Path,
    spacer_len: int = 24,
) -> Iterator[Tuple[int, List[Tuple[str,str]]]]:
    """Yield (locus_id, [(name, sequence), ...]) from mapping+summary TSVs.
    """
    # 1) Load mapping fully
    mapping = pd.read_csv(mapping_tsv, sep="\t", dtype=str)
    mapping = mapping[["locus", "core"]].dropna()
    mapping["locus"] = mapping["locus"].astype(int)
    mapping["core"] = mapping["core"].astype(str)

    # 2) Stream summary to build a minimal seed->sequence dict
    name_to_seq = {}
    with open(summary_tsv, 'rt') as indata:
        # skip first line header
        _header = indata.readline()
        for line in indata:
            _, _, core, _, _, _, _, seq = line.strip().split("\t")
            name_to_seq[core] = seq

    # 4) Group by locus and yield one locus at a time
    for cluster_id, sub in mapping.groupby("locus", sort=False):
        seeds = sub["core"].tolist()
        recs: List[Tuple[str, str]] = []
        for name in seeds:
            seq = name_to_seq[name]
            seq = strip_joined_spacer(seq, spacer_len)
            recs.append((name, seq))
        yield int(cluster_id), recs


def iter_mafft_jobs(
    mapping_tsv: Path,
    summary_tsv: Path,
    spacer_len: int,
    min_prop: float,
    threads: int,
) -> Iterator[Tuple[int, Tuple[Callable[..., Any], Dict[str, Any]]]]:
    """Yield MAFFT jobs in mapping order for ordered streaming.

    Yields
    ------
    (submit_index, (worker_build_consensus, kwargs))
    """
    submit_idx = 0
    for locus_id, record in iter_locus(mapping_tsv, summary_tsv, spacer_len=spacer_len):
        kwargs = dict(
            locus_id=locus_id,
            record=record,
            min_prop=min_prop,
            threads=threads,
        )
        yield submit_idx, (worker_build_consensus, kwargs)
        submit_idx += 1


def write_ordered_consensus_stream_to_file(
    outdir: Path,
    log_level: str = "INFO",
    spacer_len: int = 24,
) -> None:
    """Run MAFFT consensuses in parallel and write FASTA in mapping order.

    Parameters
    ----------
    out_fa
        Output FASTA path.
    log_level
        Log level passed to the pool initializer.
    spacer_len
        Length of N-joiner to strip before alignment.
    """
    # paths
    mapping_tsv = outdir / "loci.mapping.tsv"
    summary_tsv = outdir / "tmpdir" / "concat.summary.tsv"
    out_fa = outdir / "denovo_reference.fa"

    # Hack for getting the number of loci. Read in chunks and only parse the
    # first column, makes it faster.
    nloci = sum(len(chunk) for chunk in pd.read_csv(mapping_tsv, chunksize=10000, usecols=[0]))
    prog = ProgressBar(nloci, 0, "Writing denovo reference sequence")
    prog.finished = 0
    prog.update()

    # open outfile and write loci
    with open(out_fa, "wt") as fh:
        for lid, records in iter_locus(mapping_tsv, summary_tsv, spacer_len=spacer_len):
            # Take the longest sequence in this cluster as the representative
            cons = max([x[1] for x in records], key=len)
            fh.write(f">locus_{lid}\n{cons}\n")
        prog.finished += 1
        prog.update()
    print("")
    logger.info(f"wrote denovo reference to {out_fa}")


# Old-style function that produced consensus sequence from MAFFT aligned
# loci. We are skipping all the overhead and just taking the longest
# record per locus as the 'representative'. Keeping this here as a all-back
# but can be deleted when we are comfortable with the simplified function.
# IAO 12/11/25
#def write_ordered_consensus_stream_to_file(
#    outdir: Path,
#    log_level: str = "INFO",
#    max_workers: int = 6,
#    threads: int = 1,
#    spacer_len: int = 24,
#    min_prop: float = 0.5,
#) -> None:
#    """Run MAFFT consensuses in parallel and write FASTA in mapping order.
#
#    Parameters
#    ----------
#    mapping_tsv
#        Path to mapping table (cols: locus, core).
#    summary_tsv
#        Path to summary table; assumes known column layout with core/seq.
#    out_fa
#        Output FASTA path.
#    log_level
#        Log level passed to the pool initializer.
#    max_workers
#        Number of worker processes.
#    threads
#        Threads per MAFFT subprocess (keep max_workers*threads ≤ cores).
#    spacer_len
#        Length of N-joiner to strip before alignment.
#    min_prop
#        Majority threshold per column for consensus.
#    """
#    # paths
#    mapping_tsv = outdir / "loci.mapping.tsv"
#    summary_tsv = outdir / "tmpdir" / "concat.summary.tsv"
#    out_fa = outdir / "denovo_reference.fa"
#
#    # generator of: (id, (func, data))
#    jobs_it = iter_mafft_jobs(mapping_tsv, summary_tsv, spacer_len, min_prop, threads)
#    buffer: Dict[int, Tuple[str, str]] = {}
#    next_key = 0
#
#    # open outfile and write loci
#    with open(out_fa, "wt") as fh:
#        # feed jobs from generator to pool and write as they finish
#        kwargs = dict(jobs_iter=jobs_it, log_level=log_level, max_workers=max_workers)
#        for key, (locus_id, consensus) in run_with_pool_iter(**kwargs):
#            buffer[key] = (locus_id, consensus)
#            # Flush in-order results as far as we can
#            while next_key in buffer:
#                lid, cons = buffer.pop(next_key)
#                fh.write(f">locus_{lid}\n{cons}\n")
#                next_key += 1
#    logger.info(f"wrote denovo reference to {out_fa}")


if __name__ == "__main__":

    import pandas as pd
    outdir = Path("/home/deren/Documents/ipyrad-tests/WMERGE_DENOVO/tmpdir")
    mapping_tsv = outdir.parent / "loci.mapping.tsv"
    summary_tsv = outdir / "concat.summary.tsv"
    out = outdir.parent / "denovo_reference.fa"
    write_ordered_consensus_stream_to_file(mapping_tsv, summary_tsv, out, max_workers=11)

    # ilocus = iter_locus(mapping_tsv, summary_tsv)

    # for i, locus in ilocus:
    #     print(f"\n{locus}\n")

    # # 1) Load your big summary table (seed -> sequence)
    # summary_df = pd.read_csv(outdir / "concat.summary.tsv", sep="\t").set_index("seed", drop=False)

    # # 2) mapping_df from your make_global_tables() call
    # # mapping_df = ...  # already in memory, or read it back:
    # mapping_df = pd.read_csv(outdir.parent / "loci.mapping.tsv", sep="\t")

    # # 3) Pack per-locus records: { locus_id: [(seed, sequence), ...], ... }
    # loci_to_records = {}
    # for locus, sub in mapping_df.groupby("locus", sort=True):
    #     recs = []
    #     for seed in sub["core"].astype(str):
    #         seq = summary_df.at[seed, "consensus"] if seed in summary_df.index else None
    #         if isinstance(seq, str) and seq:
    #             recs.append((seed, seq))
    #     if recs:
    #         loci_to_records[str(locus)] = recs

    # run_consensus_pool(loci_to_records, outdir.parent / "denovo_reference.fa", )
