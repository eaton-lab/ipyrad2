#!/usr/bin/env python

"""Find the most likely restriction overhang from kmer analysis.

Finds the most common kmers of size 3-8 for the first 20 bp of ~10K
reads in an input data file.
"""


from typing import List, Iterator
from pathlib import Path
from collections import Counter
import gzip
from loguru import logger
from .exceptions import IPyradError
from .parallel import run_with_pool



def iter_reads(fastq: Path, max_len: int, max_reads: int) -> Iterator[bytes]:
    """Generator of sequences (line 2/4) from a fastq file gzipped or not."""
    fastq = Path(fastq)
    open_func = gzip.open if fastq.suffix == ".gz" else open
    with open_func(fastq, 'rb') as inline:
        quart = zip(inline, inline, inline, inline)
        bound = range(max_reads)
        for _, q in zip(bound, quart):
            yield q[1].strip().upper()[:max_len]


def get_kmer_counts(fastq: Path, kmer_size: int, max_len: int, max_reads: int):
    """Return kmer counts for N reads."""
    counts = Counter()
    for read in iter_reads(fastq, max_len, max_reads):
        for idx in range(0, len(read) - kmer_size + 1):
            kmer = read[idx: idx + kmer_size]
            counts[kmer] += 1
    return counts


def get_overhang_from_kmers(fastqs: List[Path], max_len: int, max_reads: int, workers: int, log_level: str):
    """..."""
    # get N kmer counts sampled evenly across files
    max_reads_per_file = int(max_reads / len(fastqs))
    jobs = {}
    for kmer_size in range(3, 9):
        for fastq in fastqs:
            kwargs = dict(
                fastq=fastq,
                kmer_size=kmer_size,
                max_len=max_len,
                max_reads=max_reads_per_file,
            )
            jobs[(kmer_size, fastq.name)] = (get_kmer_counts, kwargs)

    # fetch counters in parallel
    kcounts = run_with_pool(jobs, log_level, workers)

    # combine results and store top 10 at each size
    top_counts = {}
    for kmer_size in range(3, 9):
        kc = Counter()
        for i in fastqs:
            kc.update(kcounts[(kmer_size, i.name)])
        top_counts[kmer_size] = kc.most_common(10)

    # if every kmer_size has only 1 observed kmer, return that kmer from the largest k
    single_k = [k for k, lst in top_counts.items() if len(lst) == 1]
    if single_k:
        k = max(single_k)
        return top_counts[k][0][0].decode()

    # Compare ratios the most common to next most to find the ...
    # If True is AAA but kmer_size=2 then
    # XXAA, AAXX will create 32 equally likely codes
    # If True is AAA but kmer_size=4 then
    # XAAA, AAAX will create 8 equally likely codes
    # So find the kmer_size that minimizes the ratios of alternative
    # kmers to the most common kmer to find the optimal size. THen
    # return the most frequent kmer at that size.
    ratios = {}
    for k, lst in top_counts.items():
        if not lst:
            ratios[k] = float("inf")
            continue
        max_count = lst[0][1]
        # sum of the next two (if present), normalized by the max
        denom = max(1, max_count)
        ratios[k] = sum(x[1] / denom for x in lst[1:3])
        logger.trace(f"K={k}: {[round(i[1] / denom, 3) for i in lst]}")

    logger.trace(f"best Ks: {ratios}")
    best_k = min(ratios, key=ratios.get)

    # FIX: guard when there is only one entry for best_k
    entries = top_counts[best_k]
    if not entries:
        return ""  # nothing counted at all
    if len(entries) == 1:
        return entries[0][0].decode()

    # ambiguous barcodes will have top 2 at near equal frequency
    top_count = entries[0][1]
    sec_count = entries[1][1]
    if top_count and (sec_count / top_count) > 0.90:
        # TODO: handle IUPAC
        raise IPyradError("ambiguous support for restriction overhang. Disable auto-infer-re-overhangs and set manually.")
        return entries[0][0].decode(), entries[1][0].decode()
    return entries[0][0].decode()


if __name__ == "__main__":

    DIR = Path("/home/deren/Documents/ipyrad-tests/examples/Ama-PE-ddRAD/")
    R1s = list(DIR.glob("SLH_AL*_R2*"))
    x = get_overhang_from_kmers(R1s, 18, 500_000, 10)
    print(x)