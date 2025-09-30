#!/usr/bin/env python

"""Find the most likely cut-site overhang.

Finds the most common kmers of size 3-8 for the first 20 bp of 10K
reads in an input data file.
"""

from typing import List, Iterator
from pathlib import Path
from collections import Counter
import gzip
from loguru import logger



def iter_reads(fastq: Path, max_len: int, max_reads: int) -> Iterator[bytes]:
    """Generator of sequences (line 2/4) from a fastq file gzipped or not."""
    fastq = Path(fastq)
    open_func = gzip.open if fastq.suffix == ".gz" else open
    with open_func(fastq, 'rb') as inline:
        quart = zip(inline, inline, inline, inline)
        bound = range(max_reads)
        for _, q in zip(bound, quart):
            yield q[1].strip().upper()[:max_len]


def infer_overhang(fastqs: List[Path], max_len: int = 20, max_reads: int = 50_000, anchored: bool = False) -> str:
    """Returns the restriction overhang site in a RAD-seq file.

    This function assumes the fastq data being examined is derived from
    a collection of samples (e.g., undemultiplexed or re-mixed after;
    thus it takes the input as a List of paths) and uses kmers
    (size 3-9) to search for the most common substring in the first
    max_len sites of the first max_reads reads in the input fastqs.
    """
    # If no reads (None or []; e.g., no R2 paths) then return no cutter.
    if not fastqs:
        return ""

    # ensure fastqs is a list
    if isinstance(fastqs, (str, Path)):
        fastqs = [fastqs]

    # {3: [(TCA, 100), (TCG, 100), ...]}
    # {4: [(TCAG, 100), (TCGA, 100), ...]}
    top_counts = {}
    for kmer_size in range(3, 9):
        counts = Counter()
        for fastq in fastqs:
            for read in iter_reads(fastq, max_len, max_reads):
                if anchored:
                    kmer = read[:kmer_size]
                    counts[kmer] += 1
                else:
                    for idx in range(0, len(read) - kmer_size + 1):
                        kmer = read[idx: idx + kmer_size]
                        counts[kmer] += 1
        top_counts[kmer_size] = counts.most_common(10)

    # if every k has exactly 1 observed kmer, return that kmer from the largest k
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
        # (you had a TODO to emit IUPAC; keeping your simple tuple behavior)
        return entries[0][0].decode(), entries[1][0].decode()
    return entries[0][0].decode()


if __name__ == "__main__":

    # ...
    import ipyrad as ip
    ip.set_log_level("TRACE")

    # R1 = "../../examples/Pedic-PE-ddRAD/bella-JJ85-plate_J2_R1.fastq.gz"
    # R2 = "../../examples/Pedic-PE-ddRAD/bella-JJ85-plate_J2_R2.fastq.gz"
    # print(R1, infer_overhang(R1))
    # print(R2, infer_overhang(R2))
    R1s = [
        "../../examples/Pedic-PE-ddRAD/debilis-var-debilior-DE619-1-plate_J2_R1.fastq.gz",
        "../../examples/Pedic-PE-ddRAD/delavayi-JJ18-plate_J2_R1.fastq.gz",
        "../../examples/Pedic-PE-ddRAD/densispica-DE273-plate_J2_R1.fastq.gz",
        "../../examples/Pedic-PE-ddRAD/delavayi-DE545-plate_J2_R1.fastq.gz",
        "../../examples/Pedic-PE-ddRAD/densispica-41491-plate_J2_R1.fastq.gz",
    ]
    # R1 = "../../examples/Pedic-PE-ddRAD/delavayi-DE545-plate_J2_R1.fastq.gz"
    # R2 = "../../examples/Pedic-PE-ddRAD/delavayi-DE545-plate_J2_R2.fastq.gz"
    print('a', infer_overhang(R1s))
    print('x', infer_overhang(R1s, max_reads=50_000, anchored=True))
    # print(R2, infer_overhang(R2))

    # Harder example, overhang includes ambiguous character
    # R1 = "../../pedtest/small_tmp_R1.fastq.gz"
    # R2 = "../../pedtest/small_tmp_R1.fastq.gz"

