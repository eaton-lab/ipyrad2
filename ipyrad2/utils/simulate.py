#!/usr/bin/env python

from __future__ import annotations
from typing import List, Tuple, Dict, Iterable, Optional
from dataclasses import dataclass
import re
import math
import numpy as np
from ipyrad2.utils.seqs import revcomp


@dataclass(frozen=True)
class Enzyme:
    """..."""
    name: str
    site: str           # recognition sequence (e.g., "CTGCAG")
    cut_offset: int     # 0-based index from the start of the site on the forward strand
                        # (i.e., where the phosphodiester cut occurs along the reference strand)


def generate_scaffold(length: int, seed: Optional[int] = None) -> str:
    """Generate a random ATCG scaffold using NumPy."""
    rng = np.random.default_rng(seed)
    bases = np.array(list("ATCG"))
    return "".join(rng.choice(bases, size=length))


def _find_cut_positions(seq: str, enz: Enzyme) -> List[int]:
    """
    Return genomic cut positions (0..len(seq)) for an enzyme.
    We search for the motif on the forward strand; if the motif is non-palindromic,
    we also search for the reverse-complement motif and map its cut to forward coords.
    """
    L = len(seq)
    motif = enz.site.upper()
    rc = revcomp(motif)
    mlen = len(motif)

    # compile simple literal regex (escape just in case)
    pat_fwd = re.compile(re.escape(motif))
    pat_rev = re.compile(re.escape(rc)) if rc != motif else None

    cuts: List[int] = []

    # forward-strand motif: cut at start + offset
    for m in pat_fwd.finditer(seq):
        cuts.append(m.start() + enz.cut_offset)

    # reverse-strand motif: map cut to forward-coordinates.
    # If the enzyme cuts at offset k from the *forward* motif start,
    # then on the reverse motif starting at r, the corresponding
    # forward coordinate is: r + (mlen - enz.cut_offset)
    if pat_rev:
        for m in pat_rev.finditer(seq):
            # position on forward strand aligned to reverse motif's cut
            pos = m.start() + (mlen - enz.cut_offset)
            cuts.append(pos)

    # Keep only cuts strictly within the molecule (avoid negative/out-of-range)
    cuts = [p for p in cuts if 0 < p < L]
    return cuts


def digest_and_collect(
    seq: str,
    enzymes: Iterable[Enzyme],
    size_range: Tuple[int, int],
    require_mixed_ends: bool = True,
) -> List[Dict]:
    """
    Digest the sequence with the given enzymes and return fragments in size_range.
    If require_mixed_ends=True, only keep fragments whose left and right cuts
    come from *different* enzymes (ddRAD-style).
    Returns a list of dicts with coordinates, length, end enzyme names, and sequence.
    """
    # find cuts for each enzyme, annotate with enzyme index/name
    all_cuts: List[Tuple[int, str]] = []
    for e in enzymes:
        for p in _find_cut_positions(seq, e):
            all_cuts.append((p, e.name))

    # Sort and de-duplicate cuts (some motifs can overlap and yield same cut)
    all_cuts = sorted(set(all_cuts), key=lambda x: x[0])

    # Build fragments only between *enzyme* cuts (exclude 0 and L unless they are cuts)
    frags: List[Dict] = []
    for (left_pos, left_e), (right_pos, right_e) in zip(all_cuts[:-1], all_cuts[1:]):
        length = right_pos - left_pos
        if size_range[0] <= length <= size_range[1]:
            if (not require_mixed_ends) or (left_e != right_e):
                frags.append({
                    "start": left_pos,
                    "end": right_pos,
                    "length": length,
                    "left_enzyme": left_e,
                    "right_enzyme": right_e,
                    "seq": seq[left_pos:right_pos],
                })
    return frags

# ---------- convenience wrapper for ddRAD simulation ----------

def simulate_ddrad(
    genome_len: int,
    enzymes: Iterable[Enzyme],
    size_range: Tuple[int, int] = (200, 500),
    seed: Optional[int] = None,
    require_mixed_ends: bool = True,
) -> Tuple[str, List[Dict]]:
    """
    1) Generate a random genome/scaffold.
    2) Digest with enzymes.
    3) Return (scaffold, filtered_fragments).
    """
    scaffold = generate_scaffold(genome_len, seed=seed)
    frags = digest_and_collect(
        scaffold, enzymes, size_range=size_range, require_mixed_ends=require_mixed_ends
    )
    return scaffold, frags



# ---------------- adapters/barcodes + read simulation ----------------

@dataclass(frozen=True)
class Adapters:
    """Literal DNA sequences for adapters. If you include {index} in i7/i5,
    it will be formatted with the index string (e.g., 'ACGT')."""
    P5_i5: str  # sequence for the P5/i5 side (R1 side)
    P7_i7: str  # sequence for the P7/i7 side (R2 side)

def ligate_adapters(
    fragments: List[Dict],
    *,
    p5_is_left: bool,
    adapters: Adapters,
    inline_barcodes: Optional[List[str]] = None,
    i7_indices: Optional[List[str]] = None,
    i5_indices: Optional[List[str]] = None,
    rng: Optional[np.random.Generator] = None,
) -> List[Dict]:
    """
    Add adapters (and optional inline barcodes) to each fragment.
    - p5_is_left: True if the library prep places P5 (R1) on the *left* cut.
                  For ddRAD, this depends on which enzyme end got which adapter.
    - inline_barcodes: if provided, an inline tag will be inserted on the R1 side
      *between* the P5 adapter and the insert (i.e., in-read barcode).
    - i7/i5 indices: if provided, they format the adapter sequences by replacing
      '{index}' with the chosen index. This is mostly to carry metadata and allow
      read-through to hit index if a read overruns the insert.
    Returns new fragment dicts with keys:
      'seq_adapter_plus', 'r1_side', 'r2_side', 'inline_bc', 'i7', 'i5'
    """
    if rng is None:
        rng = np.random.default_rng()

    out: List[Dict] = []
    for frag in fragments:
        insert = frag["seq"]

        # Choose indices / inline bc if lists provided; otherwise empty.
        inline_bc = (rng.choice(inline_barcodes) if inline_barcodes else "")
        i7 = (rng.choice(i7_indices) if i7_indices else "")
        i5 = (rng.choice(i5_indices) if i5_indices else "")

        P5 = adapters.P5_i5.format(index=i5)
        P7 = adapters.P7_i7.format(index=i7)

        # Decide which end is R1/R2 (P5=R1, P7=R2)
        if p5_is_left:
            r1_side, r2_side = "left", "right"
            seq_with_adapters = P5 + inline_bc + insert + P7
        else:
            r1_side, r2_side = "right", "left"
            # If P5 on the right, inline barcode still belongs to R1 side,
            # so it goes adjacent to insert on the right edge.
            seq_with_adapters = P7 + insert + inline_bc + P5

        out.append({
            **frag,
            "seq_adapter_plus": seq_with_adapters,
            "r1_side": r1_side,
            "r2_side": r2_side,
            "inline_bc": inline_bc,
            "i7": i7,
            "i5": i5,
        })
    return out

def _introduce_substitution_errors(seq: str, p_error: float, rng: np.random.Generator) -> str:
    if p_error <= 0:
        return seq
    if p_error >= 1:
        # pathological, but handle
        return "".join({'A':'C','C':'A','G':'T','T':'G'}[b] for b in seq)

    bases = np.array(list(seq))
    n = len(bases)
    # Flip error coins
    flips = rng.random(n) < p_error
    if not flips.any():
        return seq

    # For error positions, choose a base != original
    alphabet = np.array(list("ACGT"))
    for i in np.where(flips)[0]:
        choices = alphabet[alphabet != bases[i]]
        bases[i] = rng.choice(choices)
    return "".join(bases.tolist())

def _introduce_homopolymer_slips(seq: str, p_slip: float, rng: np.random.Generator) -> str:
    """
    With probability p_slip *per homopolymer run*, randomly insert or delete
    a single base within runs of length >= 3. Keeps things simple and rare,
    as Illumina mainly shows substitutions; this adds occasional small INDELs.
    """
    if p_slip <= 0:
        return seq
    out = []
    i = 0
    L = len(seq)
    while i < L:
        j = i + 1
        while j < L and seq[j] == seq[i]:
            j += 1
        run_base = seq[i]
        run_len = j - i
        # copy the run
        run = [run_base] * run_len
        if run_len >= 3 and (rng.random() < p_slip):
            if rng.random() < 0.5:
                # deletion of one base from the run (if >1)
                if run_len > 1:
                    run.pop(rng.integers(0, len(run)))
            else:
                # insertion of one more base of the same type
                ins_pos = rng.integers(0, len(run)+1)
                run.insert(ins_pos, run_base)
        out.extend(run)
        i = j
    return "".join(out)

def _extract_reads_from_oriented_fragment(
    frag: Dict,
    read_len: int,
    rng: np.random.Generator,
) -> Tuple[str, str]:
    """
    Given an adapter-ligated fragment (seq_adapter_plus, r1_side/r2_side),
    return (R1, R2) sequences of fixed length read_len.
    - If read overruns the insert, it will continue into adapter/inline-bc.
    - R2 is reverse-complemented of the second read (Illumina convention).
    """
    s = frag["seq_adapter_plus"]
    n = len(s)
    if frag["r1_side"] == "left":
        r1 = s[:read_len]
        r2_template = s[-read_len:]  # take from right side
    else:
        r1 = s[-read_len:]
        r2_template = s[:read_len]
    # R2 must be the reverse-complement of the template piece
    r2 = revcomp(r2_template)
    # Pad if needed (very short fragments)
    if len(r1) < read_len:
        r1 = r1 + "N" * (read_len - len(r1))
    if len(r2) < read_len:
        r2 = r2 + "N" * (read_len - len(r2))
    return r1, r2

def simulate_illumina_reads(
    ligated_fragments: List[Dict],
    *,
    read_len: int = 150,
    target_coverage: float = 5.0,
    # error model
    substitution_error_rate: float = 0.002,
    homopolymer_slip_rate: float = 0.0005,
    rng: Optional[np.random.Generator] = None,
) -> List[Dict]:
    """
    Evenly sample fragments to achieve approximately `target_coverage`
    over the *insert* sequences (pre-adapter length in each fragment dict).

    Returns a list of dicts with: {'R1','R2','frag_idx','i7','i5','inline_bc'}.
    """
    if rng is None:
        rng = np.random.default_rng()

    if not ligated_fragments:
        return []

    # Total insert bp (not counting adapters)
    total_insert_bp = sum(f["length"] for f in ligated_fragments)
    if total_insert_bp == 0:
        return []

    # Pairs needed to reach coverage ≈ (total pairs * 2 * read_len) / total_insert_bp
    n_pairs = math.ceil(target_coverage * total_insert_bp / (2 * read_len))
    # Sample fragments roughly uniformly (could weight by length; DD by default is near-uniform after size-select)
    idxs = rng.integers(0, len(ligated_fragments), size=n_pairs)

    reads: List[Dict] = []
    for idx in idxs:
        frag = ligated_fragments[idx]
        r1, r2 = _extract_reads_from_oriented_fragment(frag, read_len, rng)

        # Apply errors (R1 then R2 independently)
        if substitution_error_rate > 0:
            r1 = _introduce_substitution_errors(r1, substitution_error_rate, rng)
            r2 = _introduce_substitution_errors(r2, substitution_error_rate, rng)
        if homopolymer_slip_rate > 0:
            r1 = _introduce_homopolymer_slips(r1, homopolymer_slip_rate, rng)
            r2 = _introduce_homopolymer_slips(r2, homopolymer_slip_rate, rng)

        reads.append({
            "R1": r1,
            "R2": r2,
            "frag_idx": idx,
            "i7": frag.get("i7", ""),
            "i5": frag.get("i5", ""),
            "inline_bc": frag.get("inline_bc", ""),
        })
    return reads



if __name__ == "__main__":

    # ---------------- convenience: example adapters and end-mapping ----------------

    # Common minimal stubs. You can drop in your exact kit sequences if you prefer.
    # Use '{index}' to splice in an i7/i5 index sequence if you want to carry it.
    ILLUMINA_P5_MIN = "AATGATACGGCGACCACCGAGATCTACAC{index}ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    ILLUMINA_P7_MIN = "CAAGCAGAAGACGGCATACGAGAT{index}GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    DEFAULT_ADAPTERS = Adapters(P5_i5=ILLUMINA_P5_MIN, P7_i7=ILLUMINA_P7_MIN)

    # ClaI (ATCGAT) cuts between ATCG | AT; set cut_offset=4
    claI  = Enzyme("ClaI",  "ATCGAT", 4)

    # BamHI (GGATCC) cuts between G | GATCC; set cut_offset=1
    bamHI = Enzyme("BamHI", "GGATCC", 1)

    # Simulate a 2 Mb scaffold, keep 300–500 bp ddRAD-style mixed-end fragments
    scaff, dd_frags = simulate_ddrad(
        genome_len=20_000_000,
        enzymes=[claI, bamHI],
        size_range=(150, 550),
        seed=42,
        require_mixed_ends=True
    )

    print(len(dd_frags), "fragments")
    print(dd_frags[0])


    # 2) Ligation: suppose P5 (R1) was put on the PstI side.
    #    If your earlier mapping says left_enzyme == 'PstI' for the “left” cut,
    #    then p5_is_left=True.
    ligated = ligate_adapters(
        dd_frags,
        p5_is_left=True,  # set based on your lab’s adapter scheme
        adapters=DEFAULT_ADAPTERS,
        inline_barcodes=["ACGTAC", "TGCATG", "GGTACC"],  # optional
        i7_indices=["ACGTACGT"],                         # optional
        i5_indices=["TGCATGCA"],                         # optional
    )

    # 3) Simulate reads to ~8× coverage over the inserts, 150 bp PE
    reads = simulate_illumina_reads(
        ligated,
        read_len=150,
        target_coverage=8.0,
        substitution_error_rate=0.003,    # ~0.3%
        homopolymer_slip_rate=0.0005,     # rare single-bp indels in runs >=3
    )

    print(f"Fragments: {len(dd_frags)}  Read pairs: {len(reads)}")
    print(reads[0]["R1"][:60], "...")
    print(reads[0]["R2"][:60], "...")

