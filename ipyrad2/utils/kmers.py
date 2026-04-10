#!/usr/bin/env python

"""Find likely restriction junctions from anchored read-start kmer counts."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple
import gzip

from loguru import logger

from .exceptions import IPyradError
from .parallel import run_with_pool


MIN_KMER_SIZE = 3
MAX_KMER_SIZE = 10
BRANCH_SEED_FRACTION = 10
BRANCH_EXTEND_NUMERATOR = 3
BRANCH_EXTEND_DENOMINATOR = 4
BRANCH_SPLIT_NUMERATOR = 1
BRANCH_SPLIT_DENOMINATOR = 2
MAX_ACCEPTED_MOTIFS = 3
TOP_LOGGED_KMERS = 5
SHIFTED_OFFSET_MARGIN_PERCENT = 10
IUPAC_BARCODE_MATCHES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "M": frozenset("AC"),
    "R": frozenset("AG"),
    "W": frozenset("AT"),
    "S": frozenset("CG"),
    "Y": frozenset("CT"),
    "K": frozenset("GT"),
    "N": frozenset("ACGT"),
}


@dataclass(frozen=True)
class InferredJunction:
    """Metadata for the primary inferred read-start junction."""

    sequence: str
    offset: int
    k: int
    winner_count: int
    runner_up_count: int
    candidate_offsets: Tuple[int, ...]

    @property
    def trim_length(self) -> int:
        """Number of leading bases to remove before the insert starts."""
        return self.offset + len(self.sequence)


@dataclass(frozen=True)
class InferredJunctionSet:
    """Metadata for one read-end junction inference result."""

    motifs: Tuple[str, ...]
    motif_counts: Tuple[int, ...]
    offset: int
    total_support: int
    runner_up_offset_support: int
    candidate_offsets: Tuple[int, ...]
    position_mode: str = "offset"
    sampled_reads: int = 0
    accepted_reads: int = 0
    skipped_no_match_reads: int = 0
    skipped_ambiguous_reads: int = 0
    boundary_supports: Tuple[Tuple[int, int, int], ...] = ()

    @property
    def trim_length(self) -> int:
        """Longest leading trim implied by the accepted motifs."""
        if not self.motifs:
            return self.offset
        return self.offset + max(len(i) for i in self.motifs)

    @property
    def primary_motif(self) -> str:
        """Return the strongest accepted motif."""
        return self.motifs[0] if self.motifs else ""

    @property
    def primary_count(self) -> int:
        """Return the support count for the strongest accepted motif."""
        return self.motif_counts[0] if self.motif_counts else 0

    @property
    def motif_support_fractions(self) -> Tuple[float, ...]:
        """Return motif support fractions within the accepted motif set."""
        if not self.total_support:
            return tuple(0.0 for _ in self.motifs)
        return tuple(count / self.total_support for count in self.motif_counts)

    @property
    def position_summary(self) -> str:
        """Return a short human-readable description of the inferred positions."""
        if self.position_mode == "barcode_boundary" and self.boundary_supports:
            return ", ".join(
                f"{length}+{slack}:{count}"
                for length, slack, count in self.boundary_supports
            )
        return f"offset {self.offset}"

    def as_primary(self) -> InferredJunction:
        """Convert to the legacy single-junction metadata shape."""
        runner_up_count = 0
        if len(self.motif_counts) > 1:
            runner_up_count = self.motif_counts[1]
        elif self.runner_up_offset_support:
            runner_up_count = self.runner_up_offset_support
        return InferredJunction(
            sequence=self.primary_motif,
            offset=self.offset,
            k=len(self.primary_motif),
            winner_count=self.primary_count,
            runner_up_count=runner_up_count,
            candidate_offsets=self.candidate_offsets,
        )


@dataclass(frozen=True)
class _OffsetInference:
    """Internal representation of one offset's retained motif set."""

    offset: int
    motifs: Tuple[str, ...]
    motif_counts: Tuple[int, ...]
    total_support: int


def _offset_candidate_sort_key(candidate: _OffsetInference) -> Tuple[int, int, int]:
    """Return a stable sort key for offset candidates."""
    return (-candidate.total_support, len(candidate.motifs), candidate.offset)


def _primary_support(candidate: _OffsetInference) -> int:
    """Return the strongest retained motif support for one offset."""
    return candidate.motif_counts[0] if candidate.motif_counts else 0


def _meets_shift_margin(shifted_support: int, base_support: int) -> bool:
    """Return True when offset=1 beats offset=0 by the required margin."""
    if base_support <= 0:
        return shifted_support > 0
    return shifted_support * 100 >= base_support * (100 + SHIFTED_OFFSET_MARGIN_PERCENT)


def _select_best_offset_candidate(
    candidates: Sequence[_OffsetInference],
    candidate_offsets: Tuple[int, ...],
) -> _OffsetInference:
    """Choose the best candidate offset, requiring strong evidence for +1 shifts."""
    ordered = sorted(candidates, key=_offset_candidate_sort_key)
    best = ordered[0]
    if tuple(candidate_offsets) != (0, 1) or best.offset != 1:
        return best

    candidates_by_offset = {candidate.offset: candidate for candidate in ordered}
    offset0 = candidates_by_offset.get(0)
    offset1 = candidates_by_offset.get(1)
    if offset0 is None or offset1 is None:
        return best
    if len(offset0.motifs) > MAX_ACCEPTED_MOTIFS and len(offset1.motifs) <= MAX_ACCEPTED_MOTIFS:
        logger.debug(
            "generic kmer inference accepted shifted offset=1 because offset=0 retained {} strong motifs, exceeding the supported maximum {}",
            len(offset0.motifs),
            MAX_ACCEPTED_MOTIFS,
        )
        return offset1

    offset0_primary = _primary_support(offset0)
    offset1_primary = _primary_support(offset1)
    total_margin_ok = _meets_shift_margin(offset1.total_support, offset0.total_support)
    primary_margin_ok = _meets_shift_margin(offset1_primary, offset0_primary)
    logger.trace(
        "generic kmer inference evaluating shifted offset=1 against offset=0: total_support {} vs {} primary_support {} vs {} required_margin={}%",
        offset1.total_support,
        offset0.total_support,
        offset1_primary,
        offset0_primary,
        SHIFTED_OFFSET_MARGIN_PERCENT,
    )
    if total_margin_ok and primary_margin_ok:
        logger.debug(
            "generic kmer inference accepted shifted offset=1 because total_support {} vs {} and primary_support {} vs {} met the required {}% margin",
            offset1.total_support,
            offset0.total_support,
            offset1_primary,
            offset0_primary,
            SHIFTED_OFFSET_MARGIN_PERCENT,
        )
        return offset1

    logger.debug(
        "generic kmer inference rejected shifted offset=1 because total_support {} vs {} and primary_support {} vs {} did not meet the required {}% margin; forcing offset=0",
        offset1.total_support,
        offset0.total_support,
        offset1_primary,
        offset0_primary,
        SHIFTED_OFFSET_MARGIN_PERCENT,
    )
    return offset0


@dataclass(frozen=True)
class _BoundaryCandidateInference:
    """One inferred motif family for one exact barcode-boundary class."""

    barcode_pattern: str
    barcode_length: int
    slack: int
    support_reads: int
    motifs: Tuple[str, ...]
    motif_counts: Tuple[int, ...]
    total_support: int


@dataclass(frozen=True)
class _BarcodeKmerCounts:
    """Per-file barcode-aware kmer counts and diagnostics."""

    counts: Dict[Tuple[int, int], Counter]
    candidate_counts: Dict[Tuple[str, int, int, int], Counter]
    boundary_counts: Dict[Tuple[int, int, int], Counter]
    boundary_supports: Tuple[Tuple[int, int, int], ...]
    candidate_supports: Tuple[Tuple[str, int, int, int], ...]
    sampled_reads: int
    accepted_reads: int
    skipped_no_match_reads: int
    skipped_ambiguous_reads: int
    reads_with_multiple_boundary_matches: int = 0


def is_informative_motif(seq: str | bytes) -> bool:
    """Return True for motifs that are not monomorphic after dropping Ns."""
    if isinstance(seq, bytes):
        seq = seq.decode()
    informative = {base for base in seq.upper() if base != "N"}
    return len(informative) >= 2


def validate_named_motif(
    seq: str | None,
    label: str,
    *,
    allow_empty: bool = False,
    allowed_chars: Iterable[str] | None = None,
) -> str:
    """Validate one user-provided motif."""
    seq = (seq or "").strip().upper()
    if not seq:
        if allow_empty:
            return ""
        raise IPyradError(f"{label} cannot be empty.")
    if allowed_chars is not None:
        allowed_chars = set(allowed_chars)
        if any(base not in allowed_chars for base in seq):
            raise IPyradError(f"{label} contains unexpected characters.")
    if not is_informative_motif(seq):
        raise IPyradError(f"{label} cannot be monomorphic or all N.")
    return seq


def validate_named_motif_list(
    seqs: str | Sequence[str] | None,
    label: str,
    *,
    allow_empty: bool = False,
    allowed_chars: Iterable[str] | None = None,
) -> Tuple[str, ...]:
    """Validate one or more comma-delimited motifs and preserve user order."""
    if seqs is None:
        if allow_empty:
            return ()
        raise IPyradError(f"{label} cannot be empty.")

    if isinstance(seqs, str):
        parts = [i.strip() for i in seqs.split(",")]
    else:
        parts = []
        for seq in seqs:
            parts.extend(part.strip() for part in str(seq).split(","))

    if not parts or all(not part for part in parts):
        if allow_empty:
            return ()
        raise IPyradError(f"{label} cannot be empty.")

    motifs: list[str] = []
    seen = set()
    for idx, part in enumerate(parts, start=1):
        if not part:
            raise IPyradError(f"{label} contains an empty motif entry.")
        motif = validate_named_motif(
            part,
            f"{label} motif {idx}",
            allow_empty=False,
            allowed_chars=allowed_chars,
        )
        if motif not in seen:
            motifs.append(motif)
            seen.add(motif)
    return tuple(motifs)


def _normalize_candidate_offsets(candidate_offsets: Sequence[int] | None) -> Tuple[int, ...]:
    """Return sorted unique non-negative candidate offsets."""
    if candidate_offsets is None:
        return (0,)
    offsets = tuple(sorted({int(i) for i in candidate_offsets if int(i) >= 0}))
    if not offsets:
        raise IPyradError("No candidate offsets were provided for kmer inference.")
    return offsets


def _normalize_barcodes_by_length(
    barcodes_by_length: Dict[int, Sequence[str]],
) -> Dict[int, Tuple[str, ...]]:
    """Return barcode patterns grouped by length in deterministic order."""
    normalized = {}
    for length, barcodes in sorted(barcodes_by_length.items()):
        values = tuple(sorted({str(barcode).upper() for barcode in barcodes if barcode}))
        if values:
            normalized[int(length)] = values
    if not normalized:
        raise IPyradError("No barcode patterns were provided for barcode-aware kmer inference.")
    return normalized


def _barcode_matches_pattern(seq: bytes, pattern: str) -> bool:
    """Return True when a read prefix matches a barcode pattern with zero mismatches."""
    if len(seq) != len(pattern):
        return False
    for base, expected in zip(seq.decode(), pattern):
        allowed = IUPAC_BARCODE_MATCHES.get(expected, frozenset(expected))
        if base not in allowed:
            return False
    return True


def _matching_barcode_candidates(
    read: bytes,
    barcodes_by_length: Dict[int, Tuple[str, ...]],
    max_slack: int = 1,
) -> List[Tuple[str, int, int]]:
    """Return exact barcode/slack boundaries that match this read exactly."""
    matches: list[Tuple[str, int, int]] = []
    for barcode_length, patterns in barcodes_by_length.items():
        for slack in range(max_slack + 1):
            end = slack + barcode_length
            if end > len(read):
                continue
            prefix = read[slack:end]
            for pattern in patterns:
                if _barcode_matches_pattern(prefix, pattern):
                    matches.append((pattern, barcode_length, slack))
    return list(dict.fromkeys(matches))


def _validate_barcode_boundary_slack(max_slack: int) -> int:
    """Return a normalized barcode-boundary slack value."""
    max_slack = int(max_slack)
    if max_slack not in (0, 1):
        raise IPyradError("barcode boundary slack must be 0 or 1.")
    return max_slack


def iter_reads(fastq: Path, max_len: int, max_reads: int) -> Iterator[bytes]:
    """Yield FASTQ sequences as uppercase bytes."""
    fastq = Path(fastq)
    open_func = gzip.open if fastq.suffix == ".gz" else open
    with open_func(fastq, "rb") as inline:
        quart = zip(inline, inline, inline, inline)
        for _, q in zip(range(max_reads), quart):
            yield q[1].strip().upper()[:max_len]


def get_kmer_counts(
    fastq: Path,
    candidate_offsets: Sequence[int],
    max_len: int,
    max_reads: int,
) -> Dict[Tuple[int, int], Counter]:
    """Return anchored prefix kmer counts for one FASTQ."""
    candidate_offsets = _normalize_candidate_offsets(candidate_offsets)
    counts = {
        (offset, kmer_size): Counter()
        for offset in candidate_offsets
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)
    }
    needed_len = max(max_len, max(candidate_offsets, default=0) + MAX_KMER_SIZE)
    for read in iter_reads(fastq, needed_len, max_reads):
        for offset in candidate_offsets:
            for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1):
                end = offset + kmer_size
                if end <= len(read):
                    counts[(offset, kmer_size)][read[offset:end]] += 1
    return counts


def _sorted_valid_motifs(counter: Counter) -> List[Tuple[str, int]]:
    """Return informative motifs from a Counter in deterministic support order."""
    motifs = []
    for seq, count in counter.items():
        motif = seq.decode()
        if is_informative_motif(motif):
            motifs.append((motif, count))
    motifs.sort(key=lambda item: (-item[1], item[0]))
    return motifs


def _top_kmers(counter: Counter, limit: int = TOP_LOGGED_KMERS) -> Tuple[Tuple[str, int], ...]:
    """Return the top informative kmers from a counter."""
    return tuple(_sorted_valid_motifs(counter)[:limit])


def _format_motif_list(motifs: Sequence[str]) -> str:
    """Return motifs in a concise human-readable list."""
    return "[" + ", ".join(motifs) + "]" if motifs else "[]"


def _format_boundary_positions(
    boundary_supports: Sequence[Tuple[int, int, int]],
) -> str:
    """Return concise barcode-boundary positions without support counts."""
    return ", ".join(
        f"{barcode_length}+{slack}"
        for barcode_length, slack, _count in boundary_supports
    )


def _format_cutsite_motif_inference_label(label: str) -> str:
    """Return a short user-facing label for cutsite motif inference."""
    return label or "cutsite motif inference"


def _log_offset_top_kmers(
    counts: Dict[Tuple[int, int], Counter],
    candidate_offsets: Tuple[int, ...],
    label: str,
) -> None:
    """Emit TRACE logs of top kmers for anchored offset candidates."""
    for offset in candidate_offsets:
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1):
            top = _top_kmers(counts[(offset, kmer_size)])
            if top:
                logger.trace(
                    "{} offset={} k={} top_kmers={}",
                    label,
                    offset,
                    kmer_size,
                    top,
                )


def _log_boundary_top_kmers(
    boundary_counts: Dict[Tuple[int, int, int], Counter],
    boundary_supports: Tuple[Tuple[int, int, int], ...],
    label: str,
) -> None:
    """Emit TRACE logs of top kmers for barcode-boundary candidates."""
    for barcode_length, slack, matched_reads in boundary_supports:
        if not matched_reads:
            continue
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1):
            top = _top_kmers(boundary_counts[(barcode_length, slack, kmer_size)])
            if top:
                logger.trace(
                    "{} boundary={}+{} matched_reads={} k={} top_kmers={}",
                    label,
                    barcode_length,
                    slack,
                    matched_reads,
                    kmer_size,
                    top,
                )


def _dedupe_active_branches(active: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Keep the strongest copy of each active branch sequence."""
    best: Dict[str, Tuple[int, int]] = {}
    for seq, count, kmer_size in active:
        current = best.get(seq)
        if current is None or count > current[0]:
            best[seq] = (count, kmer_size)
    ordered = [
        (seq, count, kmer_size)
        for seq, (count, kmer_size) in best.items()
    ]
    ordered.sort(key=lambda item: (-item[1], -len(item[0]), item[0]))
    return ordered


def _significant_children(
    counts: Dict[Tuple[int, int], Counter],
    offset: int,
    seq: str,
    next_kmer_size: int,
    top_seed_count: int,
) -> List[Tuple[str, int]]:
    """Return significant child motifs for one prefix at the next kmer size."""
    return [
        (child_seq, child_count)
        for child_seq, child_count in _sorted_valid_motifs(counts[(offset, next_kmer_size)])
        if child_seq.startswith(seq) and child_count * BRANCH_SEED_FRACTION >= top_seed_count
    ]


def _branch_terminals(
    counts: Dict[Tuple[int, int], Counter],
    offset: int,
    seq: str,
    count: int,
    kmer_size: int,
    top_seed_count: int,
    context: str,
) -> Tuple[Tuple[str, int], ...]:
    """Return retained terminal motifs reachable from one active branch."""
    if kmer_size >= MAX_KMER_SIZE:
        logger.trace(
            "{} retained motif {} at k={} because max kmer size was reached.",
            context,
            seq,
            kmer_size,
        )
        return ((seq, count),)

    children = _significant_children(counts, offset, seq, kmer_size + 1, top_seed_count)
    if not children:
        logger.trace(
            "{} retained motif {} at k={} because no significant children survived at k={}.",
            context,
            seq,
            kmer_size,
            kmer_size + 1,
        )
        return ((seq, count),)

    if len(children) == 1:
        child_seq, child_count = children[0]
        if child_count * BRANCH_EXTEND_DENOMINATOR < count * BRANCH_EXTEND_NUMERATOR:
            logger.trace(
                "{} retained motif {} at k={} because child {} support {} dropped below {:.2f}x parent support {}.",
                context,
                seq,
                kmer_size,
                child_seq,
                child_count,
                BRANCH_EXTEND_NUMERATOR / BRANCH_EXTEND_DENOMINATOR,
                count,
            )
            return ((seq, count),)
        return _branch_terminals(
            counts,
            offset,
            child_seq,
            child_count,
            kmer_size + 1,
            top_seed_count,
            context,
        )

    best_child_count = children[0][1]
    if (
        len(children) > MAX_ACCEPTED_MOTIFS or
        best_child_count * BRANCH_SPLIT_DENOMINATOR < count * BRANCH_SPLIT_NUMERATOR
    ):
        if len(children) > MAX_ACCEPTED_MOTIFS:
            logger.trace(
                "{} retained motif {} at k={} because {} significant children survived, exceeding max accepted motifs {}.",
                context,
                seq,
                kmer_size,
                len(children),
                MAX_ACCEPTED_MOTIFS,
            )
        else:
            logger.trace(
                "{} retained motif {} at k={} because child branches were diffuse; best child support {} was below {:.2f}x parent support {}.",
                context,
                seq,
                kmer_size,
                best_child_count,
                BRANCH_SPLIT_NUMERATOR / BRANCH_SPLIT_DENOMINATOR,
                count,
            )
        return ((seq, count),)

    child_terminals = []
    for child_seq, child_count in children:
        terminals = _branch_terminals(
            counts,
            offset,
            child_seq,
            child_count,
            kmer_size + 1,
            top_seed_count,
            context,
        )
        child_terminals.append((child_seq, child_count, terminals))

    retained: list[Tuple[str, int]] = []
    for child_seq, child_count, terminals in child_terminals:
        if any(len(term_seq) > len(child_seq) for term_seq, _ in terminals):
            retained.extend(terminals)
        else:
            retained.append((child_seq, child_count))
    return tuple(retained)


def _terminal_motifs_for_offset(
    counts: Dict[Tuple[int, int], Counter],
    offset: int,
    context: str | None = None,
) -> Tuple[Tuple[str, int], ...]:
    """Infer one or more retained motifs for one candidate offset."""
    context = context or f"offset={offset}"
    seeds = _sorted_valid_motifs(counts[(offset, MIN_KMER_SIZE)])
    if not seeds:
        return ()

    top_seed_count = seeds[0][1]
    terminals: list[Tuple[str, int]] = []
    for seq, count in seeds:
        if count * BRANCH_SEED_FRACTION >= top_seed_count:
            terminals.extend(
                _branch_terminals(
                    counts,
                    offset,
                    seq,
                    count,
                    MIN_KMER_SIZE,
                    top_seed_count,
                    context,
                )
            )

    best_by_seq: Dict[str, int] = {}
    for seq, count in terminals:
        best_by_seq[seq] = max(count, best_by_seq.get(seq, 0))

    shadowed = {
        seq for seq in best_by_seq
        if any(other != seq and other.startswith(seq) for other in best_by_seq)
    }
    retained = [
        (seq, count)
        for seq, count in best_by_seq.items()
        if seq not in shadowed
    ]
    retained.sort(key=lambda item: (-item[1], -len(item[0]), item[0]))
    if not retained:
        return ()

    strongest = retained[0][1]
    retained = tuple(
        item for item in retained
        if item[1] * BRANCH_SEED_FRACTION >= strongest
    )
    logger.trace("{} retained_terminal_motifs={}", context, retained)
    return retained


def _offset_inference(
    counts: Dict[Tuple[int, int], Counter],
    offset: int,
) -> _OffsetInference | None:
    """Build one candidate motif set for one offset."""
    retained = _terminal_motifs_for_offset(counts, offset, context=f"offset={offset}")
    if not retained:
        return None
    motifs = tuple(seq for seq, _ in retained)
    motif_counts = tuple(count for _, count in retained)
    return _OffsetInference(
        offset=offset,
        motifs=motifs,
        motif_counts=motif_counts,
        total_support=sum(motif_counts),
    )


def _select_best_inference_set(
    counts: Dict[Tuple[int, int], Counter],
    candidate_offsets: Tuple[int, ...],
) -> InferredJunctionSet:
    """Choose the best anchored junction set from merged counts."""
    candidates = [
        candidate
        for candidate in (
            _offset_inference(counts, offset)
            for offset in candidate_offsets
        )
        if candidate is not None
    ]
    if not candidates:
        raise IPyradError(
            "kmer analysis found only invalid or low-information motifs. "
            "Disable auto-infer-re-overhangs and set manually."
        )

    best = _select_best_offset_candidate(candidates, candidate_offsets)
    runner_ups = sorted(
        (candidate for candidate in candidates if candidate.offset != best.offset),
        key=_offset_candidate_sort_key,
    )
    runner_up_offset_support = runner_ups[0].total_support if runner_ups else 0
    if len(best.motifs) > MAX_ACCEPTED_MOTIFS:
        raise IPyradError(
            "kmer analysis found more than 3 strong motifs at one read end. "
            "Increase max_reads_kmer/--max-reads-kmer or set manually."
        )
    return InferredJunctionSet(
        motifs=best.motifs,
        motif_counts=best.motif_counts,
        offset=best.offset,
        total_support=best.total_support,
        runner_up_offset_support=runner_up_offset_support,
        candidate_offsets=candidate_offsets,
    )


def _merge_kmer_counts(
    kcounts: Dict[int, Dict[Tuple[int, int], Counter]],
    candidate_offsets: Tuple[int, ...],
) -> Dict[Tuple[int, int], Counter]:
    """Merge per-file kmer counters into one offset/kmer table."""
    merged_counts = {
        (offset, kmer_size): Counter()
        for offset in candidate_offsets
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)
    }
    for file_counts in kcounts.values():
        for key, counter in file_counts.items():
            merged_counts[key].update(counter)
    return merged_counts


def get_barcoded_kmer_counts(
    fastq: Path,
    barcodes_by_length: Dict[int, Tuple[str, ...]],
    max_len: int,
    max_reads: int,
    max_barcode_boundary_slack: int = 1,
) -> _BarcodeKmerCounts:
    """Return barcode-boundary-aware kmer counts for one FASTQ."""
    barcodes_by_length = _normalize_barcodes_by_length(barcodes_by_length)
    max_barcode_boundary_slack = _validate_barcode_boundary_slack(max_barcode_boundary_slack)
    slack_offsets = tuple(range(max_barcode_boundary_slack + 1))
    barcode_lengths = tuple(sorted(barcodes_by_length))
    counts = {(0, kmer_size): Counter() for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)}
    candidate_counts: Dict[Tuple[str, int, int, int], Counter] = {}
    boundary_counts = {
        (barcode_length, slack, kmer_size): Counter()
        for barcode_length in barcode_lengths
        for slack in slack_offsets
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)
    }
    boundary_supports = Counter()
    candidate_supports = Counter()
    sampled_reads = 0
    accepted_reads = 0
    skipped_no_match_reads = 0
    skipped_ambiguous_reads = 0
    needed_len = max(max_len, max(barcode_lengths, default=0) + max_barcode_boundary_slack + MAX_KMER_SIZE)
    for read in iter_reads(fastq, needed_len, max_reads):
        sampled_reads += 1
        matched_candidates = _matching_barcode_candidates(
            read,
            barcodes_by_length,
            max_slack=max_barcode_boundary_slack,
        )
        if not matched_candidates:
            skipped_no_match_reads += 1
            continue
        if len(matched_candidates) > 1:
            skipped_ambiguous_reads += 1
        accepted_reads += 1
        matched_boundary_positions = {(barcode_length, slack) for _, barcode_length, slack in matched_candidates}
        for barcode_length, slack in matched_boundary_positions:
            boundary_supports[(barcode_length, slack)] += 1
        for barcode_pattern, barcode_length, slack in matched_candidates:
            candidate_supports[(barcode_pattern, barcode_length, slack)] += 1
            start = barcode_length + slack
            for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1):
                end = start + kmer_size
                if end <= len(read):
                    kmer = read[start:end]
                    counts[(0, kmer_size)][kmer] += 1
                    boundary_counts[(barcode_length, slack, kmer_size)][kmer] += 1
                    key = (barcode_pattern, barcode_length, slack, kmer_size)
                    counter = candidate_counts.get(key)
                    if counter is None:
                        counter = Counter()
                        candidate_counts[key] = counter
                    counter[kmer] += 1
    return _BarcodeKmerCounts(
        counts=counts,
        candidate_counts=candidate_counts,
        boundary_counts=boundary_counts,
        boundary_supports=tuple(
            (barcode_length, slack, boundary_supports[(barcode_length, slack)])
            for barcode_length in barcode_lengths
            for slack in slack_offsets
            if boundary_supports[(barcode_length, slack)]
        ),
        candidate_supports=tuple(
            (barcode_pattern, barcode_length, slack, candidate_supports[(barcode_pattern, barcode_length, slack)])
            for barcode_length in barcode_lengths
            for slack in slack_offsets
            for barcode_pattern in barcodes_by_length[barcode_length]
            if candidate_supports[(barcode_pattern, barcode_length, slack)]
        ),
        sampled_reads=sampled_reads,
        accepted_reads=accepted_reads,
        skipped_no_match_reads=skipped_no_match_reads,
        skipped_ambiguous_reads=skipped_ambiguous_reads,
        reads_with_multiple_boundary_matches=skipped_ambiguous_reads,
    )


def _merge_barcoded_kmer_counts(
    kcounts: Dict[int, _BarcodeKmerCounts],
    barcode_lengths: Tuple[int, ...],
    max_barcode_boundary_slack: int = 1,
) -> _BarcodeKmerCounts:
    """Merge per-file barcode-aware kmer counters into one result."""
    max_barcode_boundary_slack = _validate_barcode_boundary_slack(max_barcode_boundary_slack)
    slack_offsets = tuple(range(max_barcode_boundary_slack + 1))
    merged_counts = {(0, kmer_size): Counter() for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)}
    merged_candidate_counts: Dict[Tuple[str, int, int, int], Counter] = {}
    merged_boundary_counts = {
        (barcode_length, slack, kmer_size): Counter()
        for barcode_length in barcode_lengths
        for slack in slack_offsets
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)
    }
    boundary_support_counter = Counter()
    candidate_support_counter = Counter()
    sampled_reads = 0
    accepted_reads = 0
    skipped_no_match_reads = 0
    skipped_ambiguous_reads = 0
    for result in kcounts.values():
        for key, counter in result.counts.items():
            merged_counts[key].update(counter)
        for key, counter in result.candidate_counts.items():
            merged_counter = merged_candidate_counts.get(key)
            if merged_counter is None:
                merged_counter = Counter()
                merged_candidate_counts[key] = merged_counter
            merged_counter.update(counter)
        for key, counter in result.boundary_counts.items():
            merged_boundary_counts[key].update(counter)
        for barcode_length, slack, count in result.boundary_supports:
            boundary_support_counter[(barcode_length, slack)] += count
        for barcode_pattern, barcode_length, slack, count in result.candidate_supports:
            candidate_support_counter[(barcode_pattern, barcode_length, slack)] += count
        sampled_reads += result.sampled_reads
        accepted_reads += result.accepted_reads
        skipped_no_match_reads += result.skipped_no_match_reads
        skipped_ambiguous_reads += result.skipped_ambiguous_reads
    candidate_order = sorted(candidate_support_counter)
    return _BarcodeKmerCounts(
        counts=merged_counts,
        candidate_counts=merged_candidate_counts,
        boundary_counts=merged_boundary_counts,
        boundary_supports=tuple(
            (barcode_length, slack, boundary_support_counter[(barcode_length, slack)])
            for barcode_length in barcode_lengths
            for slack in slack_offsets
            if boundary_support_counter[(barcode_length, slack)]
        ),
        candidate_supports=tuple(
            (
                barcode_pattern,
                barcode_length,
                slack,
                candidate_support_counter[(barcode_pattern, barcode_length, slack)],
            )
            for barcode_pattern, barcode_length, slack in candidate_order
            if candidate_support_counter[(barcode_pattern, barcode_length, slack)]
        ),
        sampled_reads=sampled_reads,
        accepted_reads=accepted_reads,
        skipped_no_match_reads=skipped_no_match_reads,
        skipped_ambiguous_reads=skipped_ambiguous_reads,
        reads_with_multiple_boundary_matches=skipped_ambiguous_reads,
    )


def _candidate_inference(
    merged: _BarcodeKmerCounts,
    barcode_pattern: str,
    barcode_length: int,
    slack: int,
    support_reads: int,
) -> _BoundaryCandidateInference | None:
    """Infer one motif family for one exact barcode-boundary class."""
    counts = {
        (0, kmer_size): merged.candidate_counts.get(
            (barcode_pattern, barcode_length, slack, kmer_size),
            Counter(),
        )
        for kmer_size in range(MIN_KMER_SIZE, MAX_KMER_SIZE + 1)
    }
    try:
        inferred = _select_best_inference_set(counts, (0,))
    except IPyradError:
        return None
    return _BoundaryCandidateInference(
        barcode_pattern=barcode_pattern,
        barcode_length=barcode_length,
        slack=slack,
        support_reads=support_reads,
        motifs=inferred.motifs,
        motif_counts=inferred.motif_counts,
        total_support=inferred.total_support,
    )


def _select_best_barcode_boundary_family(
    merged: _BarcodeKmerCounts,
    label: str,
) -> InferredJunctionSet:
    """Choose the strongest motif family across exact barcode-boundary classes."""
    candidate_inferences = [
        inference
        for inference in (
            _candidate_inference(merged, barcode_pattern, barcode_length, slack, support_reads)
            for barcode_pattern, barcode_length, slack, support_reads in merged.candidate_supports
        )
        if inference is not None
    ]
    if not candidate_inferences:
        if merged.reads_with_multiple_boundary_matches:
            raise IPyradError(
                "kmer analysis found only ambiguous barcode-boundary matches while searching for restriction motifs. "
                "This usually means a restriction motif appears inside one barcode and also creates another exact barcode boundary. "
                "Set cutsite motifs manually or adjust the barcode design."
            )
        raise IPyradError(
            "kmer analysis found only invalid or low-information motifs after barcode boundary matching. "
            "Set cutsite motifs manually."
        )

    grouped = {}
    for inference in candidate_inferences:
        state = grouped.setdefault(
            inference.motifs,
            {
                "total_support": 0,
                "motif_counts": [0 for _ in inference.motif_counts],
                "boundary_supports": Counter(),
                "candidate_count": 0,
            },
        )
        state["total_support"] += inference.total_support
        state["candidate_count"] += 1
        for idx, count in enumerate(inference.motif_counts):
            state["motif_counts"][idx] += count
        state["boundary_supports"][(inference.barcode_length, inference.slack)] += inference.support_reads

    ranked = sorted(
        grouped.items(),
        key=lambda item: (
            -item[1]["total_support"],
            -item[1]["candidate_count"],
            len(item[0]),
            item[0],
        ),
    )
    winner_motifs, winner_state = ranked[0]
    runner_up_support = ranked[1][1]["total_support"] if len(ranked) > 1 else 0
    boundary_supports = tuple(
        (barcode_length, slack, winner_state["boundary_supports"][(barcode_length, slack)])
        for barcode_length, slack in sorted(winner_state["boundary_supports"])
    )
    logger.debug(
        "{}: sampled {} reads, matched barcode boundaries in {}, evaluated {} boundary classes",
        _format_cutsite_motif_inference_label(label),
        merged.sampled_reads,
        merged.accepted_reads,
        len(candidate_inferences),
    )
    logger.trace(
        "{} detailed selection: sampled_reads={} boundary_matched_reads={} "
        "reads_with_multiple_boundary_matches={} candidate_boundaries={} "
        "selected_motifs={} total_support={} boundary_supports={}",
        label,
        merged.sampled_reads,
        merged.accepted_reads,
        merged.reads_with_multiple_boundary_matches,
        len(candidate_inferences),
        winner_motifs,
        winner_state["total_support"],
        boundary_supports,
    )
    return InferredJunctionSet(
        motifs=winner_motifs,
        motif_counts=tuple(winner_state["motif_counts"]),
        offset=0,
        total_support=winner_state["total_support"],
        runner_up_offset_support=runner_up_support,
        candidate_offsets=(0,),
        position_mode="barcode_boundary",
        sampled_reads=merged.sampled_reads,
        accepted_reads=merged.accepted_reads,
        skipped_no_match_reads=merged.skipped_no_match_reads,
        skipped_ambiguous_reads=merged.skipped_ambiguous_reads,
        boundary_supports=boundary_supports,
    )


def get_overhangs_from_kmers(
    fastqs: List[Path],
    max_len: int,
    max_reads: int,
    workers: int,
    log_level: str,
    candidate_offsets: Sequence[int] | None = None,
    *,
    label: str = "cutsite motif inference",
) -> InferredJunctionSet:
    """Infer one or more read-start junction motifs from anchored prefix counts."""
    if not fastqs:
        raise IPyradError("No FASTQ files were provided for kmer inference.")

    candidate_offsets = _normalize_candidate_offsets(candidate_offsets)
    max_reads_per_file = max(1, ceil(max_reads / len(fastqs)))
    jobs = {}
    for idx, fastq in enumerate(fastqs):
        jobs[idx] = (
            get_kmer_counts,
            dict(
                fastq=fastq,
                candidate_offsets=candidate_offsets,
                max_len=max_len,
                max_reads=max_reads_per_file,
            ),
        )

    kcounts = run_with_pool(jobs, log_level, workers, msg="Counting kmers")
    merged_counts = _merge_kmer_counts(kcounts, candidate_offsets)
    _log_offset_top_kmers(merged_counts, candidate_offsets, label)
    best = _select_best_inference_set(merged_counts, candidate_offsets)
    logger.debug(
        "{} summary: motifs {}, offset {}, support {}, runner-up-offset-support {}",
        label,
        _format_motif_list(best.motifs),
        best.offset,
        best.total_support,
        best.runner_up_offset_support,
    )
    logger.trace(
        "{} detailed summary: selected_motifs={} offset={} total_support={} runner_up_offset_support={} top5={} top6={}",
        label,
        best.motifs,
        best.offset,
        best.total_support,
        best.runner_up_offset_support,
        _top_kmers(merged_counts[(best.offset, 5)]),
        _top_kmers(merged_counts[(best.offset, 6)]),
    )
    return best


def get_overhangs_from_barcoded_reads(
    fastqs: List[Path],
    barcodes_by_length: Dict[int, Sequence[str]],
    max_len: int,
    max_reads: int,
    workers: int,
    log_level: str,
    *,
    label: str = "demux",
    max_barcode_boundary_slack: int = 1,
) -> InferredJunctionSet:
    """Infer read-start motifs after exact barcode-prefix matching."""
    if not fastqs:
        raise IPyradError("No FASTQ files were provided for barcode-aware kmer inference.")

    max_barcode_boundary_slack = _validate_barcode_boundary_slack(max_barcode_boundary_slack)
    barcodes_by_length = _normalize_barcodes_by_length(barcodes_by_length)
    barcode_lengths = tuple(sorted(barcodes_by_length))
    max_reads_per_file = max(1, ceil(max_reads / len(fastqs)))
    jobs = {}
    for idx, fastq in enumerate(fastqs):
        jobs[idx] = (
            get_barcoded_kmer_counts,
            dict(
                fastq=fastq,
                barcodes_by_length=barcodes_by_length,
                max_len=max_len,
                max_reads=max_reads_per_file,
                max_barcode_boundary_slack=max_barcode_boundary_slack,
            ),
        )

    raw_counts = run_with_pool(jobs, log_level, workers, msg="Counting kmers")
    merged = _merge_barcoded_kmer_counts(
        raw_counts,
        barcode_lengths,
        max_barcode_boundary_slack=max_barcode_boundary_slack,
    )
    if not merged.accepted_reads:
        if merged.reads_with_multiple_boundary_matches:
            raise IPyradError(
                "kmer analysis found only multiply-matching barcode boundaries and no uniquely usable barcode-boundary reads. "
                "This usually means a restriction motif appears inside one barcode and creates another exact boundary. "
                "Set cutsite motifs manually or adjust the barcode design."
            )
        raise IPyradError(
            "kmer analysis found no reads with an exact barcode-prefix match at expected boundaries. "
            "Set cutsite motifs manually."
        )
    _log_boundary_top_kmers(merged.boundary_counts, merged.boundary_supports, label)
    best = _select_best_barcode_boundary_family(merged, label)
    logger.debug(
        "{} summary: motifs {}, support {}, no-barcode-match {}, multiple-boundary-match {}",
        _format_cutsite_motif_inference_label(label),
        _format_motif_list(best.motifs),
        best.total_support,
        best.skipped_no_match_reads,
        best.skipped_ambiguous_reads,
    )
    logger.trace(
        "{} detailed summary: sampled_reads={} accepted_reads={} skipped_no_match_reads={} "
        "reads_with_multiple_boundary_matches={} boundary_supports={} boundary_positions={} "
        "selected_motifs={} total_support={} top5={} top6={}",
        label,
        best.sampled_reads,
        best.accepted_reads,
        best.skipped_no_match_reads,
        best.skipped_ambiguous_reads,
        best.boundary_supports,
        _format_boundary_positions(best.boundary_supports),
        best.motifs,
        best.total_support,
        _top_kmers(merged.counts[(0, 5)]),
        _top_kmers(merged.counts[(0, 6)]),
    )
    return best


def get_overhang_from_kmers(
    fastqs: List[Path],
    max_len: int,
    max_reads: int,
    workers: int,
    log_level: str,
    candidate_offsets: Sequence[int] | None = None,
) -> InferredJunction:
    """Infer the primary read-start junction from anchored prefix kmer counts."""
    return get_overhangs_from_kmers(
        fastqs,
        max_len,
        max_reads,
        workers,
        log_level,
        candidate_offsets=candidate_offsets,
    ).as_primary()


if __name__ == "__main__":
    DIR = Path("/home/deren/Documents/ipyrad-tests/examples/Ama-PE-ddRAD/")
    R1s = list(DIR.glob("SLH_AL*_R2*"))
    x = get_overhangs_from_kmers(R1s, 18, 500_000, 10, "INFO")
    print(x)
