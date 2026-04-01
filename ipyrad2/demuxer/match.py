#!/usr/bin/env python

"""Barcode-matching classes and serial helpers for demux."""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, Iterator, List, Sequence, Tuple
import io
from pathlib import Path
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from loguru import logger
from ..utils.exceptions import IPyradError
from ..utils.parallel import ParallelJobError
from ..utils.parallel.pool import _ManagedProcessPool
from .sample_names import (
    final_output_sample_name,
    is_technical_replicate_name,
    technical_replicate_base_name,
)


DEFAULT_PIPELINE_BATCH_BYTES = 64 * 1024 * 1024
DEFAULT_WRITER_FLUSH_BYTES = 1 * 1024 * 1024
DEFAULT_PROGRESS_REPORT_READS = 10_000
DEFAULT_QUEUE_PUT_TIMEOUT = 0.25


@dataclass(frozen=True)
class DemuxRunConfig:
    """Static configuration shared across demux reader and writer workers."""

    barcodes_to_names: Dict[bytes, str]
    barcode_lengths1: Tuple[int, ...]
    barcode_lengths2: Tuple[int, ...]
    cuts1: List[bytes]
    cuts2: List[bytes]
    merge_technical_replicates: bool
    outdir: Path
    chunksize: int
    max_reads: int | None
    i7: bool
    log_level: str
    barcodes_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    barcode1_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    barcode2_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    barcode1_candidates_by_length: Dict[int, frozenset[bytes]] = field(default_factory=dict)
    barcode2_candidates_by_length: Dict[int, frozenset[bytes]] = field(default_factory=dict)
    barcode1_mismatch_by_barcode: Dict[bytes, int] = field(default_factory=dict)
    barcode2_mismatch_by_barcode: Dict[bytes, int] = field(default_factory=dict)
    pigz: bool = False
    batch_bytes: int = DEFAULT_PIPELINE_BATCH_BYTES
    writer_flush_bytes: int = DEFAULT_WRITER_FLUSH_BYTES
    queue_put_timeout: float = DEFAULT_QUEUE_PUT_TIMEOUT


@dataclass(frozen=True)
class BarcodeBoundaryCandidate:
    """One exact barcode boundary candidate found on a read end."""

    barcode: bytes
    slack: int
    trim_start: int
    mismatch_distance: int = 0

@dataclass
class BarMatching:
    """Base class for barcode matching.

    See subclasses which have different versions of the function
    `_iter_matched_barcode` to find barcode matches based on i7,
    combinatorial, or single inline barcodes. The subclasses all
    share the functions of this class, which includes iterating
    over the fastq(s), storing stats, and writing to tmp files.
    """
    fastqs: Tuple[Path, Path | None]
    """: A tuple with paired R1 and R2 fastq files."""
    barcodes_to_names: Dict[bytes, str]
    """: Dict matching barcodes to sample names."""
    barcode_lengths1: Tuple[int, ...]
    """: Expected barcode lengths for read1."""
    barcode_lengths2: Tuple[int, ...]
    """: Expected barcode lengths for read2."""
    cuts1: List[bytes]
    """: List of RE overhangs to match on R1."""
    cuts2: List[bytes]
    """: List of RE overhangs to match on R2."""
    merge_technical_replicates: bool
    """: ..."""
    outdir: Path
    """: ..."""
    log_level: str
    """: Log level to use for helper worker pools."""
    workers: int
    """: ..."""
    chunksize: int
    """: Number of reads to store in memory before writing to disk."""
    barcodes_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    """: Runtime barcode combinations mapped to all matching sample names."""
    barcode1_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    """: Runtime R1 barcode candidates mapped to all matching sample names."""
    barcode2_to_samples: Dict[bytes, Tuple[str, ...]] = field(default_factory=dict)
    """: Runtime R2 barcode candidates mapped to all matching sample names."""
    barcode1_candidates_by_length: Dict[int, frozenset[bytes]] = field(default_factory=dict)
    """: Runtime R1 barcode candidates grouped by length for exact boundary matching."""
    barcode2_candidates_by_length: Dict[int, frozenset[bytes]] = field(default_factory=dict)
    """: Runtime R2 barcode candidates grouped by length for exact boundary matching."""
    barcode1_mismatch_by_barcode: Dict[bytes, int] = field(default_factory=dict)
    """: Minimum mismatch distance for each acceptable R1 barcode candidate."""
    barcode2_mismatch_by_barcode: Dict[bytes, int] = field(default_factory=dict)
    """: Minimum mismatch distance for each acceptable R2 barcode candidate."""
    max_reads: int | None = int(1e20)
    """: Only sample this many reads from a file (mainly used in testing)."""
    progress_callback: Callable[[int, int], None] | None = None
    """: Optional callback receiving absolute (raw_reads, matched_reads)."""
    progress_interval_reads: int = DEFAULT_PROGRESS_REPORT_READS
    """: Minimum raw-read interval between progress callback updates."""

    # stats counters
    barcode_misses: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that don't match."""
    barcode_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record observed barcodes that match."""
    sample_hits: Dict[str, int] = field(default_factory=dict)
    """: Dict to record number of hits per sample."""
    barcode_boundary_ambiguities: Dict[bytes, int] = field(default_factory=dict)
    """: Dict to record ambiguous barcode-boundary candidate sets."""
    reads_seen: int = 0
    """: Absolute raw reads examined."""
    matched_seen: int = 0
    """: Absolute reads matched to a sample."""
    _last_progress_reads: int = 0
    """: Last raw-read count reported to progress_callback."""

    def __post_init__(self):
        self._format_check()

    def _ensure_runtime_barcode_maps(self) -> None:
        """Backfill end-specific barcode candidate maps from legacy combo maps when needed."""
        if self.barcode1_candidates_by_length:
            return
        combo_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        barcode1_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        barcode2_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        for barcode, sample_name in self.barcodes_to_names.items():
            combo_to_samples[barcode].add(sample_name)
            if b"_" in barcode:
                barcode1, barcode2 = barcode.split(b"_", 1)
                barcode1_to_samples[barcode1].add(sample_name)
                barcode2_to_samples[barcode2].add(sample_name)
            else:
                barcode1_to_samples[barcode].add(sample_name)
        self.barcodes_to_samples = {
            barcode: tuple(sorted(samples))
            for barcode, samples in combo_to_samples.items()
        }
        self.barcode1_to_samples = {
            barcode: tuple(sorted(samples))
            for barcode, samples in barcode1_to_samples.items()
        }
        self.barcode2_to_samples = {
            barcode: tuple(sorted(samples))
            for barcode, samples in barcode2_to_samples.items()
        }
        self.barcode1_candidates_by_length = {}
        for barcode in self.barcode1_to_samples:
            self.barcode1_candidates_by_length.setdefault(len(barcode), set()).add(barcode)
        self.barcode1_candidates_by_length = {
            length: frozenset(values)
            for length, values in self.barcode1_candidates_by_length.items()
        }
        self.barcode2_candidates_by_length = {}
        for barcode in self.barcode2_to_samples:
            self.barcode2_candidates_by_length.setdefault(len(barcode), set()).add(barcode)
        self.barcode2_candidates_by_length = {
            length: frozenset(values)
            for length, values in self.barcode2_candidates_by_length.items()
        }
        self.barcode1_mismatch_by_barcode = {barcode: 0 for barcode in self.barcode1_to_samples}
        self.barcode2_mismatch_by_barcode = {barcode: 0 for barcode in self.barcode2_to_samples}

    def _format_check(self):
        """Check that data is appropriate for selected format."""
        pass

    def _maybe_report_progress(self, force: bool = False) -> None:
        """Report absolute raw/matched counters to the optional callback."""
        if self.progress_callback is None:
            return
        if force or (self.reads_seen - self._last_progress_reads) >= self.progress_interval_reads:
            self.progress_callback(self.reads_seen, self.matched_seen)
            self._last_progress_reads = self.reads_seen

    def _iter_fastq_reads(
        self,
    ) -> Iterator[
        Tuple[
            Tuple[bytes, bytes, bytes, bytes],
            Tuple[bytes, bytes, bytes, bytes] | tuple[()],
        ]
    ]:
        """Yields fastq quartets of lines from fastqs (gzip OK)."""
        # create first read iterator for paired data
        opener = gzip.open if self.fastqs[0].suffix == ".gz" else io.open
        # ofile1 = opener(self.fastqs[0], 'rt', encoding="utf-8")
        ofile1 = opener(self.fastqs[0], 'rb')
        quart1 = zip(ofile1, ofile1, ofile1, ofile1)

        # create second read iterator for paired data
        if self.fastqs[1]:
            # ofile2 = opener(self.fastqs[1], 'rt', encoding="utf-8")
            ofile2 = opener(self.fastqs[1], 'rb')
            quart2 = zip(ofile2, ofile2, ofile2, ofile2)
        else:
            quart2 = iter(int, 1)

        # yield from iterators as 4 items as a time (fastq)
        ridx = 0
        for read1, read2 in zip(quart1, quart2):
            # stop if max_reads is reached
            ridx += 1
            if self.max_reads is not None and ridx > self.max_reads:
                return
            self.reads_seen += 1
            self._maybe_report_progress()
            yield read1, read2

    def _iter_matched_barcode(self):
        """SUBCLASSES REPLACE THIS FUNCTION."""
        raise NotImplementedError("See subclasses.")

    def _output_sample_name(self, sample_name: str) -> str:
        """Return the output sample name after optional replicate merging."""
        return final_output_sample_name(sample_name, self.merge_technical_replicates)

    def iter_output_records(self) -> Iterator[Tuple[str, bytes, bytes | None]]:
        """Yield matched FASTQ payloads ready for downstream writers."""
        for read1, read2, match in self._iter_matched_barcode():
            self.matched_seen += 1
            self._maybe_report_progress()
            yield (
                self._output_sample_name(match),
                b"".join(read1),
                b"".join(read2) if read2 else None,
            )

    def _iter_matched_chunks(
        self,
    ) -> Iterator[Tuple[Dict[str, List[bytes]], Dict[str, List[bytes]]]]:
        """Stores matched reads until N then writes to file."""
        read1s = {}
        read2s = {}
        nstored = 0

        # iterate over matched reads
        for read1, read2, match in self._iter_matched_barcode():
            # store r1 as 4-line string
            fastq1 = b"".join(read1)
            if match in read1s:
                read1s[match].append(fastq1)
            else:
                read1s[match] = [fastq1]

            # store r2 as 4-line string
            if read2:
                fastq2 = b"".join(read2)
                if match in read2s:
                    read2s[match].append(fastq2)
                else:
                    read2s[match] = [fastq2]

            # write to file when size is big enough and reset.
            nstored += 1
            if nstored > self.chunksize:
                yield read1s, read2s
                read1s = {}
                read2s = {}
                nstored = 0

        # write final chunk if data
        yield read1s, read2s

    def _build_write_jobs(
        self,
        read1s: Dict[str, List[bytes]],
        read2s: Dict[str, List[bytes]],
    ) -> Dict[str, Tuple[Callable[..., Any], Dict[str, object]]]:
        """Build per-file write jobs for one buffered chunk."""
        jobs = {}
        for name in read1s:
            jobs[f"{name}_R1"] = (
                write,
                dict(path=self.outdir / f"{name}_R1.fastq.gz", data=read1s[name]),
            )
            if read2s:
                jobs[f"{name}_R2"] = (
                    write,
                    dict(path=self.outdir / f"{name}_R2.fastq.gz", data=read2s[name]),
                )
        return jobs

    def run(self) -> None:
        """Multiprocessed writing is much faster, especially on HPC.

        Some overhead from i/o limitations, but most time here is spent
        on the string concatenation and gzip compression, which can
        happen in parallel on different engines.

        TODO
        ----
        Allow writers to use read1s or read2s dict as shared memory
        to prevent duplication of the Memory usage here on every CPU.
        Example here: https://stackoverflow.com/questions/65980183/processpoolexecutor-on-shared-dataset-and-multiple-arguments
        """
        pool = _ManagedProcessPool(log_level=self.log_level, max_workers=self.workers)
        close_wait = True
        try:
            total = 0
            for read1s, read2s in self._iter_matched_chunks():
                if not read1s:
                    continue
                nprocessed = min(self.chunksize, sum(len(i) for i in read1s.values()))
                total += nprocessed
                logger.info(
                    f"writing/compressing {nprocessed:.0f} matched reads (total={total:.0f})")

                # parallel workers cannot write to the same file so lists
                # assigned to technical replicates need to be grouped
                if self.merge_technical_replicates:
                    _merge_technical_replicate_chunks(read1s, read2s)

                jobs = self._build_write_jobs(read1s, read2s)
                for _key, _result in pool.iter_results(jobs.items(), max_inflight=self.workers):
                    pass
        except ParallelJobError:
            raise
        except KeyboardInterrupt:
            pool.abort(fast=True)
            close_wait = False
            raise
        finally:
            pool.close(wait=close_wait)


@dataclass
class BarMatchingI7(BarMatching):
    """Subclass of Barmatching that matches barcode in i7 header.

    Example 3RAD R1 file with i7+i5 tag in header
    ---------------------------------------------
    >>> # asterisk part is the i7 --->                  ********
    >>> @NB551405:60:H7T2GAFXY:4:21612:8472:20380 1:N:0:TATCGGTC+ACCAGGGA
    >>> ATCGGTATGCTGGAGGTGGTGGTGGTGGAGGTGGACGTTACAAGGGTTCTGGTGGTAGCCGATCAG...
    >>> +
    >>> EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE...
    """
    # check that i7's exist for this data
    def _format_check(self):
        read1, _ = next(self._iter_fastq_reads())
        barcode = read1[0].strip().rsplit(b":", 1)[-1].split(b"+")[0]
        if not barcode:
            raise IPyradError(
                "No i7 index exists in this data. Example read1:\n{read1}")

    def _iter_matched_barcode(self) -> Iterator[Tuple[str, str, str]]:
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():
            # pull barcode from header
            barcode = read1[0].strip().rsplit(b":", 1)[-1].split(b"+")[0]
            # look for match
            match = self.barcodes_to_names.get(barcode)

            # record stats and yield the reads if matched.
            if match:
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1
                yield read1, read2, match
            else:
                self.barcode_misses[barcode] = self.barcode_misses.get(barcode, 0) + 1


@dataclass
class BarMatchingSingleInline(BarMatching):
    """Subclass of Barmatching SE or PE data w/ inline barcodes only on R1.

    Example R1 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= cutsite motif.
    >>>
    >>> ********-----
    >>> @E00526:227:H53YNCCX2:8:1202:7710:23354 1:N:0:
    >>> CTGCAACTATCGGAGCGAATGAAAC........GACTCAACATAACGGGTCTGATCATTGAG
    >>> +
    >>> AA<FFJJJJJJJJJJJJJJJJJJJJ........JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    """
    maxlen1: int = 0
    """: Max len of the read1 inline barcodes."""

    def __post_init__(self):
        self._ensure_runtime_barcode_maps()
        self.maxlen1 = _match_window_length(
            self.barcode1_candidates_by_length,
            self.barcode_lengths1,
            self.cuts1,
        )

    def _iter_matched_barcode(self) -> Iterator[Tuple[str, str, str]]:
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        for read1, read2 in self._iter_fastq_reads():
            candidates = match_barcode_candidates(
                read1[1][:self.maxlen1],
                self.barcode1_candidates_by_length,
                self.cuts1,
                self.barcode1_mismatch_by_barcode,
            )
            sample_candidates: Dict[str, List[BarcodeBoundaryCandidate]] = {}
            for candidate in candidates:
                for sample_name in self.barcode1_to_samples.get(candidate.barcode, ()):
                    sample_candidates.setdefault(sample_name, []).append(candidate)

            if len(sample_candidates) == 1:
                match = next(iter(sample_candidates))
                chosen = _best_candidate(sample_candidates[match])
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[chosen.barcode] = self.barcode_hits.get(chosen.barcode, 0) + 1
                read1 = _trim_fastq_record(read1, chosen.trim_start)
                yield read1, read2, match
                continue

            if sample_candidates:
                _record_boundary_ambiguity(self.barcode_boundary_ambiguities, sample_candidates)
                continue

            _record_barcode_miss(self.barcode_misses, candidates)


@dataclass
class BarMatchingCombinatorialInline(BarMatching):
    """Subclass of Barmatching for combinatorial inline barcodes.

    Example R1 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= cutsite motif.
    >>>
    >>> ********-----
    >>> @E00526:227:H53YNCCX2:8:1202:7710:23354 1:N:0:
    >>> CTGCAACTATCGGAGCGAATGAAAC........GACTCAACATAACGGGTCTGATCATTGAG
    >>> +
    >>> AA<FFJJJJJJJJJJJJJJJJJJJJ........JJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

    Example R2 with inline barcodes
    -------------------------------
    >>> # '*'=inline barcode, '-'= cutsite motif.
    >>>
    >>> ********----
    >>> @E00526:227:H53YNCCX2:8:1202:7446:23354 2:N:0:CGAACTGT+ACAACAGT
    >>> ATGCTGTCGATCCCAACCACCACGC........TTTTTTTCTATCTCAACTATTTACAACAA
    >>> +
    >>> AAFFFJJJJJJJJJFJFJJJJJJ-F........AFJ<JFJJJJAJFFAA-F<A-AAF-AFFJ
    """
    maxlen1: int = 0
    """: Max len of the read1 inline barcodes + re."""
    maxlen2: int = 0
    """: Max len of the read2 inline barcodes + re."""

    def __post_init__(self):
        self._ensure_runtime_barcode_maps()
        self.maxlen1 = _match_window_length(
            self.barcode1_candidates_by_length,
            self.barcode_lengths1,
            self.cuts1,
        )
        self.maxlen2 = _match_window_length(
            self.barcode2_candidates_by_length,
            self.barcode_lengths2,
            self.cuts2,
        )

    def _iter_matched_barcode(self):
        """Find barcode in read and check for match.

        In i7 matching there is nothing to be trimmed from the reads.
        """
        # get a list of cutters and off-by-one's
        for read1, read2 in self._iter_fastq_reads():
            candidates_r1 = match_barcode_candidates(
                read1[1][:self.maxlen1],
                self.barcode1_candidates_by_length,
                self.cuts1,
                self.barcode1_mismatch_by_barcode,
            )
            candidates_r2 = match_barcode_candidates(
                read2[1][:self.maxlen2],
                self.barcode2_candidates_by_length,
                self.cuts2,
                self.barcode2_mismatch_by_barcode,
            )

            sample_candidates: Dict[str, List[Tuple[BarcodeBoundaryCandidate, BarcodeBoundaryCandidate, bytes]]] = {}
            for candidate_r1 in candidates_r1:
                for candidate_r2 in candidates_r2:
                    barcode = candidate_r1.barcode + b"_" + candidate_r2.barcode
                    for sample_name in self.barcodes_to_samples.get(barcode, ()):
                        sample_candidates.setdefault(sample_name, []).append(
                            (candidate_r1, candidate_r2, barcode)
                        )

            if len(sample_candidates) == 1:
                match = next(iter(sample_candidates))
                chosen_r1, chosen_r2, barcode = min(
                    sample_candidates[match],
                    key=lambda item: (
                        item[0].mismatch_distance + item[1].mismatch_distance,
                        item[0].slack + item[1].slack,
                        -(len(item[0].barcode) + len(item[1].barcode)),
                        item[2],
                    ),
                )
                self.sample_hits[match] = self.sample_hits.get(match, 0) + 1
                self.barcode_hits[barcode] = self.barcode_hits.get(barcode, 0) + 1
                read1 = _trim_fastq_record(read1, chosen_r1.trim_start)
                read2 = _trim_fastq_record(read2, chosen_r2.trim_start)
                yield read1, read2, match
                continue

            if sample_candidates:
                _record_boundary_ambiguity(
                    self.barcode_boundary_ambiguities,
                    {
                        sample_name: [pair[0] for pair in pairs] + [pair[1] for pair in pairs]
                        for sample_name, pairs in sample_candidates.items()
                    }
                )
                continue

            _record_combined_barcode_miss(self.barcode_misses, candidates_r1, candidates_r2)

def cut_matcher(
    read: bytes,
    barcode_lengths: Tuple[int, ...],
    cutters: List[bytes],
) -> Tuple[bytes, int] | None:
    """Return the matched barcode and trim start at valid barcode boundaries."""
    matches = _match_boundary_candidates(
        read,
        barcode_lengths,
        cutters,
        barcode_candidates_by_length=None,
    )
    if len(matches) == 1:
        candidate = matches[0]
        return candidate.barcode, candidate.trim_start
    if matches:
        min_slack = min(match.slack for match in matches)
        preferred = [match for match in matches if match.slack == min_slack]
        if len(preferred) == 1:
            candidate = preferred[0]
            return candidate.barcode, candidate.trim_start
    return None


def _match_boundary_candidates(
    read: bytes,
    barcode_lengths: Iterable[int],
    cutters: Sequence[bytes],
    barcode_candidates_by_length: Dict[int, frozenset[bytes]] | None,
    mismatch_by_barcode: Dict[bytes, int] | None = None,
) -> List[BarcodeBoundaryCandidate]:
    """Return deterministic barcode-boundary candidates supported by an immediate cutsite match."""
    mismatch_by_barcode = mismatch_by_barcode or {}
    matches: Dict[Tuple[bytes, int, int], BarcodeBoundaryCandidate] = {}
    for barcode_length in sorted(barcode_lengths):
        allowed_barcodes = None
        if barcode_candidates_by_length is not None:
            allowed_barcodes = barcode_candidates_by_length.get(barcode_length)
            if not allowed_barcodes:
                continue
        for slack in (0, 1):
            cut_start = slack + barcode_length
            if cut_start <= 0 or cut_start > len(read):
                continue
            barcode = read[slack:cut_start]
            if allowed_barcodes is not None and barcode not in allowed_barcodes:
                continue
            for cut in cutters:
                cut_end = cut_start + len(cut)
                if cut_end <= len(read) and read[cut_start:cut_end] == cut:
                    key = (barcode, slack, cut_start)
                    matches[key] = BarcodeBoundaryCandidate(
                        barcode=barcode,
                        slack=slack,
                        trim_start=cut_start,
                        mismatch_distance=mismatch_by_barcode.get(barcode, 0),
                    )
                    break
    return sorted(
        matches.values(),
        key=lambda item: (item.mismatch_distance, item.slack, -len(item.barcode), item.barcode),
    )


def match_barcode_candidates(
    read: bytes,
    barcode_candidates_by_length: Dict[int, frozenset[bytes]],
    cutters: List[bytes],
    mismatch_by_barcode: Dict[bytes, int] | None = None,
) -> List[BarcodeBoundaryCandidate]:
    """Return exact barcode-boundary candidates supported by an immediate cutsite match."""
    return _match_boundary_candidates(
        read,
        barcode_candidates_by_length,
        cutters,
        barcode_candidates_by_length=barcode_candidates_by_length,
        mismatch_by_barcode=mismatch_by_barcode,
    )


def _best_candidate(candidates: List[BarcodeBoundaryCandidate]) -> BarcodeBoundaryCandidate:
    """Choose one deterministic candidate from one sample-equivalent candidate list."""
    return min(
        candidates,
        key=lambda item: (item.mismatch_distance, item.slack, -len(item.barcode), item.barcode),
    )


def _match_window_length(
    barcode_candidates_by_length: Dict[int, frozenset[bytes]],
    fallback_lengths: Tuple[int, ...],
    cutters: Sequence[bytes],
) -> int:
    """Return the maximum read prefix length needed for barcode+cut matching."""
    barcode_lengths = tuple(barcode_candidates_by_length) or fallback_lengths
    return max(barcode_lengths) + 1 + max(len(cut) for cut in cutters)


def _trim_fastq_record(
    read: Tuple[bytes, bytes, bytes, bytes] | tuple[()],
    trim_start: int,
) -> List[bytes] | tuple[()]:
    """Trim one FASTQ record at the matched barcode boundary."""
    if not read:
        return read
    return [read[0], read[1][trim_start:], read[2], read[3][trim_start:]]


def _format_boundary_ambiguity_label(sample_candidates: Dict[str, List[BarcodeBoundaryCandidate]]) -> bytes:
    """Encode one deterministic ambiguity label for stats/reporting."""
    parts = []
    for sample_name, candidates in sorted(sample_candidates.items()):
        barcodes = sorted({candidate.barcode.decode() for candidate in candidates})
        parts.append(f"{sample_name}:{'|'.join(barcodes)}")
    return ("boundary_ambiguous:" + ";".join(parts)).encode()


def _record_boundary_ambiguity(
    ambiguities: Dict[bytes, int],
    sample_candidates: Dict[str, List[BarcodeBoundaryCandidate]],
) -> None:
    """Increment one deterministic boundary-ambiguity observation."""
    label = _format_boundary_ambiguity_label(sample_candidates)
    ambiguities[label] = ambiguities.get(label, 0) + 1


def _record_barcode_miss(
    barcode_misses: Dict[bytes, int],
    candidates: Sequence[BarcodeBoundaryCandidate],
) -> None:
    """Increment one single-end barcode miss observation."""
    barcode = candidates[0].barcode if candidates else b"XXX"
    barcode_misses[barcode] = barcode_misses.get(barcode, 0) + 1


def _record_combined_barcode_miss(
    barcode_misses: Dict[bytes, int],
    candidates_r1: Sequence[BarcodeBoundaryCandidate],
    candidates_r2: Sequence[BarcodeBoundaryCandidate],
) -> None:
    """Increment one paired-inline barcode miss observation."""
    match_r1 = candidates_r1[0].barcode if candidates_r1 else b"XXX"
    match_r2 = candidates_r2[0].barcode if candidates_r2 else b"XXX"
    barcode = match_r1 + b"_" + match_r2
    barcode_misses[barcode] = barcode_misses.get(barcode, 0) + 1


def _merge_technical_replicate_chunks(
    read1s: Dict[str, List[bytes]],
    read2s: Dict[str, List[bytes]],
) -> None:
    """Merge buffered technical-replicate chunks into their base sample names."""
    replicate_groups: Dict[str, List[str]] = {}
    for sample_name in sorted(read1s):
        if not is_technical_replicate_name(sample_name):
            continue
        base_name = technical_replicate_base_name(sample_name)
        replicate_groups.setdefault(base_name, []).append(sample_name)

    for base_name, replicate_names in replicate_groups.items():
        merged_r1 = read1s.setdefault(base_name, [])
        for replicate_name in replicate_names:
            merged_r1.extend(read1s.pop(replicate_name))
        if not read2s:
            continue
        merged_r2 = read2s.setdefault(base_name, [])
        for replicate_name in replicate_names:
            replicate_reads = read2s.pop(replicate_name, None)
            if replicate_reads:
                merged_r2.extend(replicate_reads)


def write(path: Path, data: List[bytes]) -> None:
    with gzip.open(path, "ab") as out:
        out.write(b"".join(data))


def _has_combinatorial_barcodes(barcodes_to_names: Dict[bytes, str]) -> bool:
    """Return True when barcode keys encode paired inline barcodes."""
    try:
        first = next(iter(barcodes_to_names))
    except StopIteration:
        return False
    return b"_" in first


def get_demux_mode_label(config: DemuxRunConfig) -> str:
    """Return a human-readable demux mode description."""
    if config.i7:
        return "i7 index"
    if _has_combinatorial_barcodes(config.barcodes_to_names):
        return "R1+R2 inline barcodes"
    return "R1 inline barcodes"


def build_matcher(
    fastq_tuple: Tuple[Path, Path | None],
    config: DemuxRunConfig,
    workers: int = 1,
) -> BarMatching:
    """Construct the correct barcode matcher for a raw FASTQ tuple."""
    kwargs = dict(
        fastqs=fastq_tuple,
        barcodes_to_names=config.barcodes_to_names,
        barcodes_to_samples=config.barcodes_to_samples,
        barcode1_to_samples=config.barcode1_to_samples,
        barcode2_to_samples=config.barcode2_to_samples,
        barcode1_candidates_by_length=config.barcode1_candidates_by_length,
        barcode2_candidates_by_length=config.barcode2_candidates_by_length,
        barcode1_mismatch_by_barcode=config.barcode1_mismatch_by_barcode,
        barcode2_mismatch_by_barcode=config.barcode2_mismatch_by_barcode,
        barcode_lengths1=config.barcode_lengths1,
        barcode_lengths2=config.barcode_lengths2,
        cuts1=config.cuts1,
        cuts2=config.cuts2,
        merge_technical_replicates=config.merge_technical_replicates,
        outdir=config.outdir,
        log_level=config.log_level,
        chunksize=config.chunksize,
        max_reads=config.max_reads,
        workers=workers,
    )
    if config.i7:
        return BarMatchingI7(**kwargs)
    if _has_combinatorial_barcodes(config.barcodes_to_names):
        return BarMatchingCombinatorialInline(**kwargs)
    return BarMatchingSingleInline(**kwargs)


def run_serial_demux(
    fastq_tuple: Tuple[Path, Path | None],
    config: DemuxRunConfig,
    workers: int,
) -> Tuple[Dict[bytes, int], Dict[bytes, int], Dict[str, int], Dict[bytes, int]]:
    """Run one serial demux job through the shared matcher builder."""
    barmatcher = build_matcher(fastq_tuple, config, workers=workers)
    try:
        barmatcher.run()
    except KeyboardInterrupt:
        logger.warning("interrupted by user. Shutting down.")
        raise
    return (
        barmatcher.barcode_misses,
        barmatcher.barcode_hits,
        barmatcher.sample_hits,
        barmatcher.barcode_boundary_ambiguities,
    )
