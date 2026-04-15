#!/usr/bin/env python

"""Write the final denovo pseudoreference from split locus tables."""

from __future__ import annotations

import csv
from collections import OrderedDict
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
import time
from typing import Iterator, Sequence, TextIO

from loguru import logger

from ipyrad2.utils.progress import ProgressBar
from ipyrad2.utils.parallel import PipelineTimeoutError, run_pipeline
from ipyrad2.utils.parallel import pipeline as pipeline_module

from .common import (
    DENOVO_MAPPING_FILENAME,
    OUTPUT_JOINED_SPACER_LEN,
    split_cluster_sequence_at_boundary,
)


DEFAULT_MAFFT_TIMEOUT_SECONDS = 900.0
STALL_HEARTBEAT_SECONDS = 60.0


@dataclass(frozen=True)
class AlignmentRunSummary:
    """Summary of how denovo alignment work was scheduled and completed."""

    total_loci: int
    single_sequence_loci: int
    identical_sequence_loci: int
    mafft_required_loci: int
    mafft_threads_per_job: int
    mafft_worker_processes: int
    alignment_mode: str
    mafft_timeout_seconds: int
    joined_spacer_loci: int
    mixed_reconciled_spacer_loci: int
    stripped_output_loci: int
    output_spacer_length: int


@dataclass(frozen=True)
class SummaryRecord:
    """One raw consensus record loaded from concat.summary.tsv."""

    seed: str
    sample: str
    record_type: str
    cluster_sequence: str
    left_arm: str
    right_arm: str


@dataclass(frozen=True)
class LocusMember:
    """One mapping row enriched with summary metadata."""

    locus_id: int
    locus_name: str
    contract_group: str
    core: str
    sample: str
    record_type: str
    cluster_sequence: str
    left_arm: str
    right_arm: str
    reconcile_mode: str
    output_form: str


@dataclass(frozen=True)
class CollapsedGroup:
    """One contract-group reduced to a single sample-local spaced/stripped record."""

    contract_group: str
    sample: str
    joined_only: bool
    stripped: str
    left_arm: str
    right_arm: str


@dataclass(frozen=True)
class _AlignmentPlanItem:
    """Planned handling mode for one locus in submit order."""

    submit_idx: int
    locus_id: int
    locus_name: str
    requires_mafft: bool
    consensus: str
    uses_output_spacer: bool


@dataclass(frozen=True)
class _ThreadJobInfo:
    """Metadata tracked for one in-flight thread worker."""

    submit_idx: int
    locus_id: int
    started_at: float


@dataclass(frozen=True)
class _PlannedAlignment:
    """Alignment plan plus the summary counters derived from it."""

    plans: tuple[_AlignmentPlanItem, ...]
    single_sequence_loci: int
    identical_sequence_loci: int
    joined_spacer_loci: int
    mixed_reconciled_spacer_loci: int
    stripped_output_loci: int


LegacyRecord = Sequence[tuple[str, str]]
ConsensusRecord = Sequence[LocusMember] | LegacyRecord
ThreadedConsensusResult = tuple[int, str, str, bool]


def fasta_text(record: LegacyRecord) -> str:
    """Render a FASTA string from `(name, seq)` pairs with no line wrapping."""
    return "".join(f">{header}\n{seq}\n" for header, seq in record)


def choose_longest_sequence(record: LegacyRecord) -> str:
    """Return the longest sequence in input order, using first-seen tie breaks."""
    if not record:
        raise RuntimeError("record is empty.")
    best_seq = record[0][1]
    best_len = len(best_seq)
    for _name, seq in record[1:]:
        if len(seq) > best_len:
            best_seq = seq
            best_len = len(seq)
    return best_seq


def consensus_from_aligned(aligned: LegacyRecord, min_prop: float = 0.5) -> str:
    """Build a gap-aware majority-base consensus from an alignment."""
    if not aligned:
        return ""
    length = len(aligned[0][1])
    if any(len(seq) != length for _name, seq in aligned):
        raise ValueError("Aligned sequences differ in length.")
    out: list[str] = []
    for idx in range(length):
        bases = [
            seq[idx].upper()
            for _name, seq in aligned
            if seq[idx] != "-" and seq[idx].upper() in "ACGT"
        ]
        if not bases:
            continue
        counts = {base: bases.count(base) for base in "ACGT"}
        base, _nbase = max(counts.items(), key=lambda item: item[1])
        out.append(base)
    return "".join(out)


def _record_length_summary(record: LegacyRecord) -> tuple[int, int, int]:
    """Return `(nseq, min_len, max_len)` for one sequence record."""
    lengths = [len(seq) for _name, seq in record]
    return len(lengths), min(lengths), max(lengths)


def mafft_align_one(
    record: LegacyRecord,
    mafft_binary: str,
    threads: int = 1,
    *,
    locus_id: int | None = None,
    timeout_s: float = DEFAULT_MAFFT_TIMEOUT_SECONDS,
) -> list[tuple[str, str]]:
    """Run MAFFT on one record via stdin and return aligned sequences."""
    if not record:
        raise RuntimeError("record is empty.")
    argv = [
        mafft_binary,
        "--quiet",
        "--auto",
        "--thread",
        str(threads),
        "--adjustdirection",
        "-",
    ]
    try:
        rc, out, err = run_pipeline(
            [argv],
            outfile=None,
            stdin_text=fasta_text(record),
            timeout_s=timeout_s,
        )
    except PipelineTimeoutError as exc:
        nseq, min_len, max_len = _record_length_summary(record)
        locus_label = f"locus_{locus_id}" if locus_id is not None else "unknown_locus"
        raise RuntimeError(
            f"mafft timed out for {locus_label} after {timeout_s:.0f}s "
            f"(nseq={nseq}, min_len={min_len}, max_len={max_len})"
        ) from exc
    if rc != 0 or not out:
        raise RuntimeError(
            f"mafft failed rc={rc} stderr={err.decode('utf-8', 'replace')}"
        )

    aligned: list[tuple[str, str]] = []
    name, buf = None, []
    for line in out.decode("utf-8", "replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                aligned.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
            continue
        buf.append(line.strip())
    if name is not None:
        aligned.append((name, "".join(buf)))
    return aligned


def _validate_alignment_mode(alignment_mode: str) -> str:
    """Validate the final-reference alignment mode."""
    if alignment_mode not in {"mafft", "none"}:
        raise ValueError(f"Unsupported alignment_mode: {alignment_mode}")
    return alignment_mode


def _all_sequences_identical(record: LegacyRecord) -> bool:
    """Return True when every sequence in one record is identical."""
    if not record:
        return False
    first = record[0][1]
    return all(seq == first for _name, seq in record[1:])


def _collapse_sequence_record_without_mafft(
    record: LegacyRecord,
) -> str | None:
    """Return a sequence only when MAFFT is provably unnecessary."""
    if not record:
        raise RuntimeError("record is empty.")
    if len(record) == 1:
        return record[0][1]
    if _all_sequences_identical(record):
        return record[0][1]
    return None


def _collapse_sequence_record(
    record: LegacyRecord,
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
) -> str:
    """Collapse one sequence record using the current alignment policy."""
    alignment_mode = _validate_alignment_mode(alignment_mode)
    if alignment_mode == "none":
        return choose_longest_sequence(record)
    consensus = _collapse_sequence_record_without_mafft(record)
    if consensus is not None:
        return consensus
    aligned = mafft_align_one(
        record,
        mafft_binary=mafft_binary,
        threads=threads,
        locus_id=locus_id,
        timeout_s=timeout_s,
    )
    return consensus_from_aligned(aligned, min_prop=min_prop)


def _load_summary_records(summary_tsv: Path) -> dict[str, SummaryRecord]:
    """Return summary records keyed by seed/core name."""
    out: dict[str, SummaryRecord] = {}
    with open(summary_tsv, "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        required = {"seed", "sample", "record_type", "cluster_sequence", "arm_boundary"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            joined = ", ".join(sorted(required))
            raise RuntimeError(
                f"concat.summary.tsv is missing required columns: {joined}"
            )
        for row in reader:
            seed = str(row["seed"])
            sample = str(row["sample"])
            record_type = str(row["record_type"])
            cluster_sequence = str(row["cluster_sequence"]).upper()
            left_arm, right_arm = split_cluster_sequence_at_boundary(
                cluster_sequence,
                int(row["arm_boundary"]),
            )
            out[seed] = SummaryRecord(
                seed=seed,
                sample=sample,
                record_type=record_type,
                cluster_sequence=cluster_sequence,
                left_arm=str(left_arm).upper(),
                right_arm=str(right_arm).upper(),
            )
    return out


def _iter_locus_members(
    mapping_tsv: Path,
    summary_records: dict[str, SummaryRecord],
) -> Iterator[tuple[int, str, list[LocusMember]]]:
    """Yield one enriched locus at a time in mapping order."""
    current_locus: int | None = None
    current_name: str | None = None
    current_members: list[LocusMember] = []

    with open(mapping_tsv, "rt", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if (
            reader.fieldnames is None
            or "locus" not in reader.fieldnames
            or "core" not in reader.fieldnames
        ):
            raise RuntimeError(
                f"{DENOVO_MAPPING_FILENAME} is missing required columns: locus, core"
            )
        for row in reader:
            locus_id = int(row["locus"])
            core = str(row["core"])
            record = summary_records[core]
            locus_name = str(row.get("locus_name") or f"locus_{locus_id}")
            contract_group = str(row.get("contract_group") or core)
            if current_locus is None:
                current_locus = locus_id
                current_name = locus_name
            if locus_id != current_locus:
                assert current_name is not None
                yield current_locus, current_name, current_members
                current_locus = locus_id
                current_name = locus_name
                current_members = []
            current_members.append(
                LocusMember(
                    locus_id=locus_id,
                    locus_name=locus_name,
                    contract_group=contract_group,
                    core=core,
                    sample=record.sample,
                    record_type=record.record_type,
                    cluster_sequence=record.cluster_sequence,
                    left_arm=record.left_arm,
                    right_arm=record.right_arm,
                    reconcile_mode=str(row.get("reconcile_mode") or "none"),
                    output_form=str(row.get("output_form") or ""),
                )
            )

    if current_locus is not None:
        assert current_name is not None
        yield current_locus, current_name, current_members


def _group_members_by_contract_group(
    members: Sequence[LocusMember],
) -> list[list[LocusMember]]:
    """Return members grouped by contract-group id in first-seen order."""
    groups: OrderedDict[str, list[LocusMember]] = OrderedDict()
    for member in members:
        groups.setdefault(member.contract_group, []).append(member)
    return list(groups.values())


def _group_is_joined_only(group: Sequence[LocusMember]) -> bool:
    """Return True when every member is a joined record with recoverable arms."""
    return bool(group) and all(
        member.record_type == "joined" and bool(member.right_arm) for member in group
    )


def _locus_output_form(members: Sequence[LocusMember]) -> str:
    """Return the configured output form for one locus."""
    if not members:
        return "stripped"
    forms = {member.output_form for member in members if member.output_form}
    if not forms:
        return (
            "spaced"
            if all(
                member.record_type == "joined" and bool(member.right_arm)
                for member in members
            )
            else "stripped"
        )
    if len(forms) != 1:
        raise RuntimeError("locus rows disagree on output_form")
    return next(iter(forms))


def _locus_is_mixed_spaced(members: Sequence[LocusMember]) -> bool:
    """Return True when one locus is a mixed joined/merged spaced output."""
    return _locus_output_form(members) == "spaced" and not all(
        member.record_type == "joined" and bool(member.right_arm) for member in members
    )


def _collapse_member_group_without_mafft(
    group: Sequence[LocusMember],
) -> CollapsedGroup | None:
    """Collapse one contract group only when MAFFT is unnecessary."""
    if not group:
        raise RuntimeError("contract group is empty")
    joined_only = _group_is_joined_only(group)
    if joined_only:
        left = _collapse_sequence_record_without_mafft(
            [(m.core, m.left_arm) for m in group]
        )
        right = _collapse_sequence_record_without_mafft(
            [(m.core, m.right_arm) for m in group]
        )
        if left is None or right is None:
            return None
        return CollapsedGroup(
            contract_group=group[0].contract_group,
            sample=group[0].sample,
            joined_only=True,
            stripped=left + right,
            left_arm=left,
            right_arm=right,
        )
    stripped = _collapse_sequence_record_without_mafft(
        [(m.core, m.cluster_sequence) for m in group]
    )
    if stripped is None:
        return None
    return CollapsedGroup(
        contract_group=group[0].contract_group,
        sample=group[0].sample,
        joined_only=False,
        stripped=stripped,
        left_arm=stripped,
        right_arm="",
    )


def _collapse_member_group(
    group: Sequence[LocusMember],
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
) -> CollapsedGroup:
    """Collapse one contract group under the current alignment policy."""
    if not group:
        raise RuntimeError("contract group is empty")
    joined_only = _group_is_joined_only(group)
    if joined_only:
        left = _collapse_sequence_record(
            [(m.core, m.left_arm) for m in group],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        right = _collapse_sequence_record(
            [(m.core, m.right_arm) for m in group],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        return CollapsedGroup(
            contract_group=group[0].contract_group,
            sample=group[0].sample,
            joined_only=True,
            stripped=left + right,
            left_arm=left,
            right_arm=right,
        )
    stripped = _collapse_sequence_record(
        [(m.core, m.cluster_sequence) for m in group],
        locus_id=locus_id,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        timeout_s=timeout_s,
    )
    return CollapsedGroup(
        contract_group=group[0].contract_group,
        sample=group[0].sample,
        joined_only=False,
        stripped=stripped,
        left_arm=stripped,
        right_arm="",
    )


def _locus_uses_output_spacer(members: Sequence[LocusMember]) -> bool:
    """Return True when the final locus should be written with the output spacer."""
    return _locus_output_form(members) == "spaced"


def _aligned_breakpoint_column(aligned_seq: str, left_len: int) -> int:
    """Return the aligned column index immediately after one left arm."""
    if left_len <= 0:
        return 0
    seen = 0
    for idx, char in enumerate(aligned_seq):
        if char != "-":
            seen += 1
            if seen >= left_len:
                return idx + 1
    return len(aligned_seq)


def _infer_mixed_locus_boundary(
    members: Sequence[LocusMember],
    aligned_map: dict[str, str],
) -> int | None:
    """Infer one coherent arm boundary for a mixed spaced locus."""
    boundaries = [
        _aligned_breakpoint_column(aligned_map[member.core], len(member.left_arm))
        for member in members
        if member.record_type == "joined" and bool(member.right_arm)
    ]
    if not boundaries:
        return None
    boundaries = sorted(boundaries)
    return boundaries[len(boundaries) // 2]


def _strip_alignment_segment(segment: str) -> str:
    """Remove gaps from one aligned segment while preserving base order."""
    return "".join(char.upper() for char in segment if char != "-")


def _collapse_optional_sequence_record(
    record: LegacyRecord,
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
) -> str:
    """Collapse one sequence record, returning an empty string when no bases exist."""
    nonempty = [(name, seq) for name, seq in record if seq]
    if not nonempty:
        return ""
    return _collapse_sequence_record(
        nonempty,
        locus_id=locus_id,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        timeout_s=timeout_s,
    )


def _consensus_for_mixed_spaced_locus(
    members: Sequence[LocusMember],
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
    output_spacer_len: int,
) -> tuple[str, bool]:
    """Return a spaced consensus for one mixed joined/merged reconciled locus."""
    joined_members = [
        member
        for member in members
        if member.record_type == "joined" and bool(member.right_arm)
    ]
    if not joined_members:
        stripped = _collapse_sequence_record(
            [(member.sample, member.cluster_sequence) for member in members],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        return stripped, False

    if alignment_mode == "none":
        left = choose_longest_sequence(
            [(member.core, member.left_arm) for member in joined_members]
        )
        right = choose_longest_sequence(
            [(member.core, member.right_arm) for member in joined_members]
        )
        return left + ("N" * output_spacer_len) + right, True

    aligned_map = dict(
        mafft_align_one(
            [(member.core, member.cluster_sequence) for member in members],
            mafft_binary=mafft_binary,
            threads=threads,
            locus_id=locus_id,
            timeout_s=timeout_s,
        )
    )
    boundary = _infer_mixed_locus_boundary(members, aligned_map)
    if boundary is None:
        stripped = _collapse_sequence_record(
            [(member.sample, member.cluster_sequence) for member in members],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        return stripped, False

    grouped_members = _group_members_by_contract_group(members)
    groups: list[CollapsedGroup] = []
    for group in grouped_members:
        left_record = []
        right_record = []
        for member in group:
            aligned_seq = aligned_map[member.core]
            left_record.append(
                (member.core, _strip_alignment_segment(aligned_seq[:boundary]))
            )
            right_record.append(
                (member.core, _strip_alignment_segment(aligned_seq[boundary:]))
            )
        left = _collapse_optional_sequence_record(
            left_record,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        right = _collapse_optional_sequence_record(
            right_record,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        groups.append(
            CollapsedGroup(
                contract_group=group[0].contract_group,
                sample=group[0].sample,
                joined_only=True,
                stripped=left + right,
                left_arm=left,
                right_arm=right,
            )
        )

    left = _collapse_optional_sequence_record(
        [(group.sample, group.left_arm) for group in groups],
        locus_id=locus_id,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        timeout_s=timeout_s,
    )
    right = _collapse_optional_sequence_record(
        [(group.sample, group.right_arm) for group in groups],
        locus_id=locus_id,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        timeout_s=timeout_s,
    )
    return left + ("N" * output_spacer_len) + right, True


def _plan_locus_without_mafft(
    members: Sequence[LocusMember],
    *,
    output_spacer_len: int,
) -> tuple[str, str, bool] | None:
    """Return `(reason, consensus, uses_output_spacer)` when MAFFT is unnecessary."""
    if _locus_is_mixed_spaced(members):
        return None
    groups = _group_members_by_contract_group(members)
    collapsed = []
    for group in groups:
        reduced = _collapse_member_group_without_mafft(group)
        if reduced is None:
            return None
        collapsed.append(reduced)

    uses_output_spacer = bool(collapsed) and all(
        group.joined_only for group in collapsed
    )
    if uses_output_spacer:
        left = _collapse_sequence_record_without_mafft(
            [(group.sample, group.left_arm) for group in collapsed]
        )
        right = _collapse_sequence_record_without_mafft(
            [(group.sample, group.right_arm) for group in collapsed]
        )
        if left is None or right is None:
            return None
        reason = "single" if len(collapsed) == 1 else "identical"
        return reason, left + ("N" * output_spacer_len) + right, True

    stripped = _collapse_sequence_record_without_mafft(
        [(group.sample, group.stripped) for group in collapsed]
    )
    if stripped is None:
        return None
    reason = "single" if len(collapsed) == 1 else "identical"
    return reason, stripped, False


def _consensus_for_locus_members(
    members: Sequence[LocusMember],
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
    output_spacer_len: int,
) -> tuple[str, bool]:
    """Return one final locus sequence plus whether the output spacer was used."""
    if _locus_is_mixed_spaced(members):
        return _consensus_for_mixed_spaced_locus(
            members,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
            output_spacer_len=output_spacer_len,
        )
    groups = [
        _collapse_member_group(
            group,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        for group in _group_members_by_contract_group(members)
    ]
    uses_output_spacer = bool(groups) and all(group.joined_only for group in groups)
    if uses_output_spacer:
        left = _collapse_sequence_record(
            [(group.sample, group.left_arm) for group in groups],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        right = _collapse_sequence_record(
            [(group.sample, group.right_arm) for group in groups],
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
        return left + ("N" * output_spacer_len) + right, True
    stripped = _collapse_sequence_record(
        [(group.sample, group.stripped) for group in groups],
        locus_id=locus_id,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        timeout_s=timeout_s,
    )
    return stripped, False


def _consensus_for_legacy_record(
    record: LegacyRecord,
    *,
    locus_id: int,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    timeout_s: float,
) -> tuple[str, bool]:
    """Return the old stripped-sequence consensus path for direct tests."""
    return (
        _collapse_sequence_record(
            record,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        ),
        False,
    )


def worker_build_consensus(
    locus_id: int,
    record: ConsensusRecord,
    mafft_binary: str,
    min_prop: float = 0.5,
    threads: int = 1,
    alignment_mode: str = "mafft",
    timeout_s: float = DEFAULT_MAFFT_TIMEOUT_SECONDS,
    locus_name: str | None = None,
    output_spacer_len: int = OUTPUT_JOINED_SPACER_LEN,
) -> tuple[int, str, str, bool]:
    """Build one final locus sequence for the denovo reference."""
    if record and isinstance(record[0], LocusMember):
        consensus, uses_output_spacer = _consensus_for_locus_members(
            record,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
            output_spacer_len=output_spacer_len,
        )
    else:
        consensus, uses_output_spacer = _consensus_for_legacy_record(
            record,
            locus_id=locus_id,
            mafft_binary=mafft_binary,
            min_prop=min_prop,
            threads=threads,
            alignment_mode=alignment_mode,
            timeout_s=timeout_s,
        )
    return locus_id, (locus_name or f"locus_{locus_id}"), consensus, uses_output_spacer


def _plan_alignment(
    mapping_tsv: Path,
    summary_records: dict[str, SummaryRecord],
    *,
    alignment_mode: str,
    output_spacer_len: int,
) -> _PlannedAlignment:
    """Plan how each locus will be handled and count trivial cases."""
    plans: list[_AlignmentPlanItem] = []
    single_sequence_loci = 0
    identical_sequence_loci = 0
    mixed_reconciled_spacer_loci = 0
    joined_spacer_loci = 0
    stripped_output_loci = 0

    for submit_idx, (locus_id, locus_name, members) in enumerate(
        _iter_locus_members(mapping_tsv, summary_records)
    ):
        uses_output_spacer = _locus_uses_output_spacer(members)
        if uses_output_spacer:
            if _locus_is_mixed_spaced(members):
                mixed_reconciled_spacer_loci += 1
            else:
                joined_spacer_loci += 1
        else:
            stripped_output_loci += 1

        if alignment_mode == "none":
            consensus, _uses_output_spacer = _consensus_for_locus_members(
                members,
                locus_id=locus_id,
                mafft_binary="",
                min_prop=0.5,
                threads=1,
                alignment_mode="none",
                timeout_s=DEFAULT_MAFFT_TIMEOUT_SECONDS,
                output_spacer_len=output_spacer_len,
            )
            plans.append(
                _AlignmentPlanItem(
                    submit_idx=submit_idx,
                    locus_id=locus_id,
                    locus_name=locus_name,
                    requires_mafft=False,
                    consensus=consensus,
                    uses_output_spacer=_uses_output_spacer,
                )
            )
            continue

        planned = _plan_locus_without_mafft(
            members, output_spacer_len=output_spacer_len
        )
        if planned is None:
            plans.append(
                _AlignmentPlanItem(
                    submit_idx=submit_idx,
                    locus_id=locus_id,
                    locus_name=locus_name,
                    requires_mafft=True,
                    consensus="",
                    uses_output_spacer=uses_output_spacer,
                )
            )
            continue

        reason, consensus, planned_output_spacer = planned
        if reason == "single":
            single_sequence_loci += 1
        else:
            identical_sequence_loci += 1
        plans.append(
            _AlignmentPlanItem(
                submit_idx=submit_idx,
                locus_id=locus_id,
                locus_name=locus_name,
                requires_mafft=False,
                consensus=consensus,
                uses_output_spacer=planned_output_spacer,
            )
        )
    return _PlannedAlignment(
        plans=tuple(plans),
        single_sequence_loci=single_sequence_loci,
        identical_sequence_loci=identical_sequence_loci,
        joined_spacer_loci=joined_spacer_loci,
        mixed_reconciled_spacer_loci=mixed_reconciled_spacer_loci,
        stripped_output_loci=stripped_output_loci,
    )


def _iter_alignment_jobs_from_plan(
    mapping_tsv: Path,
    summary_records: dict[str, SummaryRecord],
    plans: Sequence[_AlignmentPlanItem],
    *,
    mafft_binary: str,
    min_prop: float,
    threads: int,
    alignment_mode: str,
    mafft_timeout_s: float,
    output_spacer_len: int,
) -> Iterator[tuple[int, dict[str, object]]]:
    """Yield per-locus MAFFT jobs reusing a previously computed plan."""
    if _validate_alignment_mode(alignment_mode) == "none":
        return

    plan_iter = iter(plans)
    for locus_id, locus_name, members in _iter_locus_members(
        mapping_tsv, summary_records
    ):
        plan = next(plan_iter)
        if plan.locus_id != locus_id:
            raise RuntimeError(
                f"alignment plan drifted from {DENOVO_MAPPING_FILENAME} order"
            )
        if not plan.requires_mafft:
            continue
        yield (
            plan.submit_idx,
            dict(
                locus_id=locus_id,
                locus_name=locus_name,
                record=members,
                mafft_binary=mafft_binary,
                min_prop=min_prop,
                threads=threads,
                alignment_mode=alignment_mode,
                timeout_s=mafft_timeout_s,
                output_spacer_len=output_spacer_len,
            ),
        )


def iter_alignment_jobs(
    mapping_tsv: Path,
    summary_tsv: Path,
    mafft_binary: str,
    min_prop: float = 0.5,
    threads: int = 1,
    alignment_mode: str = "mafft",
    mafft_timeout_s: float = DEFAULT_MAFFT_TIMEOUT_SECONDS,
    output_spacer_len: int = OUTPUT_JOINED_SPACER_LEN,
) -> Iterator[tuple[int, dict[str, object]]]:
    """Yield ordered per-locus MAFFT jobs for loci that require alignment."""
    summary_records = _load_summary_records(summary_tsv)
    planned = _plan_alignment(
        mapping_tsv,
        summary_records,
        alignment_mode=alignment_mode,
        output_spacer_len=output_spacer_len,
    )
    yield from _iter_alignment_jobs_from_plan(
        mapping_tsv,
        summary_records,
        planned.plans,
        mafft_binary=mafft_binary,
        min_prop=min_prop,
        threads=threads,
        alignment_mode=alignment_mode,
        mafft_timeout_s=mafft_timeout_s,
        output_spacer_len=output_spacer_len,
    )


def _iter_threaded_alignment_results(
    jobs_iter: Iterator[tuple[int, dict[str, object]]],
    *,
    max_workers: int,
    heartbeat_s: float,
) -> Iterator[tuple[int, ThreadedConsensusResult]]:
    """Yield completed threaded locus results."""
    inflight: dict[Future[ThreadedConsensusResult], _ThreadJobInfo] = {}
    job_it = iter(jobs_iter)
    executor = ThreadPoolExecutor(
        max_workers=max_workers, thread_name_prefix="denovo-align"
    )
    close_wait = True

    def _submit_next() -> bool:
        try:
            submit_idx, kwargs = next(job_it)
        except StopIteration:
            return False
        future = executor.submit(worker_build_consensus, **kwargs)
        inflight[future] = _ThreadJobInfo(
            submit_idx=submit_idx,
            locus_id=int(kwargs["locus_id"]),
            started_at=time.monotonic(),
        )
        return True

    while len(inflight) < max_workers and _submit_next():
        pass

    next_heartbeat = time.monotonic() + heartbeat_s
    try:
        while inflight:
            timeout = max(0.0, next_heartbeat - time.monotonic())
            done, _pending = wait(
                tuple(inflight), timeout=timeout, return_when=FIRST_COMPLETED
            )
            if not done:
                oldest = min(inflight.values(), key=lambda info: info.started_at)
                logger.debug(
                    f"alignment heartbeat: completed=0 inflight={len(inflight)} "
                    f"oldest_locus=locus_{oldest.locus_id}"
                )
                next_heartbeat = time.monotonic() + heartbeat_s
                continue

            for future in done:
                info = inflight.pop(future)
                try:
                    result = future.result()
                except BaseException:
                    pipeline_module._kill_all_children()
                    executor.shutdown(wait=False, cancel_futures=True)
                    close_wait = False
                    raise
                yield info.submit_idx, result
                while len(inflight) < max_workers and _submit_next():
                    pass
            next_heartbeat = time.monotonic() + heartbeat_s
    except KeyboardInterrupt:
        pipeline_module._kill_all_children()
        executor.shutdown(wait=False, cancel_futures=True)
        close_wait = False
        raise
    finally:
        if close_wait:
            executor.shutdown(wait=True, cancel_futures=False)


def _choose_mafft_runtime(
    cores: int,
    mafft_required_loci: int,
) -> tuple[int, int]:
    """Choose MAFFT worker and thread counts from actual alignable loci."""
    if mafft_required_loci <= 0:
        return 0, 0
    workers = min(cores, mafft_required_loci)
    threads = max(1, cores // workers)
    return workers, threads


def _build_alignment_run_summary(
    planned: _PlannedAlignment,
    *,
    alignment_mode: str,
    mafft_required_loci: int,
    mafft_workers: int,
    mafft_threads: int,
    mafft_timeout_s: float,
    output_spacer_len: int,
) -> AlignmentRunSummary:
    """Return the final alignment scheduling summary."""
    return AlignmentRunSummary(
        total_loci=len(planned.plans),
        single_sequence_loci=planned.single_sequence_loci,
        identical_sequence_loci=planned.identical_sequence_loci,
        mafft_required_loci=mafft_required_loci,
        mafft_threads_per_job=mafft_threads,
        mafft_worker_processes=mafft_workers,
        alignment_mode=alignment_mode,
        mafft_timeout_seconds=int(mafft_timeout_s) if mafft_required_loci else 0,
        joined_spacer_loci=planned.joined_spacer_loci,
        mixed_reconciled_spacer_loci=planned.mixed_reconciled_spacer_loci,
        stripped_output_loci=planned.stripped_output_loci,
        output_spacer_length=(
            output_spacer_len
            if (planned.joined_spacer_loci or planned.mixed_reconciled_spacer_loci)
            else 0
        ),
    )


def _log_alignment_summary(summary: AlignmentRunSummary) -> None:
    """Log the alignment scheduling summary in one compact line."""
    logger.debug(
        "alignment scheduling: "
        f"total_loci={summary.total_loci} "
        f"single_sequence_loci={summary.single_sequence_loci} "
        f"identical_sequence_loci={summary.identical_sequence_loci} "
        f"mafft_required_loci={summary.mafft_required_loci} "
        f"mafft_worker_processes={summary.mafft_worker_processes} "
        f"mafft_threads_per_job={summary.mafft_threads_per_job} "
        f"mafft_timeout_seconds={summary.mafft_timeout_seconds} "
        f"joined_spacer_loci={summary.joined_spacer_loci} "
        f"mixed_reconciled_spacer_loci={summary.mixed_reconciled_spacer_loci} "
        f"stripped_output_loci={summary.stripped_output_loci}"
    )


def _progress_message(total_loci: int, alignment_mode: str) -> str:
    """Return the progress-bar label for final reference writing."""
    action = "Writing loci" if alignment_mode == "none" else "Aligning loci"
    return f"{action} - total jobs: {total_loci}"


def _write_planned_locus(
    fh: TextIO,
    *,
    locus_name: str,
    consensus: str,
) -> None:
    """Write one FASTA record for the final denovo reference."""
    fh.write(f">{locus_name}\n{consensus}\n")


def write_ordered_consensus_stream_to_file(
    mapping_tsv: Path,
    summary_tsv: Path,
    out_fa: Path,
    mafft_binary: str,
    cores: int = 1,
    min_prop: float = 0.5,
    alignment_mode: str = "mafft",
    mafft_timeout_s: float = DEFAULT_MAFFT_TIMEOUT_SECONDS,
    output_spacer_len: int = OUTPUT_JOINED_SPACER_LEN,
) -> AlignmentRunSummary:
    """Write the denovo reference FASTA in mapping order and return runtime summary."""
    alignment_mode = _validate_alignment_mode(alignment_mode)
    if cores < 1:
        raise ValueError("cores must be >= 1")
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    summary_records = _load_summary_records(summary_tsv)
    planned = _plan_alignment(
        mapping_tsv,
        summary_records,
        alignment_mode=alignment_mode,
        output_spacer_len=output_spacer_len,
    )
    mafft_required_loci = sum(1 for plan in planned.plans if plan.requires_mafft)
    mafft_workers, mafft_threads = _choose_mafft_runtime(cores, mafft_required_loci)
    summary = _build_alignment_run_summary(
        planned,
        alignment_mode=alignment_mode,
        mafft_required_loci=mafft_required_loci,
        mafft_workers=mafft_workers,
        mafft_threads=mafft_threads,
        mafft_timeout_s=mafft_timeout_s,
        output_spacer_len=output_spacer_len,
    )
    _log_alignment_summary(summary)

    if summary.total_loci == 0:
        out_fa.write_text("", encoding="utf-8")
        logger.info("wrote denovo reference")
        logger.debug("wrote denovo reference to {}", out_fa)
        return summary

    with open(out_fa, "wt", encoding="utf-8") as fh:
        prog = ProgressBar(
            summary.total_loci,
            None,
            _progress_message(summary.total_loci, alignment_mode),
        )
        prog.update()
        try:
            if alignment_mode == "none" or mafft_required_loci == 0:
                for plan in planned.plans:
                    _write_planned_locus(
                        fh, locus_name=plan.locus_name, consensus=plan.consensus
                    )
                    prog.finished += 1
                    prog.update()
                logger.info("wrote denovo reference")
                logger.debug("wrote denovo reference to {}", out_fa)
                return summary

            jobs_it = _iter_alignment_jobs_from_plan(
                mapping_tsv=mapping_tsv,
                summary_records=summary_records,
                plans=planned.plans,
                mafft_binary=mafft_binary,
                min_prop=min_prop,
                threads=mafft_threads,
                alignment_mode=alignment_mode,
                mafft_timeout_s=mafft_timeout_s,
                output_spacer_len=output_spacer_len,
            )
            result_buffer: dict[int, ThreadedConsensusResult] = {}
            next_idx = 0

            def _flush_ready() -> None:
                nonlocal next_idx
                # MAFFT results can finish out of order; buffer them until every
                # earlier locus has been written so the final FASTA stays stable.
                while next_idx < len(planned.plans):
                    plan = planned.plans[next_idx]
                    if not plan.requires_mafft:
                        _write_planned_locus(
                            fh, locus_name=plan.locus_name, consensus=plan.consensus
                        )
                        next_idx += 1
                        prog.finished += 1
                        prog.update()
                        continue
                    if next_idx not in result_buffer:
                        break
                    _locus_id, locus_name, consensus, _uses_output_spacer = (
                        result_buffer.pop(next_idx)
                    )
                    _write_planned_locus(fh, locus_name=locus_name, consensus=consensus)
                    next_idx += 1
                    prog.finished += 1
                    prog.update()

            _flush_ready()
            for submit_idx, result in _iter_threaded_alignment_results(
                jobs_it,
                max_workers=mafft_workers,
                heartbeat_s=STALL_HEARTBEAT_SECONDS,
            ):
                result_buffer[submit_idx] = result
                _flush_ready()
        except KeyboardInterrupt:
            logger.warning("interrupted by user. Cleaning up.")
            raise SystemExit(130)
        finally:
            prog.close()

    logger.info("wrote denovo reference")
    logger.debug("wrote denovo reference to {}", out_fa)
    return summary
