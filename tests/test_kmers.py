import gzip
from pathlib import Path

import pytest

import ipyrad2.utils.kmers as kmers
from ipyrad2.utils.exceptions import IPyradError


@pytest.fixture
def sequential_pool(monkeypatch: pytest.MonkeyPatch) -> None:
    def run_jobs_sequentially(
        jobs,
        log_level,
        max_workers=None,
        max_inflight=None,
        msg="Processing",
    ):
        results = {}
        for key, (func, kwargs) in jobs.items():
            results[key] = func(**kwargs)
        return results

    monkeypatch.setattr(kmers, "run_with_pool", run_jobs_sequentially)


def _write_fastq(path: Path, reads: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8") as out:
        for idx, read in enumerate(reads):
            out.write(f"@r{idx}\n{read}\n+\n{'I' * len(read)}\n")
    return path


def _junction_reads(prefix: str, tail_prefixes: tuple[str, ...]) -> list[str]:
    return [prefix + tail + "ACGTACGT" for tail in tail_prefixes * 5]


def _reads_from_barcodes(
    barcodes_by_length: dict[int, tuple[str, ...]],
    motifs: tuple[str, ...],
    *,
    include_mismatch: bool = False,
) -> list[str]:
    reads: list[str] = []
    for barcodes in barcodes_by_length.values():
        for barcode in barcodes:
            for motif in motifs:
                for tail in ("A", "C", "G", "T") * 4:
                    reads.append(barcode + motif + tail + "ACGTACGT")
    if include_mismatch:
        reads.append("ACGTTCATCGGACGTACGT")
    return reads


def _offset_candidate(
    offset: int,
    motifs: tuple[str, ...],
    counts: tuple[int, ...],
) -> kmers._OffsetInference:
    return kmers._OffsetInference(
        offset=offset,
        motifs=motifs,
        motif_counts=counts,
        total_support=sum(counts),
    )


def test_get_overhang_from_kmers_finds_offset_zero_junction(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        _junction_reads("ATCGG", ("A", "C", "G", "T")),
    )

    inferred = kmers.get_overhang_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0, 1),
    )

    assert inferred.sequence == "ATCGG"
    assert inferred.offset == 0
    assert inferred.trim_length == 5


def test_get_overhang_from_kmers_supports_plain_fastq_input(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq",
        _junction_reads("ATCGG", ("A", "C", "G", "T")),
    )

    inferred = kmers.get_overhang_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0, 1),
    )

    assert inferred.sequence == "ATCGG"
    assert inferred.offset == 0


def test_get_overhang_from_kmers_finds_offset_one_junction(
    tmp_path: Path,
    sequential_pool,
) -> None:
    reads = []
    for lead in ("A", "C", "G", "T") * 5:
        for tail in ("A", "C", "G", "T"):
            reads.append(lead + "ATCGG" + tail + "ACGTACGT")
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    inferred = kmers.get_overhang_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0, 1),
    )

    assert inferred.sequence == "ATCGG"
    assert inferred.offset == 1
    assert inferred.trim_length == 6


def test_get_overhang_from_kmers_discards_monomorphic_candidate_offsets(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        _junction_reads("AAAAAAATCGG", ("A", "C", "G", "T")),
    )

    inferred = kmers.get_overhang_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0, 6),
    )

    assert inferred.sequence == "ATCGG"
    assert inferred.offset == 6


def test_get_overhang_from_kmers_raises_when_only_invalid_motifs_exist(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", ["AAAAAAAAAAAA"] * 20)

    with pytest.raises(IPyradError, match="only invalid or low-information motifs"):
        kmers.get_overhang_from_kmers(
            [fastq],
            20,
            100,
            1,
            "ERROR",
            candidate_offsets=(0, 1),
        )


def test_get_overhang_from_kmers_aggregates_duplicate_basenames(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq1 = _write_fastq(
        tmp_path / "a" / "shared.fastq.gz",
        _junction_reads("ATCGG", ("A", "C", "G", "T")),
    )
    fastq2 = _write_fastq(tmp_path / "b" / "shared.fastq.gz", ["AAAAAAAAAAAA"] * 20)

    inferred = kmers.get_overhang_from_kmers(
        [fastq1, fastq2],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0,),
    )

    assert inferred.sequence == "ATCGG"
    assert inferred.offset == 0


def test_get_overhang_from_kmers_samples_at_least_one_read_per_file(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastqs = [
        _write_fastq(tmp_path / f"sample{i}.fastq.gz", [f"ATCGG{base}ACGTACGT"])
        for i, base in enumerate(("A", "C", "G"))
    ]

    inferred = kmers.get_overhang_from_kmers(
        fastqs,
        20,
        2,
        1,
        "ERROR",
        candidate_offsets=(0,),
    )

    assert inferred.sequence == "ATCGG"


def test_get_overhangs_from_kmers_accepts_multiple_motifs_at_one_offset(
    tmp_path: Path,
    sequential_pool,
) -> None:
    reads = (
        _junction_reads("ATCGG", ("A", "C", "G", "T")) +
        _junction_reads("ATCGAT", ("A", "C", "G", "T"))
    )
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    inferred = kmers.get_overhangs_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0, 1),
    )

    assert inferred.offset == 0
    assert inferred.motifs == ("ATCGAT", "ATCGG")
    assert inferred.trim_length == 6


def test_select_best_offset_candidate_rejects_shifted_near_tie(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    debug_messages: list[str] = []
    trace_messages: list[str] = []
    stub_logger = type(
        "LoggerStub",
        (),
        {
            "debug": staticmethod(lambda *args: debug_messages.append(args[0].format(*args[1:]))),
            "trace": staticmethod(lambda *args: trace_messages.append(args[0].format(*args[1:]))),
        },
    )
    monkeypatch.setattr(kmers, "logger", stub_logger)

    best = kmers._select_best_offset_candidate(
        [
            _offset_candidate(0, ("ATCGG", "ATCGAT"), (233008, 102085)),
            _offset_candidate(1, ("TCGG", "TCGAT"), (233341, 102207)),
        ],
        (0, 1),
    )

    assert best.offset == 0
    assert any("required_margin=10%" in message for message in trace_messages)
    assert any("rejected shifted offset=1" in message for message in debug_messages)


def test_select_best_offset_candidate_accepts_shifted_clear_win(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    debug_messages: list[str] = []
    trace_messages: list[str] = []
    stub_logger = type(
        "LoggerStub",
        (),
        {
            "debug": staticmethod(lambda *args: debug_messages.append(args[0].format(*args[1:]))),
            "trace": staticmethod(lambda *args: trace_messages.append(args[0].format(*args[1:]))),
        },
    )
    monkeypatch.setattr(kmers, "logger", stub_logger)

    best = kmers._select_best_offset_candidate(
        [
            _offset_candidate(0, ("ATCGG", "ATCGAT"), (200, 90)),
            _offset_candidate(1, ("XATCGG", "XATCGAT"), (230, 105)),
        ],
        (0, 1),
    )

    assert best.offset == 1
    assert any("required_margin=10%" in message for message in trace_messages)
    assert any("accepted shifted offset=1" in message for message in debug_messages)


def test_get_overhangs_from_kmers_rejects_low_support_noisy_branch(
    tmp_path: Path,
    sequential_pool,
) -> None:
    reads = _junction_reads("ATCGG", ("A", "C", "G", "T")) + ["ATCGATACGTACGT"] * 1
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    inferred = kmers.get_overhangs_from_kmers(
        [fastq],
        20,
        100,
        1,
        "ERROR",
        candidate_offsets=(0,),
    )

    assert inferred.motifs == ("ATCGG",)


def test_get_overhangs_from_kmers_raises_when_too_many_strong_motifs_survive(
    tmp_path: Path,
    sequential_pool,
) -> None:
    reads = []
    for prefix in ("ATCGAT", "ATCGAC", "ATCTAT", "ATCTAC"):
        reads.extend(_junction_reads(prefix, ("A", "C", "G", "T")))
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    with pytest.raises(IPyradError, match="more than 3 strong motifs"):
        kmers.get_overhangs_from_kmers(
            [fastq],
            20,
            200,
            1,
            "ERROR",
            candidate_offsets=(0,),
        )


def test_get_overhangs_from_barcoded_reads_handles_mixed_barcode_lengths(
    tmp_path: Path,
    sequential_pool,
) -> None:
    barcodes_by_length = {
        6: ("ACGTAC", "TGCATG"),
        7: ("GATTACA",),
        8: ("CCTGAATC",),
        9: ("AACCGGTTA",),
    }
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        _reads_from_barcodes(barcodes_by_length, ("ATCGG", "ATCGAT")),
    )

    inferred = kmers.get_overhangs_from_barcoded_reads(
        [fastq],
        barcodes_by_length,
        20,
        1_000,
        1,
        "ERROR",
        label="R1 cutsite motif inference",
    )

    assert inferred.position_mode == "barcode_boundary"
    assert inferred.motifs == ("ATCGAT", "ATCGG")
    assert inferred.accepted_reads == 160
    assert inferred.skipped_no_match_reads == 0
    assert inferred.skipped_ambiguous_reads == 0
    assert inferred.boundary_supports == (
        (6, 0, 64),
        (7, 0, 32),
        (8, 0, 32),
        (9, 0, 32),
    )


def test_get_overhangs_from_barcoded_reads_respects_boundary_slack_zero(
    tmp_path: Path,
    sequential_pool,
) -> None:
    barcodes_by_length = {5: ("TAGAG",), 6: ("TTAGAG",)}
    reads = ["TTAGAGTGCAG" + tail + "ACGTACGT" for tail in ("A", "C", "G", "T") * 2]
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    inferred = kmers.get_overhangs_from_barcoded_reads(
        [fastq],
        barcodes_by_length,
        20,
        100,
        1,
        "ERROR",
        label="R1 cutsite motif inference",
        max_barcode_boundary_slack=0,
    )

    assert inferred.motifs == ("TGCAG",)
    assert inferred.skipped_ambiguous_reads == 0
    assert inferred.boundary_supports == ((6, 0, 8),)


def test_get_overhangs_from_barcoded_reads_skips_non_exact_barcode_matches(
    tmp_path: Path,
    sequential_pool,
) -> None:
    barcodes_by_length = {6: ("ACGTAC",)}
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        _reads_from_barcodes(barcodes_by_length, ("ATCGG",), include_mismatch=True),
    )

    inferred = kmers.get_overhangs_from_barcoded_reads(
        [fastq],
        barcodes_by_length,
        20,
        1_000,
        1,
        "ERROR",
        label="R1 cutsite motif inference",
    )

    assert inferred.motifs == ("ATCGG",)
    assert inferred.accepted_reads == 16
    assert inferred.skipped_no_match_reads == 1


def test_get_overhangs_from_barcoded_reads_recovers_motif_from_multi_boundary_reads(
    tmp_path: Path,
    sequential_pool,
) -> None:
    barcodes_by_length = {4: ("ACGT",), 9: ("ACGTATCGG",)}
    reads = []
    for tail in ("A", "C", "G", "T") * 5:
        reads.append("ACGT" + "ATCGG" + tail + "ACGTACGT")
        reads.append("ACGTATCGG" + "ATCGG" + tail + "ACGTACGT")
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", reads)

    inferred = kmers.get_overhangs_from_barcoded_reads(
        [fastq],
        barcodes_by_length,
        24,
        1_000,
        1,
        "ERROR",
        label="R1 cutsite motif inference",
    )

    assert inferred.motifs == ("ATCGG",)
    assert inferred.accepted_reads == 40
    assert inferred.skipped_ambiguous_reads == 40
    assert inferred.boundary_supports == ((4, 0, 40),)


def test_get_overhangs_from_barcoded_reads_emits_debug_and_trace_logs(
    tmp_path: Path,
    sequential_pool,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    barcodes_by_length = {6: ("ACGTAC",)}
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        _reads_from_barcodes(barcodes_by_length, ("ATCGG", "ATCGAT")),
    )
    debug_messages: list[str] = []
    trace_messages: list[str] = []
    stub_logger = type(
        "LoggerStub",
        (),
        {
            "debug": staticmethod(lambda *args: debug_messages.append(args[0].format(*args[1:]))),
            "trace": staticmethod(lambda *args: trace_messages.append(args[0].format(*args[1:]))),
        },
    )
    monkeypatch.setattr(kmers, "logger", stub_logger)

    inferred = kmers.get_overhangs_from_barcoded_reads(
        [fastq],
        barcodes_by_length,
        20,
        1_000,
        1,
        "TRACE",
        label="R1 cutsite motif inference",
    )

    assert inferred.motifs == ("ATCGAT", "ATCGG")
    assert any(
        "R1 cutsite motif inference: sampled 32 reads, matched barcode boundaries in 32, evaluated 1 boundary classes"
        in message
        for message in debug_messages
    )
    assert any(
        "R1 cutsite motif inference summary: motifs [ATCGAT, ATCGG], support 32, no-barcode-match 0, multiple-boundary-match 0"
        in message
        for message in debug_messages
    )
    assert any("boundary=6+0" in message for message in trace_messages)
    assert any("retained motif ATCGG" in message for message in trace_messages)
    assert any("selected_motifs=('ATCGAT', 'ATCGG')" in message for message in trace_messages)


def test_validate_named_motif_list_dedupes_and_preserves_order() -> None:
    motifs = kmers.validate_named_motif_list("ATCGG,ATCGAT,ATCGG", "R1 cutsite motif")

    assert motifs == ("ATCGG", "ATCGAT")
