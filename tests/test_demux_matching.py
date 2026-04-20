import gzip
import queue
import shutil
from collections import Counter
from pathlib import Path

import pytest

import ipyrad2.demuxer.demux as demux_module
import ipyrad2.demuxer.demux_pipeline as pipeline_module
import ipyrad2.demuxer.demux_report as report_module
import ipyrad2.demuxer.match as match_module
from ipyrad2.demuxer.demux import Demux
from ipyrad2.demuxer.demux_pipeline import _demux_spool_dir
from ipyrad2.demuxer.demux_pipeline import _put_with_supervision
from ipyrad2.demuxer.demux_pipeline import _put_with_timeout
from ipyrad2.demuxer.demux_pipeline import _reader_writer_counts
from ipyrad2.demuxer.demux_pipeline import run_demux_pipeline
from ipyrad2.demuxer.match import BarMatchingSingleInline
from ipyrad2.demuxer.match import DemuxRunConfig
from ipyrad2.demuxer.match import cut_matcher
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.kmers import InferredJunctionSet


def _write_fastq(path: Path, reads: list[str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for idx, read in enumerate(reads):
            out.write(f"@r{idx}\n{read}\n+\n{'I' * len(read)}\n")
    return path


def _write_fastq_records(path: Path, records: list[tuple[str, str, str]]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for header, seq, qual in records:
            out.write(f"{header}\n{seq}\n+\n{qual}\n")
    return path


def _read_fastq_sequences(path: Path) -> list[str]:
    seqs: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as infile:
        while True:
            header = infile.readline()
            if not header:
                break
            seqs.append(infile.readline().strip())
            infile.readline()
            infile.readline()
    return seqs


def _junction_set(
    motifs: tuple[str, ...],
    *,
    offset: int = 0,
    counts: tuple[int, ...] | None = None,
    runner_up_offset_support: int = 0,
    candidate_offsets: tuple[int, ...] = (4, 5),
) -> InferredJunctionSet:
    counts = counts or tuple(100 for _ in motifs)
    return InferredJunctionSet(
        motifs=motifs,
        motif_counts=counts,
        offset=offset,
        total_support=sum(counts),
        runner_up_offset_support=runner_up_offset_support,
        candidate_offsets=candidate_offsets,
    )


def test_format_logged_motif_set_hides_boundary_support_counts() -> None:
    junction = InferredJunctionSet(
        motifs=("ATCGG", "ATCGAT"),
        motif_counts=(10, 6),
        offset=0,
        total_support=16,
        runner_up_offset_support=0,
        candidate_offsets=(0,),
        position_mode="barcode_boundary",
        boundary_supports=((6, 0, 10), (7, 0, 6)),
    )

    assert report_module.format_logged_motif_set(junction) == (
        "[ATCGG, ATCGAT] inferred from barcode boundaries"
    )


class _DemuxLoggerStub:
    def __init__(self) -> None:
        self.messages: list[str] = []

    def _record(self, message: str, *args) -> None:
        self.messages.append(message.format(*args) if args else message)

    def info(self, message: str, *args) -> None:
        self._record(message, *args)

    def warning(self, message: str, *args) -> None:
        self._record(message, *args)

    def debug(self, message: str, *args) -> None:
        self._record(message, *args)

    def error(self, message: str, *args) -> None:
        self._record(message, *args)


def test_cut_matcher_accepts_expected_and_plus_one_offsets() -> None:
    assert cut_matcher(b"ACGTATCGGAAAA", (4,), [b"ATCGG"]) == (b"ACGT", 4)
    assert cut_matcher(b"NACGTATCGGAAAA", (4,), [b"ATCGG"]) == (b"ACGT", 5)


def test_cut_matcher_rejects_internal_occurrence() -> None:
    assert cut_matcher(b"ACGTGGGGATCGGAAAA", (4,), [b"ATCGG"]) is None


def test_reader_writer_counts_cap_single_input_and_multi_input_writers() -> None:
    assert _reader_writer_counts(1, 4, pigz=False) == (1, 1)
    assert _reader_writer_counts(1, 4, pigz=True) == (1, 2)
    assert _reader_writer_counts(2, 4, pigz=False) == (2, 2)
    assert _reader_writer_counts(3, 8, pigz=False) == (3, 2)


def test_put_with_supervision_raises_if_writer_is_dead() -> None:
    class _AlwaysFullQueue:
        def put(self, payload, timeout=None):
            raise queue.Full

    class _DeadWriter:
        name = "demux-writer-0"
        exitcode = 1

    with pytest.raises(RuntimeError, match="demux writer worker"):
        _put_with_supervision(
            _AlwaysFullQueue(),
            None,
            queue.Queue(),
            [_DeadWriter()],
            timeout=0.0,
        )


def test_put_with_timeout_raises_queued_worker_error() -> None:
    class _AlwaysFullQueue:
        def put(self, payload, timeout=None):
            raise queue.Full

    errors = queue.Queue()
    errors.put(("writer", 0, "RuntimeError", "boom", "traceback"))

    with pytest.raises(RuntimeError, match="demux writer worker 0 failed"):
        _put_with_timeout(
            _AlwaysFullQueue(),
            None,
            timeout=0.0,
            error_queue=errors,
        )


def test_bounded_barcode_counter_retains_frequent_keys_with_bounded_state() -> None:
    counter = match_module.BoundedBarcodeCounter(capacity=4)
    for _ in range(100):
        counter.add(("R1", b"AAAA"))
    for _ in range(30):
        counter.add(("R1", b"CCCC"))
    for idx in range(100):
        counter.add(("R1", f"{idx:04d}".encode()))

    summary = counter.summary()

    assert len(summary) <= 4
    assert ("R1", b"AAAA") in summary
    for estimate, error in summary.values():
        assert estimate >= error
    assert len(counter._heap) <= max(counter.capacity * 4, counter.capacity + 128)


def test_bar_matching_progress_callback_reports_raw_and_matched_reads(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC"],
    )
    reports: list[tuple[int, int]] = []
    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={b"ACGT": "sample1", b"TGCA": "sample2"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="WARNING",
        workers=1,
        chunksize=10,
        max_reads=10,
        progress_callback=lambda raw_reads, matched_reads: reports.append((raw_reads, matched_reads)),
        progress_interval_reads=1,
    )

    list(matcher.iter_output_records())
    matcher._maybe_report_progress(force=True)

    assert reports[0] == (1, 0)
    assert reports[-1] == (2, 2)


def test_bar_matching_single_inline_records_suspected_unknown_barcode(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["TTTTATCGGAAAA"])
    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={b"ACGT": "sample1"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="WARNING",
        workers=1,
        chunksize=10,
        max_reads=10,
    )

    assert list(matcher.iter_output_records()) == []
    assert matcher.barcode_misses == {b"XXX": 1}
    assert matcher.suspected_barcode_summary() == {("R1", b"TTTT"): (1, 0)}


def test_bar_matching_boundary_ambiguous_reads_are_not_output_by_default(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["TTAGAGTGCAGAAAA"])
    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={
            b"TTAGAG": "multiplex-DC15-4",
            b"TAGAG": "occidentale-CO18-01",
        },
        barcode_lengths1=(5, 6),
        barcode_lengths2=(),
        cuts1=[b"TGCAG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="WARNING",
        workers=1,
        chunksize=10,
        max_reads=10,
    )

    assert list(matcher.iter_output_records()) == []
    assert matcher.sample_hits == {}
    assert matcher.barcode_boundary_ambiguities == {
        b"boundary_ambiguous:multiplex-DC15-4:TTAGAG;occidentale-CO18-01:TAGAG": 1,
    }


def test_bar_matching_boundary_slack_zero_recovers_position_zero_candidate(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["TTAGAGTGCAGAAAA"])
    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={
            b"TTAGAG": "multiplex-DC15-4",
            b"TAGAG": "occidentale-CO18-01",
        },
        barcode_lengths1=(5, 6),
        barcode_lengths2=(),
        cuts1=[b"TGCAG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="WARNING",
        workers=1,
        chunksize=10,
        max_reads=10,
        barcode_boundary_slack=0,
    )

    records = list(matcher.iter_output_records())

    assert len(records) == 1
    assert records[0][0] == "multiplex-DC15-4"
    assert b"\nTGCAGAAAA\n" in records[0][1]
    assert matcher.sample_hits == {"multiplex-DC15-4": 1}
    assert matcher.barcode_boundary_ambiguities == {}


def test_demux_rejects_monomorphic_barcode1(tmp_path: Path) -> None:
    r1 = _write_fastq(tmp_path / "reads.fastq.gz", ["AAAAAAATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 AAAAAA\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="cannot be monomorphic or all N"):
        Demux(
            fastqs=[r1],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="ERROR",
        )


def test_demux_rejects_monomorphic_barcode2(tmp_path: Path) -> None:
    r1 = _write_fastq(tmp_path / "reads_R1.fastq.gz", ["ACGTATCGGAAAA"])
    r2 = _write_fastq(tmp_path / "reads_R2.fastq.gz", ["TGCAATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT AAAAAA\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="cannot be monomorphic or all N"):
        Demux(
            fastqs=[r1, r2],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2="CGATCC",
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="ERROR",
        )


def test_demux_rejects_monomorphic_user_overhang(tmp_path: Path) -> None:
    r1 = _write_fastq(tmp_path / "reads.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="cannot be monomorphic or all N"):
        Demux(
            fastqs=[r1],
            barcodes=barcodes,
            cutsite_1="AAAAAA",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="ERROR",
        )


def test_demux_creates_missing_outdir_parents(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outdir = tmp_path / "nested" / "parents" / "out"

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=outdir,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert outdir.is_dir()
    assert _read_fastq_sequences(outdir / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]
    assert (outdir / "ipyrad_demux_stats_0.txt").exists()


def test_demux_rejects_outdir_when_path_is_existing_file(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outpath = tmp_path / "out"
    outpath.write_text("not a directory", encoding="utf-8")

    with pytest.raises(IPyradError, match="exists and is not a directory"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=outpath,
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_rejects_barcode_glob_matching_multiple_files(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    (tmp_path / "barcodes_a.tsv").write_text("sample1 ACGT\n", encoding="utf-8")
    (tmp_path / "barcodes_b.tsv").write_text("sample2 TGCA\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="matches multiple files"):
        Demux(
            fastqs=[raw],
            barcodes=tmp_path / "barcodes_*.tsv",
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_rejects_incomplete_barcode_rows(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="contains incomplete rows"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_rejects_sanitized_sample_name_collisions(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample+1 ACGT\nsample?1 TGCA\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="Sanitized sample names would collide"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_raises_on_existing_current_run_artifact_without_force(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outdir = tmp_path / "out"
    outdir.mkdir()
    _write_fastq(outdir / "sample1_R1.fastq.gz", ["OLDREAD"])

    with pytest.raises(
        IPyradError,
        match="One or more files matching the expected output names exist in outdir. Use --force to overwrite.",
    ):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=outdir,
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_force_preserves_unrelated_files_and_warns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outdir = tmp_path / "out"
    outdir.mkdir()
    _write_fastq(outdir / "sample1_R1.fastq.gz", ["OLDREAD"])
    _write_fastq(outdir / "other.fastq.gz", ["UNRELATED"])
    (outdir / "ipyrad_demux_stats_0.txt").write_text("old stats\n", encoding="utf-8")
    (outdir / "notes.txt").write_text("keep me\n", encoding="utf-8")

    logger = _DemuxLoggerStub()
    monkeypatch.setattr(demux_module, "logger", logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=outdir,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        force=True,
    )
    tool.run()

    assert _read_fastq_sequences(outdir / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]
    assert _read_fastq_sequences(outdir / "other.fastq.gz") == ["UNRELATED"]
    assert (outdir / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8") == "old stats\n"
    assert (outdir / "ipyrad_demux_stats_1.txt").exists()
    assert (outdir / "notes.txt").read_text(encoding="utf-8") == "keep me\n"
    assert any("demux stats files are present" in message for message in logger.messages)
    assert any("FASTQ.gz files are present" in message for message in logger.messages)
    assert not any("sample1_R1.fastq.gz" in message and "will not be overwritten" in message for message in logger.messages)


def test_demux_warns_about_partial_outputs_on_generic_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    logger = _DemuxLoggerStub()
    monkeypatch.setattr(demux_module, "logger", logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    monkeypatch.setattr(tool, "_demultiplex", lambda: (_ for _ in ()).throw(RuntimeError("boom")))

    with pytest.raises(RuntimeError, match="boom"):
        tool.run()

    assert any("may contain partial files" in message for message in logger.messages)


def test_demux_warns_for_multi_motif_auto_inference_message(monkeypatch: pytest.MonkeyPatch) -> None:
    messages = []
    stub_logger = type(
        "LoggerStub",
        (),
        {
            "warning": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
        },
    )
    monkeypatch.setattr(report_module, "logger", stub_logger)

    report_module.warn_multi_motif_inference(
        "R1",
        _junction_set(("ATCGG", "ATCGAT"), counts=(80, 20)),
        100_000,
    )

    assert "multiple motifs" in messages[0]
    assert "--max-reads-kmer" in messages[0]
    assert "3RAD" in messages[0]
    assert "using all detected motifs" in messages[0]


def test_demux_pipeline_accepts_multiple_manual_r1_motifs(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "ACGTATCGATCCCC", "ACGTATCGGGGGG"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG,ATCGAT",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGATCCCC", "ATCGGGGGG"]
    )


def test_demux_manual_r1_motif_still_runs_matching_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    calls = []
    logger = _DemuxLoggerStub()

    def fake_barcode_aware_inference(*args, **kwargs):
        calls.append((args, kwargs))
        return _junction_set(("ATCGG",), counts=(10,), candidate_offsets=(0,))

    monkeypatch.setattr(demux_module, "get_overhangs_from_barcoded_reads", fake_barcode_aware_inference)
    monkeypatch.setattr(demux_module, "logger", logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=False,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert len(calls) == 1
    assert tool._re1_motifs == ("ATCGG",)
    assert tool._re1_source == "manual"
    assert any("match detected motifs" in message for message in logger.messages)
    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]


def test_demux_manual_r1_motif_overrides_different_detected_motif(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    logger = _DemuxLoggerStub()

    monkeypatch.setattr(
        demux_module,
        "get_overhangs_from_barcoded_reads",
        lambda *args, **kwargs: _junction_set(("GGGGG",), counts=(10,), candidate_offsets=(0,)),
    )
    monkeypatch.setattr(demux_module, "logger", logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=False,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")

    assert any("do not match" in message and "overrule" in message for message in logger.messages)
    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]
    assert "detected" in stats_text
    assert "selected" in stats_text
    assert "GGGGG" in stats_text
    assert "ATCGG" in stats_text
    assert "manual motifs override detected motifs" in stats_text


def test_demux_manual_r1_motif_warns_and_proceeds_when_inference_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    logger = _DemuxLoggerStub()

    def fake_barcode_aware_inference(*args, **kwargs):
        raise IPyradError("no motifs found")

    monkeypatch.setattr(demux_module, "get_overhangs_from_barcoded_reads", fake_barcode_aware_inference)
    monkeypatch.setattr(demux_module, "logger", logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=False,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")

    assert any("inference failed" in message and "Using user-defined motifs" in message for message in logger.messages)
    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]
    assert "inference failed; using manual motifs" in stats_text


def test_demux_disable_infer_with_manual_motif_does_not_run_kmer_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")

    def fail_if_called(*args, **kwargs):
        raise AssertionError("kmer inference should not run")

    monkeypatch.setattr(demux_module, "get_overhangs_from_barcoded_reads", fail_if_called)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]


def test_demux_without_manual_motif_still_fails_when_inference_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")

    def fake_barcode_aware_inference(*args, **kwargs):
        raise IPyradError("no motifs found")

    monkeypatch.setattr(demux_module, "get_overhangs_from_barcoded_reads", fake_barcode_aware_inference)

    with pytest.raises(IPyradError, match="no motifs found"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1=None,
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=False,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_writes_multi_motif_inference_stats_and_demuxes_both_classes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "ACGTATCGATCCCC", "ACGTATCGGGGGG"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")

    monkeypatch.setattr(
        demux_module,
        "get_overhangs_from_barcoded_reads",
        lambda *args, **kwargs: _junction_set(("ATCGG", "ATCGAT"), counts=(80, 20)),
    )

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=False,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGATCCCC", "ATCGGGGGG"]
    )
    assert "# Restriction motif inference" in stats_text
    assert "ATCGG" in stats_text
    assert "ATCGAT" in stats_text
    assert "support_fraction" in stats_text


def test_demux_direct_path_allows_max_reads_none(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=None,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGTTTT"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz") == ["ATCGGCCCC"]


def test_demux_rejects_unrecoverable_single_inline_boundary_collision(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "ACGTATCGGCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("short ACGT\nlong ACGTATCGG\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="unrecoverable R1 barcode-boundary collisions"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_rejects_ambiguous_barcode_mismatch_expansion(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "ACGAATCGGCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 ACGA\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="ambiguous barcode candidates"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=1,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_combinatorial_collision_can_be_resolved_by_joint_pairing(tmp_path: Path) -> None:
    r1 = _write_fastq(
        tmp_path / "sample_R1.fastq.gz",
        [
            "ACGTATCGGAAAA",
            "ACGTATCGGCCCC",
            "ACGTATCGGATCGGAAAA",
            "ACGTATCGGATCGGCCCC",
        ],
    )
    r2 = _write_fastq(
        tmp_path / "sample_R2.fastq.gz",
        [
            "ACGTCGATCAAAA",
            "ACGTCGATCCCC",
            "TGCACGATCAAAA",
            "TGCACGATCCCC",
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("short ACGT ACGT\nlong ACGTATCGG TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[r1, r2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2="CGATC",
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert tool._sample_stats["short"] == 2
    assert tool._sample_stats["long"] == 2
    assert tool._barcode_boundary_collisions
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "short_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGCCCC"]
    )
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "long_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGCCCC"]
    )


def test_demux_auto_inference_uses_exact_barcode_matches_but_demux_respects_max_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        [
            "ACGTACATCGGAAAA",
            "GATTACAATCGATCCCC",
            "ACGTTCATCGGTTTT",
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGTAC\nsample2 GATTACA\n", encoding="utf-8")
    calls = []

    def fake_barcode_aware_inference(
        fastqs,
        barcodes_by_length,
        max_len,
        max_reads,
        workers,
        log_level,
        *,
        label="demux",
        max_barcode_boundary_slack=1,
    ):
        calls.append((barcodes_by_length, max_barcode_boundary_slack))
        return _junction_set(
            ("ATCGG", "ATCGAT"),
            counts=(80, 20),
            candidate_offsets=(0,),
        )

    monkeypatch.setattr(demux_module, "get_overhangs_from_barcoded_reads", fake_barcode_aware_inference)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=1,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=False,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        barcode_boundary_slack=0,
    )
    tool.run()

    assert calls == [({6: ("ACGTAC",), 7: ("GATTACA",)}, 0)]
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGTTTT"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz") == ["ATCGATCCCC"]


def test_demux_stats_report_barcode_boundary_classes_for_auto_inference(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTACATCGGAAAA", "GATTACAATCGATCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGTAC\nsample2 GATTACA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG,ATCGAT",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool._re1_inference = InferredJunctionSet(
        motifs=("ATCGG", "ATCGAT"),
        motif_counts=(10, 6),
        offset=0,
        total_support=16,
        runner_up_offset_support=0,
        candidate_offsets=(0,),
        position_mode="barcode_boundary",
        sampled_reads=16,
        accepted_reads=16,
        skipped_no_match_reads=0,
        skipped_ambiguous_reads=0,
        boundary_supports=((6, 0, 10), (7, 0, 6)),
    )
    tool._re1_source = "auto"
    tool._re1_motifs = tool._re1_inference.motifs
    tool.run()

    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")

    assert "position" in stats_text
    assert "6+0:10, 7+0:6" in stats_text
    assert "not assigned or written to any sample output file" in stats_text


def test_demux_stats_report_suspected_missing_inline_barcode(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["TTTTATCGGAAAA"] * 60 + ["ACGTATCGGCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")

    assert "# Suspected missing barcode statistics" in stats_text
    assert "bounded-memory estimated counts" in stats_text
    assert "TTTT" in stats_text
    assert "nearest_expected_mismatches" in stats_text


def test_demux_pipeline_single_inline_writes_expected_outputs(tmp_path: Path) -> None:
    raw1 = _write_fastq(
        tmp_path / "lane1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    raw2 = _write_fastq(
        tmp_path / "lane2.fastq.gz",
        ["TGCAATCGGGGGG", "ACGTATCGGAAAA"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw1, raw2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGTTTT", "ATCGGAAAA"]
    )
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz")) == Counter(
        ["ATCGGCCCC", "ATCGGGGGG"]
    )
    assert tool._sample_stats["sample1"] == 3
    assert tool._sample_stats["sample2"] == 2


def test_demux_accepted_mismatch_barcode_is_not_suspected(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["AGGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=1,
        cores=1,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]
    assert report_module.aggregate_suspected_barcode_stats(tool._file_stats) == {}


def test_demux_pipeline_aggregates_suspected_unknown_barcodes(tmp_path: Path) -> None:
    raw1 = _write_fastq(tmp_path / "lane1.fastq.gz", ["TTTTATCGGAAAA"] * 30)
    raw2 = _write_fastq(tmp_path / "lane2.fastq.gz", ["TTTTATCGGCCCC"] * 30)
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw1, raw2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=5,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    suspected = report_module.aggregate_suspected_barcode_stats(tool._file_stats)
    assert suspected == {("R1", b"TTTT"): (60, 0)}


def test_demux_pipeline_paired_end_writes_r1_and_r2_outputs(tmp_path: Path) -> None:
    r1 = _write_fastq(
        tmp_path / "lane_R1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    r2 = _write_fastq(
        tmp_path / "lane_R2.fastq.gz",
        ["GGGGAAAA", "CCCCGGGG", "TTTTAAAA"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[r1, r2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2="CGATCC",
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGTTTT"]
    )
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R2.fastq.gz")) == Counter(
        ["GGGGAAAA", "TTTTAAAA"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R2.fastq.gz") == ["CCCCGGGG"]


def test_demux_combinatorial_records_suspected_unknown_r2_barcode(tmp_path: Path) -> None:
    r1 = _write_fastq(tmp_path / "lane_R1.fastq.gz", ["ACGTATCGGAAAA"])
    r2 = _write_fastq(tmp_path / "lane_R2.fastq.gz", ["CCCCCGATCCAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text(
        "sample1 ACGT TTTA\nsample2 TGCA GGGA\n",
        encoding="utf-8",
    )

    tool = Demux(
        fastqs=[r1, r2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2="CGATCC",
        max_mismatch=0,
        cores=1,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    suspected = report_module.aggregate_suspected_barcode_stats(tool._file_stats)
    assert suspected == {("R2", b"CCCC"): (1, 0)}


def test_demux_pipeline_i7_writes_expected_outputs(tmp_path: Path) -> None:
    raw = _write_fastq_records(
        tmp_path / "lane.fastq.gz",
        [
            ("@r0 1:N:0:ACGT+AAAA", "AAAACCCC", "IIIIIIII"),
            ("@r1 1:N:0:TGCA+AAAA", "GGGGTTTT", "IIIIIIII"),
            ("@r2 1:N:0:ACGT+AAAA", "CCCCAAAA", "IIIIIIII"),
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=True,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["AAAACCCC", "CCCCAAAA"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz") == ["GGGGTTTT"]


def test_demux_i7_records_and_reports_suspected_unknown_index(tmp_path: Path) -> None:
    raw = _write_fastq_records(
        tmp_path / "lane.fastq.gz",
        [("@known 1:N:0:ACGT+AAAA", "AAAACCCC", "IIIIIIII")]
        + [
            (f"@unknown{i} 1:N:0:CCCC+AAAA", "GGGGTTTT", "IIIIIIII")
            for i in range(60)
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=True,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    suspected = report_module.aggregate_suspected_barcode_stats(tool._file_stats)
    assert suspected == {("i7", b"CCCC"): (60, 0)}
    stats_text = (tmp_path / "out" / "ipyrad_demux_stats_0.txt").read_text(encoding="utf-8")
    assert "i7" in stats_text
    assert "CCCC" in stats_text


def test_demux_pipeline_i7_paired_end_writes_r1_and_r2_outputs(tmp_path: Path) -> None:
    r1 = _write_fastq_records(
        tmp_path / "lane_R1.fastq.gz",
        [
            ("@r0 1:N:0:ACGT+AAAA", "AAAACCCC", "IIIIIIII"),
            ("@r1 1:N:0:TGCA+AAAA", "GGGGTTTT", "IIIIIIII"),
            ("@r2 1:N:0:ACGT+AAAA", "CCCCAAAA", "IIIIIIII"),
        ],
    )
    r2 = _write_fastq_records(
        tmp_path / "lane_R2.fastq.gz",
        [
            ("@r0 2:N:0:ACGT+AAAA", "TTTTGGGG", "IIIIIIII"),
            ("@r1 2:N:0:TGCA+AAAA", "CCCCGGGG", "IIIIIIII"),
            ("@r2 2:N:0:ACGT+AAAA", "AAAATTTT", "IIIIIIII"),
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[r1, r2],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=True,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["AAAACCCC", "CCCCAAAA"]
    )
    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R2.fastq.gz")) == Counter(
        ["TTTTGGGG", "AAAATTTT"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R2.fastq.gz") == ["CCCCGGGG"]


def test_demux_i7_ignores_barcode2_column_and_warns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    raw = _write_fastq_records(
        tmp_path / "lane.fastq.gz",
        [
            ("@r0 1:N:0:ACGT+AAAA", "AAAACCCC", "IIIIIIII"),
            ("@r1 1:N:0:TGCA+AAAA", "GGGGTTTT", "IIIIIIII"),
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text(
        "sample1 ACGT AAAAAA\nsample2 TGCA TTTTTT\n",
        encoding="utf-8",
    )
    messages: list[str] = []

    def _record(*args) -> None:
        template = str(args[0]) if args else ""
        values = args[1:]
        messages.append(template.format(*values) if values else template)

    stub_logger = type(
        "LoggerStub",
        (),
        {
            "info": staticmethod(_record),
            "warning": staticmethod(_record),
            "debug": staticmethod(_record),
            "error": staticmethod(_record),
        },
    )
    monkeypatch.setattr(demux_module, "logger", stub_logger)

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=True,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert tool._barcode_lengths2 == ()
    assert tool._names_to_barcodes == {
        "sample1": ("ACGT", ""),
        "sample2": ("TGCA", ""),
    }
    assert any("Ignoring barcode2 and any extra barcode columns" in message for message in messages)
    assert _read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz") == ["AAAACCCC"]
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz") == ["GGGGTTTT"]


def test_demux_pipeline_i7_respects_max_mismatch(tmp_path: Path) -> None:
    raw = _write_fastq_records(
        tmp_path / "lane.fastq.gz",
        [
            ("@r0 1:N:0:AGGT+AAAA", "AAAACCCC", "IIIIIIII"),
            ("@r1 1:N:0:TGCA+AAAA", "GGGGTTTT", "IIIIIIII"),
            ("@r2 1:N:0:AGGT+AAAA", "CCCCAAAA", "IIIIIIII"),
        ],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1=None,
        cutsite_2=None,
        max_mismatch=1,
        cores=1,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=True,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["AAAACCCC", "CCCCAAAA"]
    )
    assert _read_fastq_sequences(tmp_path / "out" / "sample2_R1.fastq.gz") == ["GGGGTTTT"]


def test_demux_pipeline_merges_technical_replicates(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample1 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=True,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGCCCC"]
    )
    assert tool._sample_stats["sample1"] == 2
    assert tool._technical_replicates["sample1"] == [
        "sample1-technical-replicate-0",
        "sample1-technical-replicate-1",
    ]


def test_demux_rejects_final_output_name_collisions_after_technical_replicate_merging(
    tmp_path: Path,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA", "TGCAATCGGCCCC"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text(
        "sample1 ACGT\nsample1-technical-replicate-0 TGCA\n",
        encoding="utf-8",
    )

    with pytest.raises(IPyradError, match="Final demux output sample names would collide"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=1,
            merge_technical_replicates=True,
            outdir=tmp_path / "out",
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_serial_path_merges_technical_replicates(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample1 TGCA\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=1,
        merge_technical_replicates=True,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    tool.run()

    assert Counter(_read_fastq_sequences(tmp_path / "out" / "sample1_R1.fastq.gz")) == Counter(
        ["ATCGGAAAA", "ATCGGCCCC"]
    )
    assert tool._sample_stats["sample1"] == 2
    assert tool._technical_replicates["sample1"] == [
        "sample1-technical-replicate-0",
        "sample1-technical-replicate-1",
    ]


def test_demux_pipeline_logs_single_input_writer_cap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    r1 = _write_fastq(
        tmp_path / "lane_R1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    r2 = _write_fastq(
        tmp_path / "lane_R2.fastq.gz",
        ["GGGGAAAA", "CCCCGGGG", "TTTTAAAA"],
    )
    messages: list[str] = []
    stub_logger = type(
        "LoggerStub",
        (),
        {
            "info": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
            "warning": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
            "error": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
        },
    )
    monkeypatch.setattr(pipeline_module, "logger", stub_logger)

    config = DemuxRunConfig(
        barcodes_to_names={b"ACGT": "sample1", b"TGCA": "sample2"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        chunksize=1,
        max_reads=100,
        i7=False,
        log_level="WARNING",
    )
    run_demux_pipeline({"lane": (r1, r2)}, config, cores=4)

    assert any("1 reader(s) and 1 writer(s)" in message for message in messages)


def test_demux_pipeline_cleans_up_spool_dir(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    config = DemuxRunConfig(
        barcodes_to_names={b"ACGT": "sample1", b"TGCA": "sample2"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        chunksize=1,
        max_reads=None,
        i7=False,
        log_level="WARNING",
    )

    run_demux_pipeline({"lane": (raw, None)}, config, cores=3)

    assert not _demux_spool_dir(tmp_path / "out").exists()


def test_demux_rejects_existing_spool_dir_without_force(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outdir = tmp_path / "out"
    outdir.mkdir()
    spool_dir = _demux_spool_dir(outdir)
    spool_dir.mkdir()
    (spool_dir / "stale.part").write_text("stale\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="Existing demux spool directory"):
        Demux(
            fastqs=[raw],
            barcodes=barcodes,
            cutsite_1="ATCGG",
            cutsite_2=None,
            max_mismatch=0,
            cores=1,
            chunksize=10,
            merge_technical_replicates=False,
            outdir=outdir,
            i7=False,
            disable_infer_cutsite_motifs=True,
            max_reads=100,
            max_reads_kmer=100,
            log_level="WARNING",
        )


def test_demux_force_removes_existing_spool_dir(tmp_path: Path) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\n", encoding="utf-8")
    outdir = tmp_path / "out"
    outdir.mkdir()
    spool_dir = _demux_spool_dir(outdir)
    spool_dir.mkdir()
    (spool_dir / "stale.part").write_text("stale\n", encoding="utf-8")

    tool = Demux(
        fastqs=[raw],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=10,
        merge_technical_replicates=False,
        outdir=outdir,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        force=True,
    )

    assert not spool_dir.exists()
    tool.run()
    assert _read_fastq_sequences(outdir / "sample1_R1.fastq.gz") == ["ATCGGAAAA"]


class _LoggerStub:
    def __init__(self) -> None:
        self.messages: list[str] = []

    def info(self, message: str, *args) -> None:
        self.messages.append(message.format(*args) if args else message)


def test_demux_pipeline_does_not_log_progress_before_full_chunk(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA"],
    )
    logger = _LoggerStub()
    monkeypatch.setattr(pipeline_module, "logger", logger)

    config = DemuxRunConfig(
        barcodes_to_names={b"ACGT": "sample1"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        chunksize=2,
        max_reads=1,
        i7=False,
        log_level="WARNING",
    )
    run_demux_pipeline({"lane": (raw, None)}, config, cores=2)

    assert not any(message.startswith("demux progress:") for message in logger.messages)
    assert not any("demux active:" in message for message in logger.messages)


def test_bar_matching_uses_configured_log_level_for_writer_pool(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    seen: dict[str, object] = {}

    class _PoolStub:
        def __init__(self, *, log_level: str, max_workers: int) -> None:
            seen["log_level"] = log_level
            seen["max_workers"] = max_workers

        def iter_results(self, items, max_inflight=None):
            return iter(())

        def abort(self, fast: bool = False) -> None:
            seen["aborted"] = fast

        def close(self, wait: bool = True) -> None:
            seen["close_wait"] = wait

    monkeypatch.setattr(match_module, "_ManagedProcessPool", _PoolStub)

    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={b"ACGT": "sample1"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="DEBUG",
        workers=2,
        chunksize=1,
        max_reads=10,
    )
    matcher.run()

    assert seen["log_level"] == "DEBUG"
    assert seen["max_workers"] == 2
    assert seen["close_wait"] is True


def test_bar_matching_interrupt_cleanup_does_not_print(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw = _write_fastq(tmp_path / "lane.fastq.gz", ["ACGTATCGGAAAA"])
    seen: dict[str, object] = {}

    class _PoolStub:
        def __init__(self, *, log_level: str, max_workers: int) -> None:
            seen["log_level"] = log_level
            seen["max_workers"] = max_workers

        def iter_results(self, items, max_inflight=None):
            raise KeyboardInterrupt

        def abort(self, fast: bool = False) -> None:
            seen["aborted"] = fast

        def close(self, wait: bool = True) -> None:
            seen["close_wait"] = wait

    monkeypatch.setattr(match_module, "_ManagedProcessPool", _PoolStub)
    monkeypatch.setattr(
        "builtins.print",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("print should not be called")),
    )

    matcher = BarMatchingSingleInline(
        fastqs=(raw, None),
        barcodes_to_names={b"ACGT": "sample1"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        log_level="WARNING",
        workers=2,
        chunksize=1,
        max_reads=10,
    )

    with pytest.raises(KeyboardInterrupt):
        matcher.run()

    assert seen["aborted"] is True
    assert seen["close_wait"] is False


def test_demux_pipeline_logs_cumulative_chunk_progress_across_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw1 = _write_fastq(
        tmp_path / "lane1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC"],
    )
    raw2 = _write_fastq(
        tmp_path / "lane2.fastq.gz",
        ["ACGTATCGGGGGG", "TGCAATCGGTTTT"],
    )
    logger = _LoggerStub()
    monkeypatch.setattr(pipeline_module, "logger", logger)

    config = DemuxRunConfig(
        barcodes_to_names={b"ACGT": "sample1", b"TGCA": "sample2"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        chunksize=4,
        max_reads=4,
        i7=False,
        log_level="WARNING",
    )
    run_demux_pipeline(
        {
            "lane1": (raw1, None),
            "lane2": (raw2, None),
        },
        config,
        cores=1,
    )

    progress_messages = [
        message for message in logger.messages
        if message.startswith("demux progress:")
    ]

    assert progress_messages == ["demux progress: raw_reads=4 matched_reads=4"]
    assert not any("demux active:" in message for message in logger.messages)


@pytest.mark.skipif(shutil.which("pigz") is None, reason="pigz not installed")
def test_demux_pipeline_pigz_backend_matches_python_backend(tmp_path: Path) -> None:
    raw1 = _write_fastq(
        tmp_path / "lane1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    raw2 = _write_fastq(
        tmp_path / "lane2.fastq.gz",
        ["TGCAATCGGGGGG", "ACGTATCGGAAAA"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    python_tool = Demux(
        fastqs=[raw1, raw2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "python_out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        pigz=False,
    )
    python_tool.run()

    pigz_tool = Demux(
        fastqs=[raw1, raw2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=tmp_path / "pigz_out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        pigz=True,
    )
    pigz_tool.run()

    for sample in ("sample1", "sample2"):
        assert Counter(_read_fastq_sequences(tmp_path / "python_out" / f"{sample}_R1.fastq.gz")) == Counter(
            _read_fastq_sequences(tmp_path / "pigz_out" / f"{sample}_R1.fastq.gz")
        )


def test_demux_pipeline_pigz_writer_outputs_complete_temp_fastqs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw1 = _write_fastq(
        tmp_path / "lane1.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT"],
    )
    raw2 = _write_fastq(
        tmp_path / "lane2.fastq.gz",
        ["TGCAATCGGGGGG", "ACGTATCGGAAAA"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")
    outdir = tmp_path / "pigz_out"
    outdir.mkdir()
    unrelated_fastq = outdir / "unrelated.fastq"
    unrelated_fastq.write_text("UNRELATED\n", encoding="utf-8")
    captured_paths: list[Path] = []

    def _capture_temp_fastqs(temp_fastqs, *args, **kwargs) -> None:
        captured_paths.extend(temp_fastqs)

    monkeypatch.setattr(pipeline_module, "_compress_demux_outputs_with_pigz", _capture_temp_fastqs)

    pigz_tool = Demux(
        fastqs=[raw1, raw2],
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=outdir,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=100,
        max_reads_kmer=100,
        log_level="WARNING",
        pigz=True,
    )
    pigz_tool.run()

    def _read_plain_fastq_sequences(path: Path) -> list[str]:
        lines = path.read_text(encoding="utf-8").splitlines()
        return [lines[idx] for idx in range(1, len(lines), 4)]

    assert Counter(_read_plain_fastq_sequences(outdir / "sample1_R1.fastq")) == Counter(
        ["ATCGGAAAA", "ATCGGTTTT", "ATCGGAAAA"]
    )
    assert Counter(_read_plain_fastq_sequences(outdir / "sample2_R1.fastq")) == Counter(
        ["ATCGGCCCC", "ATCGGGGGG"]
    )
    assert sorted(path.name for path in captured_paths) == ["sample1_R1.fastq", "sample2_R1.fastq"]
    assert unrelated_fastq.read_text(encoding="utf-8") == "UNRELATED\n"


def test_compress_demux_outputs_with_pigz_uses_compact_progress_label(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    temp_fastqs = [
        tmp_path / "sample1_R1.fastq",
        tmp_path / "sample2_R1.fastq",
    ]
    for path in temp_fastqs:
        path.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def _fake_run_with_pool(jobs, log_level, max_workers, msg):
        captured["job_names"] = sorted(jobs)
        captured["log_level"] = log_level
        captured["max_workers"] = max_workers
        captured["msg"] = msg
        return {}

    monkeypatch.setattr(pipeline_module, "run_with_pool", _fake_run_with_pool)

    pipeline_module._compress_demux_outputs_with_pigz(
        temp_fastqs=temp_fastqs,
        cores=8,
        log_level="INFO",
    )

    assert captured == {
        "job_names": ["sample1_R1.fastq", "sample2_R1.fastq"],
        "log_level": "INFO",
        "max_workers": 2,
        "msg": "Pigz FASTQs",
    }
