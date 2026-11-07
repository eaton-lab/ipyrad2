"""Correctness, bounded resources, and failure recovery for hardened demux."""

from collections import Counter
from dataclasses import replace
import gzip
import multiprocessing as mp
import os
from pathlib import Path
import queue
import shutil
import time

import pytest

from ipyrad2.demuxer.demux import Demux
from ipyrad2.demuxer import demux_pipeline as pipeline
from ipyrad2.demuxer.match import (
    DemuxRunConfig,
    build_matcher,
    record_bounded_observation,
    DIAGNOSTIC_OVERFLOW,
)
from ipyrad2.demuxer.transaction import DemuxTransaction, pending_runs, output_lock
from ipyrad2.utils.exceptions import IPyradError


def record(name="r0", seq="ACGTATCGGAAAA", mate=None):
    suffix = f"/{mate}" if mate else ""
    return f'@{name}{suffix}\n{seq}\n+\n{"I"*len(seq)}\n'.encode()


def write(path, data):
    if path.suffix == ".gz":
        with gzip.open(path, "wb") as out:
            out.write(data)
    else:
        path.write_bytes(data)
    return path


def config(tmp_path, **kw):
    values = dict(
        barcodes_to_names={b"ACGT": "a", b"TGCA": "b"},
        barcode_lengths1=(4,),
        barcode_lengths2=(),
        cuts1=[b"ATCGG"],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        chunksize=2,
        max_reads=None,
        i7=False,
        log_level="WARNING",
    )
    values.update(kw)
    return DemuxRunConfig(**values)


def tool(tmp_path, **kw):
    raw = tmp_path / "lane.fastq.gz"
    if not raw.exists():
        write(raw, record())
    bars = tmp_path / "bars.tsv"
    bars.write_text("a ACGT\nb TGCA\n")
    values = dict(
        fastqs=[raw],
        barcodes=bars,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=1,
        chunksize=2,
        merge_technical_replicates=False,
        outdir=tmp_path / "out",
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=None,
        max_reads_kmer=100,
        log_level="WARNING",
    )
    values.update(kw)
    return Demux(**values)


def read_records(path):
    with gzip.open(path, "rb") as source:
        lines = source.readlines()
    assert len(lines) % 4 == 0
    return [tuple(lines[i : i + 4]) for i in range(0, len(lines), 4)]


@pytest.mark.parametrize(
    "data,message",
    [
        (b"@r\nACGT\n+\n", "truncated"),
        (b"r\nACGT\n+\nIIII\n", "@ header"),
        (b"@r\nACGT\nx\nIIII\n", "separator"),
        (b"@r\nACGT\n+\nIII\n", "lengths differ"),
    ],
)
def test_invalid_fastq_has_path_and_record_number(tmp_path, data, message):
    raw = write(tmp_path / "raw.fq", record() + data)
    matcher = build_matcher((raw, None), config(tmp_path), 1)
    with pytest.raises(IPyradError, match=f"raw.fq, record 2: .*{message}"):
        list(matcher.iter_output_records())


@pytest.mark.parametrize(
    "r2,message",
    [
        (b"", "unequal mate counts"),
        (record("different", mate=2), "mismatched mate IDs"),
        (record(mate=1), "incorrect mate annotation"),
    ],
)
def test_pair_validation(tmp_path, r2, message):
    r1 = write(tmp_path / "R1.fq.gz", record(mate=1))
    r2 = write(tmp_path / "R2.fq", r2)
    with pytest.raises(IPyradError, match=message):
        list(build_matcher((r1, r2), config(tmp_path), 1).iter_output_records())


def test_mixed_compression_crlf_missing_final_newline_and_prefix_limit(tmp_path):
    r1 = write(
        tmp_path / "R1.fq.gz", record(mate=1).replace(b"\n", b"\r\n") + b"invalid tail"
    )
    r2 = write(tmp_path / "R2.fq", record(mate=2).rstrip(b"\n"))
    matcher = build_matcher((r1, r2), config(tmp_path, max_reads=1), 1)
    result = list(matcher.iter_output_records())
    assert len(result) == 1
    assert result[0][1].endswith(b"\n") and result[0][2].endswith(b"\n")
    assert matcher.reads_seen == matcher.matched_seen == 1


def test_corrupt_gzip_and_empty_i7(tmp_path):
    raw = write(tmp_path / "raw.gz", record())
    raw.write_bytes(raw.read_bytes()[:-6])
    with pytest.raises(IPyradError, match="raw.gz, record"):
        list(build_matcher((raw, None), config(tmp_path), 1).iter_output_records())
    write(raw, b"")
    matcher = build_matcher((raw, None), config(tmp_path, i7=True), 1)
    assert list(matcher.iter_output_records()) == []
    assert matcher.reads_seen == 0


def test_i7_counts_each_record_once(tmp_path):
    raw = write(tmp_path / "raw.gz", record("r0 1:N:0:ACGT"))
    matcher = build_matcher((raw, None), config(tmp_path, i7=True), 1)
    assert len(list(matcher.iter_output_records())) == 1
    assert matcher.reads_seen == matcher.matched_seen == 1


def test_exact_totals_with_bounded_observation_detail():
    counter = {}
    for i in range(10000):
        record_bounded_observation(counter, str(i).encode())
    record_bounded_observation(counter, b"0")
    assert len(counter) == 4097
    assert sum(counter.values()) == 10001
    assert counter[b"0"] == 2
    assert counter[DIAGNOSTIC_OVERFLOW] == 5904


def test_writer_buffers_handles_and_gzip_members_are_bounded(tmp_path):
    names = [f"s{i}" for i in range(96)]
    writer = pipeline.BufferedDemuxWriter(
        tmp_path, names, flush_bytes=200, buffer_bytes=500, max_handles=2
    )
    expected = {name: bytearray() for name in names}
    for cycle in range(4):
        for i, name in enumerate(names):
            data = record(f"{name}-{cycle}")
            writer.add(((i, data, None),))
            expected[name].extend(data)
            assert writer.pending_bytes < 500
            assert len(writer.handles) <= 2
    # A single oversized record is flushed immediately, never split.
    data = record("oversized", "A" * 3000)
    writer.add(((0, data, None),))
    expected[names[0]].extend(data)
    assert writer.pending_bytes < 500
    writer.finish()
    for name in names:
        assert (
            gzip.decompress((tmp_path / f"{name}_R1.fastq.gz").read_bytes())
            == expected[name]
        )
    assert writer.peak_handles == 2
    assert writer.peak_buffer_bytes <= 500 + len(data)


def test_descriptor_limit_reserves_runtime_handles(monkeypatch):
    resource = pytest.importorskip("resource")
    monkeypatch.setattr(resource, "getrlimit", lambda _: (64, 64))
    assert 1 <= pipeline._handle_limit(128) <= 32
    monkeypatch.setattr(resource, "getrlimit", lambda _: (16, 16))
    with pytest.raises(IPyradError, match="descriptor"):
        pipeline._handle_limit(128)


@pytest.mark.parametrize("failure", ["write", "close", "eviction"])
def test_writer_errors_are_not_suppressed(tmp_path, monkeypatch, failure):
    class Broken:
        def write(self, data):
            if failure == "write":
                raise OSError("disk full")
            return len(data)

        def close(self):
            if failure != "write":
                raise OSError("trailer failed")
            return None

    monkeypatch.setattr(pipeline, "_open_writer_handle", lambda *a: Broken())
    writer = pipeline.BufferedDemuxWriter(
        tmp_path, ["a", "b"], flush_bytes=1, max_handles=1
    )
    with pytest.raises(OSError):
        writer.add(((0, record(), None),))
        if failure == "eviction":
            writer.add(((1, record(), None),))
        writer.finish()
    writer.abort()


@pytest.mark.parametrize("cores", [1, 2, 4, 8])
def test_parallel_multiset_equivalence_and_pair_synchrony(tmp_path, cores):
    inputs = {}
    for lane in range(3):
        reads = []
        for index in range(12):
            barcode = "ACGT" if index % 3 else "TGCA"
            name = f"lane{lane}-read{index//2}"  # duplicate IDs and records are intentional
            r1 = record(name, barcode + "ATCGGAAAA", 1)
            r2 = record(name, "GATTACA", 2)
            reads.append((r1, r2))
        r1 = write(tmp_path / f"lane{lane}_R1.gz", b"".join(pair[0] for pair in reads))
        r2 = write(tmp_path / f"lane{lane}_R2.gz", b"".join(pair[1] for pair in reads))
        inputs[str(lane)] = (r1, r2)
    base_config = config(
        tmp_path, chunksize=3, batch_bytes=100, writer_buffer_bytes=150
    )
    pipeline.run_demux_pipeline(
        inputs, replace(base_config, outdir=tmp_path / "serial"), 1
    )
    actual_stats = pipeline.run_demux_pipeline(inputs, base_config, cores)
    for name in ("a", "b"):
        serial = zip(
            read_records(tmp_path / f"serial/{name}_R1.fastq.gz"),
            read_records(tmp_path / f"serial/{name}_R2.fastq.gz"),
        )
        actual = zip(
            read_records(tmp_path / f"out/{name}_R1.fastq.gz"),
            read_records(tmp_path / f"out/{name}_R2.fastq.gz"),
        )
        assert Counter(serial) == Counter(actual)
    assert sum(sum(stats[1].values()) for stats in actual_stats.values()) == 36
    assert not pipeline._demux_spool_dir(tmp_path / "out").exists()


def test_force_does_not_delete_outputs_on_constructor_or_processing_failure(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    previous = out / "a_R1.fastq.gz"
    previous.write_bytes(b"original")
    write(tmp_path / "lane.fastq.gz", b"@bad\nACGT\n+\n")
    demux = tool(tmp_path, force=True)
    assert previous.read_bytes() == b"original"
    with pytest.raises(IPyradError, match="truncated"):
        demux.run()
    assert previous.read_bytes() == b"original"
    assert not pending_runs(out)
    assert not list(out.glob("ipyrad_demux_stats*"))


def test_report_failure_does_not_publish_fastqs(tmp_path, monkeypatch):
    demux = tool(tmp_path)

    def fail():
        raise OSError("report failure")

    monkeypatch.setattr(demux, "_write_stats", fail)
    with pytest.raises(OSError, match="report failure"):
        demux.run()
    assert not list(demux.outdir.glob("*.fastq*"))
    assert not pending_runs(demux.outdir)


def test_success_preserves_unrelated_files_and_numbers_reports(tmp_path):
    first = tool(tmp_path)
    first.run()
    (first.outdir / "unrelated.txt").write_text("keep")
    old_report = (first.outdir / "ipyrad_demux_stats_0.txt").read_bytes()
    tool(tmp_path, force=True).run()
    assert (first.outdir / "ipyrad_demux_stats_0.txt").read_bytes() == old_report
    assert (first.outdir / "ipyrad_demux_stats_1.txt").exists()
    assert (first.outdir / "unrelated.txt").read_text() == "keep"
    assert not pending_runs(first.outdir)


@pytest.mark.parametrize("failure_position", [1, 2, 3, 4])
def test_publication_rollback_is_restartable(tmp_path, monkeypatch, failure_position):
    import ipyrad2.demuxer.transaction as transaction

    (tmp_path / "a").write_text("old-a")
    (tmp_path / "b").write_text("old-b")
    tx = DemuxTransaction(tmp_path)
    (tx.outputs / "a").write_text("new-a")
    (tx.outputs / "b").write_text("new-b")
    real_replace = os.replace
    calls = 0

    def fail(source, dest):
        nonlocal calls
        if Path(source).name != "manifest.tmp":
            calls += 1
            if calls == failure_position:
                raise OSError("rename failed")
        return real_replace(source, dest)

    monkeypatch.setattr(transaction.os, "replace", fail)
    with pytest.raises(OSError, match="rename failed"):
        tx.publish(
            [(tx.outputs / "a", tmp_path / "a"), (tx.outputs / "b", tmp_path / "b")]
        )
    monkeypatch.setattr(transaction.os, "replace", real_replace)
    # Simulate restart from the on-disk manifest, not the in-memory object.
    recovered = DemuxTransaction.load(tx.root)
    recovered.rollback()
    recovered.rollback()
    recovered.cleanup()
    assert (tmp_path / "a").read_text() == "old-a"
    assert (tmp_path / "b").read_text() == "old-b"


def test_stale_run_requires_force_and_recovers(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    tx = DemuxTransaction(out)
    with pytest.raises(IPyradError, match="Interrupted demux"):
        tool(tmp_path)
    tool(tmp_path, force=True).run()
    assert not tx.root.exists()


def test_output_lock_rejects_concurrent_runs(tmp_path):
    with output_lock(tmp_path):
        with pytest.raises(IPyradError, match="Another demux"):
            with output_lock(tmp_path):
                pass


@pytest.mark.skipif(shutil.which("pigz") is None, reason="pigz not installed")
@pytest.mark.parametrize("cores", [1, 4])
def test_pigz_external_scratch_equivalence_and_cleanup(tmp_path, cores):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "keep").write_text("keep")
    tool(tmp_path).run()
    baseline = read_records(tmp_path / "out/a_R1.fastq.gz")
    tool(tmp_path, force=True, pigz=True, tmpdir=scratch, cores=cores).run()
    assert read_records(tmp_path / "out/a_R1.fastq.gz") == baseline
    assert list(scratch.iterdir()) == [scratch / "keep"]
    assert not list((tmp_path / "out").glob("*.fastq"))


def test_tmpdir_requires_pigz(tmp_path):
    with pytest.raises(IPyradError, match="tmpdir requires pigz"):
        tool(tmp_path, tmpdir=tmp_path)


def _exit_zero_writer(*args):
    os._exit(0)


def _error_writer(*args):
    raise OSError("injected writer crash")


def _unrelated_child(event):
    event.wait(30)


@pytest.mark.parametrize("target", [_exit_zero_writer, _error_writer])
def test_worker_failure_aborts_only_owned_processes(tmp_path, monkeypatch, target):
    ctx = mp.get_context("spawn")
    event = ctx.Event()
    other = ctx.Process(target=_unrelated_child, args=(event,))
    other.start()
    before = time.monotonic()
    monkeypatch.setattr(pipeline, "_writer_process", target)
    raw = write(tmp_path / "lane.gz", record() * 100)
    try:
        with pytest.raises(RuntimeError, match="writer worker"):
            pipeline.run_demux_pipeline(
                {"lane": (raw, None)}, config(tmp_path, chunksize=1), 2
            )
        assert other.is_alive()
        assert time.monotonic() - before < 10
    finally:
        event.set()
        other.join(5)
        if other.is_alive():
            other.kill()
            other.join()
        other.close()


def test_cancelled_producer_does_not_wait_forever():
    class Full:
        def put(self, *args, **kw):
            raise queue.Full

    class Cancelled:
        def is_set(self):
            return True

    with pytest.raises(RuntimeError, match="cancelled"):
        pipeline._put_with_timeout(Full(), None, 0.001, cancel=Cancelled())


def _exit_zero_reader(*args):
    os._exit(0)


def _delayed_error_writer(*args):
    time.sleep(0.7)
    raise OSError("writer failed under backpressure")


@pytest.mark.parametrize(
    "role,target", [("reader", _exit_zero_reader), ("writer", _delayed_error_writer)]
)
def test_premature_reader_exit_and_full_queue_failure(
    tmp_path, monkeypatch, role, target
):
    monkeypatch.setattr(pipeline, f"_{role}_process", target)
    raw = write(tmp_path / "lane.gz", record() * 1000)
    start = time.monotonic()
    with pytest.raises(RuntimeError, match=f"{role} worker"):
        pipeline.run_demux_pipeline(
            {"lane": (raw, None)}, config(tmp_path, chunksize=1), 2
        )
    assert time.monotonic() - start < 10


def test_real_worker_write_failure_and_error_metadata(tmp_path):
    raw = write(tmp_path / "lane.gz", record())
    conf = config(tmp_path)
    conf.outdir.mkdir()
    (conf.outdir / "a_R1.fastq.gz").mkdir()
    with pytest.raises(RuntimeError, match="writer worker .*IsADirectoryError"):
        pipeline.run_demux_pipeline({"lane": (raw, None)}, conf, 2)


def test_heartbeat_does_not_depend_on_chunksize(tmp_path, monkeypatch):
    messages = []

    class Logger:
        def info(self, message, *args):
            messages.append(message.format(*args))

    monkeypatch.setattr(pipeline, "logger", Logger())
    monkeypatch.setattr(pipeline, "PROGRESS_INTERVAL_SECONDS", 0.01)
    raw = write(tmp_path / "lane.gz", record())
    pipeline.run_demux_pipeline(
        {"lane": (raw, None)}, config(tmp_path, chunksize=10**7), 2
    )
    assert any("demux progress:" in msg and "phase=" in msg for msg in messages)
    assert any("demux complete: raw_reads=1 matched_reads=1" in msg for msg in messages)


@pytest.mark.skipif(shutil.which("pigz") is None, reason="pigz not installed")
def test_compression_failure_removes_only_owned_scratch(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "keep").write_text("unrelated")

    def fail(*args, **kw):
        raise OSError("compression failure")

    monkeypatch.setattr(pipeline, "_compress_demux_outputs_with_pigz", fail)
    with pytest.raises(OSError, match="compression failure"):
        tool(tmp_path, pigz=True, tmpdir=scratch).run()
    assert list(scratch.iterdir()) == [scratch / "keep"]
    assert not list((tmp_path / "out").glob("*.fastq*"))
    assert not pending_runs(tmp_path / "out")


def test_recovery_of_committed_run_preserves_new_outputs(tmp_path):
    tx = DemuxTransaction(tmp_path)
    (tx.outputs / "a").write_text("new")
    tx.publish([(tx.outputs / "a", tmp_path / "a")])
    recovered = DemuxTransaction.load(tx.root)
    recovered.rollback()  # committed runs only need cleanup
    recovered.cleanup()
    assert (tmp_path / "a").read_text() == "new"


def test_recovery_rejects_unsafe_manifest(tmp_path):
    import json

    tx = DemuxTransaction(tmp_path)
    state = tx.state.copy()
    state.update(phase="publishing", entries=[{"source": "../outside", "target": "a"}])
    (tx.root / "manifest.json").write_text(json.dumps(state))
    with pytest.raises(IPyradError, match="unsafe artifact name"):
        DemuxTransaction.load(tx.root)


def test_reader_batch_byte_limit_with_large_chunksize(tmp_path):
    raw = write(tmp_path / "lane.gz", record() * 100)
    conf = config(tmp_path, chunksize=10**7, batch_bytes=100)
    matcher = build_matcher((raw, None), conf, 1)
    batches = list(pipeline._batches(matcher, conf, {"a": 0, "b": 1}))
    assert len(batches) > 1
    assert all(
        sum(len(r1) + len(r2 or b"") for _, r1, r2 in batch) <= 100 + len(record())
        for batch in batches
    )
    assert sum(r1.count(b"@r0\n") for batch in batches for _, r1, _ in batch) == 100


def _partial_message_reader(*args):
    import struct

    channel = args[6]
    # Connection framing declares a larger payload than the worker writes.
    os.write(channel.fileno(), struct.pack("!i", 100000))
    os.write(channel.fileno(), b"partial")
    os._exit(1)


@pytest.mark.skipif(os.name != "posix", reason="uses POSIX pipe framing")
def test_killed_worker_mid_status_message_does_not_hang(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline, "_reader_process", _partial_message_reader)
    raw = write(tmp_path / "lane.gz", record())
    start = time.monotonic()
    with pytest.raises(RuntimeError, match="reader worker.*incomplete message"):
        pipeline.run_demux_pipeline({"lane": (raw, None)}, config(tmp_path), 2)
    assert time.monotonic() - start < 10
