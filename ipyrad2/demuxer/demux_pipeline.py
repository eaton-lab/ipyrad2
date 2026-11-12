#!/usr/bin/env python
"""Bounded in-memory demux transport, persistent writers, and supervision."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import gzip
import multiprocessing as mp
from pathlib import Path
import queue
import shutil
import signal
import time
import traceback

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline, run_with_pool
from ..utils.logger import setup_loguru_worker
from .match import DemuxRunConfig, build_matcher, get_demux_mode_label
from .sample_names import final_output_sample_name

DEMUX_SPOOL_DIRNAME = ".ipyrad_demux_spool"
PROGRESS_INTERVAL_SECONDS = 5.0


@dataclass(frozen=True)
class DemuxBatch:
    """Immutable paired payloads: (sample ID, R1 bytes, optional R2 bytes)."""

    records: tuple[tuple[int, bytes, bytes | None], ...]


def _demux_spool_dir(outdir):
    """Legacy artifact path, recognized for recovery but never created."""
    return outdir / DEMUX_SPOOL_DIRNAME


def _reader_writer_counts(n_inputs, cores, pigz=False):
    if cores == 1:
        return 0, 0
    readers = min(n_inputs, 2, max(1, cores // 2))
    writers = min(2, cores - readers)
    if n_inputs == 1 and not pigz:
        writers = 1
    return readers, writers


def _partition_inputs(items, parts):
    return [items[idx::parts] for idx in range(parts)]


def _build_writer_route_map(
    barcodes_to_names, merge_technical_replicates, writer_count
):
    names = sorted(
        {
            final_output_sample_name(name, merge_technical_replicates)
            for name in barcodes_to_names.values()
        }
    )
    return {name: idx % writer_count for idx, name in enumerate(names)}


def _handle_limit(requested):
    limit = min(128, requested)
    try:
        import resource

        soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft != resource.RLIM_INFINITY:
            # Include descriptors already in use, in addition to runtime reserve.
            try:
                used = len(list(Path("/proc/self/fd").iterdir()))
            except OSError:
                used = 0
            limit = min(limit, int(soft) - used - 32)
    except ImportError:
        pass
    if limit < 1:
        raise IPyradError(
            "File descriptor limit is too low for demux (32 descriptors reserved)."
        )
    return limit


class BufferedDemuxWriter:
    """Own sample files with bounded pending bytes and a bounded LRU handle cache."""

    def __init__(
        self,
        outdir,
        names,
        pigz=False,
        flush_bytes=1024 * 1024,
        buffer_bytes=16 * 1024 * 1024,
        max_handles=128,
    ):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.names = names
        self.pigz = pigz
        self.flush_bytes = flush_bytes
        self.buffer_bytes = buffer_bytes
        self.max_handles = _handle_limit(max_handles)
        self.buffers = {}
        self.handles = OrderedDict()
        self.pending_bytes = 0
        self.peak_buffer_bytes = 0
        self.peak_handles = 0

    def _flush(self, key):
        buffer = self.buffers.pop(key)
        if not buffer:
            return
        handle = self.handles.pop(key, None)
        if handle is None:
            if len(self.handles) >= self.max_handles:
                _, old = self.handles.popitem(last=False)
                old.close()  # Compression/trailer errors are fatal.
            sample_id, mate = key
            path = _writer_output_path(
                self.outdir, self.names[sample_id], mate, self.pigz
            )
            handle = _open_writer_handle(path, self.pigz)
        self.handles[key] = handle
        self.peak_handles = max(self.peak_handles, len(self.handles))
        handle.write(buffer)
        self.pending_bytes -= len(buffer)

    def add(self, records):
        for sample_id, r1, r2 in records:
            for mate, data in ((1, r1), (2, r2)):
                if data is None:
                    continue
                key = (sample_id, mate)
                self.buffers.setdefault(key, bytearray()).extend(data)
                self.pending_bytes += len(data)
                self.peak_buffer_bytes = max(self.peak_buffer_bytes, self.pending_bytes)
                if len(self.buffers[key]) >= self.flush_bytes:
                    self._flush(key)
                while self.pending_bytes >= self.buffer_bytes and self.buffers:
                    self._flush(
                        max(self.buffers, key=lambda item: len(self.buffers[item]))
                    )

    def finish(self):
        for key in list(self.buffers):
            self._flush(key)
        first_error = None
        while self.handles:
            _, handle = self.handles.popitem(last=False)
            try:
                handle.close()
            except Exception as exc:
                first_error = first_error or exc
        if first_error is not None:
            raise first_error

    def abort(self):
        self.buffers.clear()
        self.pending_bytes = 0
        while self.handles:
            _, handle = self.handles.popitem(last=False)
            try:
                handle.close()
            except Exception as exc:
                logger.warning("demux cleanup: could not close output: {}", exc)


def _writer_output_path(outdir, sample_name, mate, pigz):
    suffix = ".fastq" if pigz else ".fastq.gz"
    return outdir / f"{sample_name}_R{mate}{suffix}"


def _open_writer_handle(path, pigz):
    return path.open("ab") if pigz else gzip.open(path, "ab")


def _queue_get_nowait(qobj):
    if hasattr(qobj, "get_nowait"):
        return qobj.get_nowait()
    if not qobj._reader.poll():
        raise queue.Empty
    return qobj.get()


def _raise_queued_demux_error(qobj):
    try:
        role, worker_id, exc_type, message, tb = _queue_get_nowait(qobj)
    except queue.Empty:
        return
    raise RuntimeError(
        f"demux {role} worker {worker_id} failed with {exc_type}: {message}\n{tb}"
    )


def _raise_on_bad_exit(procs, role):
    for proc in procs:
        if proc.exitcode not in (None, 0):
            raise RuntimeError(
                f"demux {role} worker '{proc.name}' exited with code {proc.exitcode}"
            )


def _put_with_timeout(qobj, payload, timeout, error_queue=None, cancel=None):
    while cancel is None or not cancel.is_set():
        try:
            qobj.put(payload, timeout=timeout)
            return
        except queue.Full:
            # Only parent-side callers may consume error messages.
            if error_queue is not None:
                _raise_queued_demux_error(error_queue)
    raise RuntimeError("demux cancelled")


def _put_with_supervision(qobj, payload, error_queue, writer_procs, timeout):
    while True:
        _raise_queued_demux_error(error_queue)
        _raise_on_bad_exit(writer_procs, "writer")
        try:
            qobj.put(payload, timeout=timeout)
            return
        except queue.Full:
            pass


def _matcher_stats(matcher):
    return (
        matcher.barcode_misses,
        matcher.barcode_hits,
        matcher.sample_hits,
        matcher.barcode_boundary_ambiguities,
        matcher.suspected_barcode_summary(),
    )


def _batches(matcher, config, sample_ids):
    buffers = {}
    nreads = nbytes = 0
    for name, r1, r2 in matcher.iter_output_records():
        sample_id = sample_ids[name]
        pair = buffers.setdefault(sample_id, [bytearray(), None])
        pair[0].extend(r1)
        if r2 is not None:
            if pair[1] is None:
                pair[1] = bytearray()
            pair[1].extend(r2)
        nreads += 1
        nbytes += len(r1) + (len(r2) if r2 is not None else 0)
        if nreads >= config.chunksize or nbytes >= config.batch_bytes:
            yield _freeze_buffers(buffers)
            nreads = nbytes = 0
    if buffers:
        yield _freeze_buffers(buffers)


def _freeze_buffers(buffers):
    records = []
    while buffers:
        sample_id, (r1, r2) = buffers.popitem()
        records.append((sample_id, bytes(r1), bytes(r2) if r2 is not None else None))
    return tuple(records)


def _report_worker_error(channel, role, worker_id, exc):
    try:
        channel.send(
            (
                "error",
                role,
                worker_id,
                type(exc).__name__,
                str(exc),
                traceback.format_exc(),
            )
        )
    except (OSError, EOFError):
        pass  # Parent will also detect a nonzero process exit.


def _reader_process(
    reader_id,
    assignments,
    config,
    queues,
    sample_ids,
    routes,
    status,
    progress,
    cancel,
):
    setup_loguru_worker(config.log_level)
    raw_total = matched_total = 0
    try:
        for fname, fastqs in assignments:
            matcher = build_matcher(fastqs, config, workers=1)

            def report(raw, matched):
                progress[2 * reader_id] = raw_total + raw
                progress[2 * reader_id + 1] = matched_total + matched
                if cancel.is_set():
                    raise RuntimeError("demux cancelled")

            matcher.progress_callback = report
            for batch in _batches(matcher, config, sample_ids):
                groups = [[] for _ in queues]
                for record in batch:
                    groups[routes[record[0]]].append(record)
                for idx, group in enumerate(groups):
                    if group:
                        _put_with_timeout(
                            queues[idx],
                            DemuxBatch(tuple(group)),
                            config.queue_put_timeout,
                            cancel=cancel,
                        )
            report(matcher.reads_seen, matcher.matched_seen)
            raw_total += matcher.reads_seen
            matched_total += matcher.matched_seen
            status.send(("stats", fname, _matcher_stats(matcher)))
        # Flush local feeder threads before acknowledging producer completion.
        for qobj in queues:
            qobj.close()
            qobj.join_thread()
        status.send(("done", "reader", reader_id))
    except BaseException as exc:
        for qobj in queues:
            qobj.cancel_join_thread()
        _report_worker_error(status, "reader", reader_id, exc)
        raise


def _writer_process(writer_id, qobj, config, names, write_dir, status):
    setup_loguru_worker(config.log_level)
    writer = None
    try:
        writer = BufferedDemuxWriter(
            write_dir,
            names,
            config.pigz,
            config.writer_flush_bytes,
            config.writer_buffer_bytes,
            config.writer_max_handles,
        )
        while True:
            batch = qobj.get()
            if batch is None:
                break
            writer.add(batch.records)
        writer.finish()
        status.send(("done", "writer", writer_id))
    except BaseException as exc:
        _report_worker_error(status, "writer", writer_id, exc)
        if writer is not None:
            writer.abort()
        raise


def _terminate_processes(procs):
    for proc in procs:
        if proc.pid is not None and proc.is_alive():
            proc.terminate()
    deadline = time.monotonic() + 1.0
    for proc in procs:
        if proc.pid is not None:
            proc.join(max(0, deadline - time.monotonic()))
    for proc in procs:
        if proc.pid is not None and proc.is_alive():
            proc.kill()
    for proc in procs:
        if proc.pid is not None:
            proc.join(timeout=0.2)


def _compress_temp_fastq(path, threads, outpath=None):
    outpath = outpath or path.with_suffix(path.suffix + ".gz")
    run_pipeline([["pigz", "-p", str(threads), "-c", str(path)]], outfile=outpath)
    path.unlink()


def _compress_demux_outputs_with_pigz(temp_fastqs, cores, log_level, outdir=None):
    paths = sorted(path for path in temp_fastqs if path.exists())
    if not paths:
        return
    workers = min(2, cores, len(paths))
    # Readers and writers have exited. Leave room for pigz I/O and job overhead.
    threads = max(1, cores // workers - 4)
    jobs = {
        path.name: (
            _compress_temp_fastq,
            {
                "path": path,
                "threads": threads,
                "outpath": outdir / (path.name + ".gz") if outdir is not None else None,
            },
        )
        for path in paths
    }
    logger.info(
        "compressing {} demultiplexed FASTQ files with up to {} pigz jobs "
        "and {} compression threads per job",
        len(paths), workers, threads,
    )
    run_with_pool(jobs, log_level=log_level, max_workers=workers, msg="Pigz FASTQs")


def run_demux_pipeline(filenames_to_fastqs, config: DemuxRunConfig, cores):
    """Run serial or supervised bounded parallel demux into private staging."""
    if config.pigz and shutil.which("pigz") is None:
        raise IPyradError("--pigz requires pigz to be installed and on PATH.")
    items = list(filenames_to_fastqs.items())
    readers_count, writers_count = _reader_writer_counts(len(items), cores, config.pigz)
    names = sorted(
        {
            final_output_sample_name(name, config.merge_technical_replicates)
            for name in config.barcodes_to_names.values()
        }
    )
    sample_ids = {name: idx for idx, name in enumerate(names)}
    writers_count = min(writers_count, len(names))
    write_dir = (
        config.tmpdir if config.pigz and config.tmpdir is not None else config.outdir
    )
    write_dir.mkdir(parents=True, exist_ok=True)
    logger.info(
        "demultiplexing on {} with {} reader(s) and {} writer(s) using {} compression",
        get_demux_mode_label(config),
        readers_count or 1,
        writers_count,
        "pigz" if config.pigz else "python",
    )
    for fname, fastqs in items:
        logger.info(
            "processing {} {}",
            fname,
            tuple(path.name if path else "" for path in fastqs),
        )
    file_stats = {}
    if cores == 1:
        writer = BufferedDemuxWriter(
            write_dir,
            names,
            config.pigz,
            config.writer_flush_bytes,
            config.writer_buffer_bytes,
            config.writer_max_handles,
        )
        last_log = time.monotonic()
        raw_total = matched_total = 0

        def report(raw, matched):
            nonlocal last_log
            if time.monotonic() - last_log >= PROGRESS_INTERVAL_SECONDS:
                logger.info(
                    "demux progress: raw_reads={} matched_reads={} phase=reading/writing",
                    raw_total + raw,
                    matched_total + matched,
                )
                last_log = time.monotonic()

        try:
            for fname, fastqs in items:
                matcher = build_matcher(fastqs, config, workers=1)
                matcher.progress_callback = report
                for batch in _batches(matcher, config, sample_ids):
                    writer.add(batch)
                raw_total += matcher.reads_seen
                matched_total += matcher.matched_seen
                file_stats[fname] = _matcher_stats(matcher)
            writer.finish()
        except BaseException:
            writer.abort()
            raise
    else:
        ctx = mp.get_context("spawn")
        queues = [ctx.Queue(maxsize=2) for _ in range(writers_count)]
        # Each worker has one exclusive send endpoint. Closing the parent's copy
        # makes a killed sender produce EOF even midway through a framed message;
        # a shared Queue/SimpleQueue can otherwise leave recv() blocked forever.
        channels = {
            (role, idx): ctx.Pipe(duplex=False)
            for role, count in (("reader", readers_count), ("writer", writers_count))
            for idx in range(count)
        }
        live_channels = {key: pair[0] for key, pair in channels.items()}
        progress = ctx.RawArray("q", 2 * readers_count)
        cancel = ctx.Event()
        routes = [idx % writers_count for idx in range(len(names))]
        writers = [
            ctx.Process(
                name=f"demux-writer-{idx}",
                target=_writer_process,
                args=(
                    idx,
                    queues[idx],
                    config,
                    names,
                    write_dir,
                    channels[("writer", idx)][1],
                ),
            )
            for idx in range(writers_count)
        ]
        readers = [
            ctx.Process(
                name=f"demux-reader-{idx}",
                target=_reader_process,
                args=(
                    idx,
                    group,
                    config,
                    queues,
                    sample_ids,
                    routes,
                    channels[("reader", idx)][1],
                    progress,
                    cancel,
                ),
            )
            for idx, group in enumerate(_partition_inputs(items, readers_count))
        ]
        procs = writers + readers
        done = set()
        last_log = time.monotonic()

        def drain_status():
            for identity, channel in list(live_channels.items()):
                while channel.poll():
                    try:
                        message = channel.recv()
                    except EOFError:
                        del live_channels[identity]
                        break
                    except OSError as exc:
                        raise RuntimeError(
                            f"demux {identity[0]} worker {identity[1]} returned an incomplete message"
                        ) from exc
                    if message[0] == "error":
                        _, role, worker_id, exc_type, text, tb = message
                        raise RuntimeError(
                            f"demux {role} worker {worker_id} failed with {exc_type}: {text}\n{tb}"
                        )
                    kind, key, value = message
                    if kind == "stats":
                        file_stats[key] = value
                    else:
                        done.add((key, value))

        def supervise(phase):
            nonlocal last_log
            drain_status()
            for role, workers in (("reader", readers), ("writer", writers)):
                for idx, proc in enumerate(workers):
                    code = proc.exitcode
                    if code is not None and (code != 0 or (role, idx) not in done):
                        # Completion was sent before exit, so it must now be readable.
                        drain_status()
                        if code != 0 or (role, idx) not in done:
                            raise RuntimeError(
                                f"demux {role} worker {idx} exited unexpectedly with code {code}"
                            )
            if time.monotonic() - last_log >= PROGRESS_INTERVAL_SECONDS:
                values = list(progress[:])
                logger.info(
                    "demux progress: raw_reads={} matched_reads={} phase={} readers={} writers={}",
                    sum(values[::2]),
                    sum(values[1::2]),
                    phase,
                    sum(p.is_alive() for p in readers),
                    sum(p.is_alive() for p in writers),
                )
                last_log = time.monotonic()

        aborted = False
        try:
            for proc in procs:
                proc.start()
            for _, sender in channels.values():
                sender.close()
            while any(proc.is_alive() for proc in readers):
                supervise("reading/writing")
                time.sleep(0.05)
            supervise("writing")
            if len(file_stats) != len(items):
                raise RuntimeError("demux missing input statistics")
            for qobj in queues:
                while True:
                    supervise("writing")
                    try:
                        qobj.put(None, timeout=config.queue_put_timeout)
                        break
                    except queue.Full:
                        pass
            while any(proc.is_alive() for proc in writers):
                supervise("writing")
                time.sleep(0.05)
            supervise("writing")
            for proc in procs:
                proc.join()
        except BaseException:
            aborted = True
            previous = signal.getsignal(signal.SIGINT)
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            try:
                cancel.set()
                for qobj in queues:
                    qobj.cancel_join_thread()
                _terminate_processes(procs)
            finally:
                signal.signal(signal.SIGINT, previous)
            raise
        finally:
            for qobj in queues:
                qobj.close()
                if not aborted:
                    qobj.join_thread()
            for receiver, sender in channels.values():
                receiver.close()
                sender.close()
            for proc in procs:
                proc.close()
    if config.pigz:
        paths = [
            _writer_output_path(write_dir, name, mate, True)
            for name in names
            for mate in (1, 2)
            if _writer_output_path(write_dir, name, mate, True).exists()
        ]
        _compress_demux_outputs_with_pigz(
            paths, cores, config.log_level, outdir=config.outdir
        )
    raw_total = sum(
        sum(stats[idx].values()) for stats in file_stats.values() for idx in (0, 1, 3)
    )
    matched_total = sum(sum(stats[1].values()) for stats in file_stats.values())
    logger.info(
        "demux complete: raw_reads={} matched_reads={}", raw_total, matched_total
    )
    return file_stats
