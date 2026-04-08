#!/usr/bin/env python

"""Process-based demux pipeline orchestration and cleanup."""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple
import gzip
import multiprocessing as mp
import queue
import shutil
import signal
import time
import traceback
from pathlib import Path
from dataclasses import dataclass

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline, run_with_pool
from ..utils.logger import setup_loguru_worker
from .match import DemuxRunConfig, build_matcher, get_demux_mode_label
from .sample_names import final_output_sample_name


DEMUX_SPOOL_DIRNAME = ".ipyrad_demux_spool"


@dataclass(frozen=True)
class DemuxShardEntry:
    """Metadata describing one sample/mate byte slice within a shard file."""

    sample_name: str
    mate: int
    offset: int
    length: int


@dataclass(frozen=True)
class DemuxShardMessage:
    """Small control message sent from readers to writers for one shard file."""

    path: Path
    size_bytes: int
    entries: Tuple[DemuxShardEntry, ...]


def _reader_writer_counts(
    n_inputs: int,
    cores: int,
    pigz: bool = False,
) -> Tuple[int, int]:
    """Split the demux core budget across readers and writers."""
    if n_inputs == 1:
        reader_workers = 1
        if not pigz:
            writer_workers = 1
        else:
            writer_workers = min(2, max(1, cores - 1))
        return reader_workers, writer_workers
    reader_workers = min(n_inputs, max(1, cores // 2))
    writer_workers = min(2, max(1, cores - reader_workers))
    return reader_workers, writer_workers


def _partition_inputs(
    items: List[Tuple[str, Tuple[Path, Path | None]]],
    parts: int,
) -> List[List[Tuple[str, Tuple[Path, Path | None]]]]:
    """Distribute input file tuples across reader workers."""
    groups: List[List[Tuple[str, Tuple[Path, Path | None]]]] = [[] for _ in range(parts)]
    for idx, item in enumerate(items):
        groups[idx % parts].append(item)
    return [group for group in groups if group]


def _build_writer_route_map(
    barcodes_to_names: Dict[bytes, str],
    merge_technical_replicates: bool,
    writer_count: int,
) -> Dict[str, int]:
    """Return a stable parent-built route map from output sample name to writer."""
    sample_names = sorted({
        final_output_sample_name(name, merge_technical_replicates)
        for name in barcodes_to_names.values()
    })
    return {
        sample_name: idx % writer_count
        for idx, sample_name in enumerate(sample_names)
    }


def _demux_spool_dir(outdir: Path) -> Path:
    """Return the internal shard spool directory for one demux run."""
    return outdir / DEMUX_SPOOL_DIRNAME


def _cancel_queue_join(qobj: Any) -> None:
    """Disable blocking queue-feeder joins during process teardown."""
    try:
        qobj.cancel_join_thread()
    except Exception:
        pass


def _close_queue(qobj: Any) -> None:
    """Close a multiprocessing queue-like object."""
    try:
        qobj.close()
    except Exception:
        pass


def _queue_get_nowait(qobj: Any) -> Any:
    """Return one queued item without blocking for Queue or SimpleQueue objects."""
    getter = getattr(qobj, "get_nowait", None)
    if callable(getter):
        try:
            return getter()
        except (EOFError, OSError, ValueError) as exc:
            raise queue.Empty from exc

    poll = getattr(qobj, "_poll", None)
    if callable(poll):
        try:
            if not poll():
                raise queue.Empty
        except Exception as exc:
            raise queue.Empty from exc
        try:
            return qobj.get()
        except (EOFError, OSError, ValueError) as exc:
            raise queue.Empty from exc

    empty = getattr(qobj, "empty", None)
    if callable(empty):
        try:
            if empty():
                raise queue.Empty
        except Exception as exc:
            raise queue.Empty from exc
    try:
        return qobj.get()
    except (EOFError, OSError, ValueError) as exc:
        raise queue.Empty from exc


def _put_with_retry(
    qobj: Any,
    payload: Any,
    timeout: float,
    *,
    error_queue: Any | None = None,
    worker_procs: List[mp.Process] | None = None,
    worker_role: str = "worker",
    sleep_seconds: float = 0.01,
) -> None:
    """Put an item on a queue while optionally supervising queued worker errors."""
    while True:
        try:
            qobj.put(payload, timeout=timeout)
            return
        except queue.Full:
            if error_queue is not None:
                _raise_queued_demux_error(error_queue)
            if worker_procs is not None:
                _raise_on_bad_exit(worker_procs, worker_role)
            time.sleep(sleep_seconds)


def _put_with_timeout(
    qobj: Any,
    payload: Any,
    timeout: float,
    error_queue: Any | None = None,
) -> None:
    """Put an item on a queue with retry-on-full semantics."""
    _put_with_retry(
        qobj,
        payload,
        timeout,
        error_queue=error_queue,
        sleep_seconds=0.01,
    )


def _put_with_supervision(
    qobj: Any,
    payload: Any,
    error_queue: Any,
    writer_procs: List[mp.Process],
    timeout: float,
) -> None:
    """Put a control message while continuing to supervise writer health."""
    _put_with_retry(
        qobj,
        payload,
        timeout,
        error_queue=error_queue,
        worker_procs=writer_procs,
        worker_role="writer",
        sleep_seconds=0.05,
    )


def _build_shard_message(
    spool_dir: Path,
    reader_id: int,
    writer_idx: int,
    shard_idx: int,
    payload: Dict[str, Tuple[bytes, bytes | None]],
) -> DemuxShardMessage:
    """Write one reader batch shard to disk and return its control message."""
    path = spool_dir / f"reader{reader_id:02d}-writer{writer_idx:02d}-shard{shard_idx:06d}.bin"
    entries: List[DemuxShardEntry] = []
    offset = 0
    with path.open("wb") as out:
        for sample_name in sorted(payload):
            read1_data, read2_data = payload[sample_name]
            out.write(read1_data)
            entries.append(DemuxShardEntry(sample_name, 1, offset, len(read1_data)))
            offset += len(read1_data)
            if read2_data is not None:
                out.write(read2_data)
                entries.append(DemuxShardEntry(sample_name, 2, offset, len(read2_data)))
                offset += len(read2_data)
    return DemuxShardMessage(path=path, size_bytes=offset, entries=tuple(entries))


def _flush_sample_buffers(
    sample_buffers: Dict[str, List[bytearray | None]],
    writer_queues: List[Any],
    writer_route_map: Dict[str, int],
    spool_dir: Path,
    reader_id: int,
    shard_idx: int,
    queue_put_timeout: float,
    error_queue: Any | None = None,
) -> int:
    """Flush one reader batch to the owning writer queues."""
    writer_payloads: Dict[int, Dict[str, Tuple[bytes, bytes | None]]] = {}
    for sample_name, (read1_buf, read2_buf) in sample_buffers.items():
        writer_idx = writer_route_map[sample_name]
        writer_payloads.setdefault(writer_idx, {})[sample_name] = (
            bytes(read1_buf),
            bytes(read2_buf) if read2_buf is not None else None,
        )
    for writer_idx, payload in writer_payloads.items():
        shard = _build_shard_message(spool_dir, reader_id, writer_idx, shard_idx, payload)
        shard_idx += 1
        _put_with_timeout(
            writer_queues[writer_idx],
            shard,
            timeout=queue_put_timeout,
            error_queue=error_queue,
        )
    return shard_idx


def _writer_output_path(
    outdir: Path,
    sample_name: str,
    mate: int,
    pigz: bool,
) -> Path:
    """Return the output path for one sample mate under the selected backend."""
    suffix = ".fastq" if pigz else ".fastq.gz"
    return outdir / f"{sample_name}_R{mate}{suffix}"


def _open_writer_handle(path: Path, pigz: bool):
    """Open a persistent binary writer handle for one sample mate."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if pigz:
        return path.open("ab")
    return gzip.open(path, "ab")


def _close_writer_handles(handles: Dict[Tuple[str, int], Any]) -> None:
    """Close all open demux writer handles."""
    for handle in handles.values():
        try:
            handle.close()
        except Exception:
            pass


def _flush_writer_buffer(
    sample_name: str,
    mate: int,
    buffer: bytearray,
    handles: Dict[Tuple[str, int], Any],
    outdir: Path,
    pigz: bool,
) -> None:
    """Flush one buffered sample/mate payload to its owned output handle."""
    if not buffer:
        return
    key = (sample_name, mate)
    handle = handles.get(key)
    if handle is None:
        handle = _open_writer_handle(
            _writer_output_path(outdir, sample_name, mate, pigz),
            pigz,
        )
        handles[key] = handle
    handle.write(buffer)
    buffer.clear()


def _reader_process(
    reader_id: int,
    assignments: List[Tuple[str, Tuple[Path, Path | None]]],
    config: DemuxRunConfig,
    writer_queues: List[Any],
    writer_route_map: Dict[str, int],
    spool_dir: Path,
    stats_queue: Any,
    progress_queue: Any,
    error_queue: Any,
) -> None:
    """Read and match assigned raw FASTQs, then stream batches to writers."""
    setup_loguru_worker(config.log_level)
    try:
        shard_idx = 0
        reader_raw_total = 0
        reader_matched_total = 0
        for fname, fastq_tuple in assignments:
            matcher = build_matcher(fastq_tuple, config, workers=1)

            def _report_progress(raw_reads: int, matched_reads: int) -> None:
                progress_queue.put(
                    (
                        "reader_progress",
                        reader_id,
                        reader_raw_total + raw_reads,
                        reader_matched_total + matched_reads,
                    )
                )

            matcher.progress_callback = _report_progress
            matcher.progress_interval_reads = max(1, config.chunksize)
            sample_buffers: Dict[str, List[bytearray | None]] = {}
            buffered_reads = 0
            buffered_bytes = 0
            for sample_name, read1_data, read2_data in matcher.iter_output_records():
                buffers = sample_buffers.setdefault(
                    sample_name,
                    [bytearray(), bytearray() if read2_data is not None else None],
                )
                buffers[0].extend(read1_data)
                if read2_data is not None:
                    if buffers[1] is None:
                        buffers[1] = bytearray()
                    buffers[1].extend(read2_data)
                buffered_reads += 1
                buffered_bytes += len(read1_data) + (len(read2_data) if read2_data else 0)
                if buffered_reads >= config.chunksize or buffered_bytes >= config.batch_bytes:
                    shard_idx = _flush_sample_buffers(
                        sample_buffers,
                        writer_queues,
                        writer_route_map,
                        spool_dir,
                        reader_id,
                        shard_idx,
                        config.queue_put_timeout,
                        error_queue,
                    )
                    sample_buffers = {}
                    buffered_reads = 0
                    buffered_bytes = 0
            if sample_buffers:
                shard_idx = _flush_sample_buffers(
                    sample_buffers,
                    writer_queues,
                    writer_route_map,
                    spool_dir,
                    reader_id,
                    shard_idx,
                    config.queue_put_timeout,
                    error_queue,
                )
            matcher._maybe_report_progress(force=True)
            reader_raw_total += matcher.reads_seen
            reader_matched_total += matcher.matched_seen
            stats_queue.put(
                (
                    fname,
                    matcher.barcode_misses,
                    matcher.barcode_hits,
                    matcher.sample_hits,
                    matcher.barcode_boundary_ambiguities,
                )
            )
    except BaseException as exc:
        error_queue.put(
            (
                "reader",
                reader_id,
                type(exc).__name__,
                str(exc),
                traceback.format_exc(),
            )
        )
        raise


def _writer_process(
    writer_id: int,
    writer_queue: Any,
    outdir: Path,
    pigz: bool,
    writer_flush_bytes: int,
    log_level: str,
    error_queue: Any,
) -> None:
    """Own a stable partition of sample outputs and write batches to disk."""
    setup_loguru_worker(log_level)
    handles: Dict[Tuple[str, int], Any] = {}
    sample_buffers: Dict[Tuple[str, int], bytearray] = {}
    try:
        while True:
            payload = writer_queue.get()
            if payload is None:
                break
            shard_bytes = payload.path.read_bytes()
            for entry in payload.entries:
                chunk = shard_bytes[entry.offset:entry.offset + entry.length]
                key = (entry.sample_name, entry.mate)
                sample_buffers.setdefault(key, bytearray()).extend(chunk)
                if len(sample_buffers[key]) >= writer_flush_bytes:
                    _flush_writer_buffer(
                        entry.sample_name,
                        entry.mate,
                        sample_buffers[key],
                        handles,
                        outdir,
                        pigz,
                    )
            try:
                payload.path.unlink()
            except FileNotFoundError:
                pass
    except BaseException as exc:
        error_queue.put(
            (
                "writer",
                writer_id,
                type(exc).__name__,
                str(exc),
                traceback.format_exc(),
            )
        )
        raise
    finally:
        for (sample_name, mate), buffer in sample_buffers.items():
            _flush_writer_buffer(
                sample_name,
                mate,
                buffer,
                handles,
                outdir,
                pigz,
            )
        _close_writer_handles(handles)


def _drain_queue_into_dict(
    qobj: Any,
    file_stats: Dict[str, Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[bytes, int]]],
) -> None:
    """Drain completed file statistics from a multiprocessing queue."""
    while True:
        try:
            fname, barcode_misses, barcode_hits, sample_hits, boundary_ambiguities = _queue_get_nowait(qobj)
        except queue.Empty:
            return
        file_stats[fname] = (barcode_misses, barcode_hits, sample_hits, boundary_ambiguities)


def _drain_progress_queue(
    qobj: Any,
    progress_by_reader: Dict[int, Tuple[int, int]],
) -> bool:
    """Drain reader progress events and return whether any were seen."""
    activity = False
    while True:
        try:
            event = _queue_get_nowait(qobj)
        except queue.Empty:
            return activity
        kind = event[0]
        activity = True
        if kind == "reader_progress":
            _kind, reader_id, raw_reads, matched_reads = event
            progress_by_reader[reader_id] = (raw_reads, matched_reads)
    return activity


def _progress_totals(progress_by_reader: Dict[int, Tuple[int, int]]) -> Tuple[int, int]:
    """Return aggregate raw and matched read counts across readers."""
    raw_total = 0
    matched_total = 0
    for raw_reads, matched_reads in progress_by_reader.values():
        raw_total += raw_reads
        matched_total += matched_reads
    return raw_total, matched_total


def _raise_queued_demux_error(qobj: Any) -> None:
    """Raise the first demux worker error queued back to the parent."""
    try:
        role, worker_id, exc_type, message, tb = _queue_get_nowait(qobj)
    except queue.Empty:
        return
    raise RuntimeError(
        f"demux {role} worker {worker_id} failed with {exc_type}: {message}\n{tb}"
    )


def _raise_on_bad_exit(procs: List[mp.Process], role: str) -> None:
    """Raise if any supervised worker exits non-zero."""
    for proc in procs:
        if proc.exitcode not in (None, 0):
            raise RuntimeError(
                f"demux {role} worker '{proc.name}' exited with code {proc.exitcode}"
            )


def _maybe_log_chunk_progress(
    progress_by_reader: Dict[int, Tuple[int, int]],
    chunk_reads: int,
    next_log_threshold: int,
) -> int:
    """Emit one chunk-based demux progress log when the next threshold is crossed."""
    raw_total, matched_total = _progress_totals(progress_by_reader)
    if raw_total < next_log_threshold:
        return next_log_threshold
    logger.info(
        "demux progress: raw_reads={} matched_reads={}",
        raw_total,
        matched_total,
    )
    return ((raw_total // chunk_reads) + 1) * chunk_reads


def _collect_demux_processes(
    readers: Sequence[mp.Process],
    writers: Sequence[mp.Process],
) -> List[mp.Process]:
    """Return tracked demux workers plus any active spawned children."""
    procs: List[mp.Process] = []
    seen: set[int] = set()
    for proc in list(writers) + list(readers) + list(mp.active_children()):
        key = proc.pid if getattr(proc, "pid", None) is not None else id(proc)
        if key in seen:
            continue
        seen.add(key)
        procs.append(proc)
    return procs


def _process_is_alive(proc: mp.Process) -> bool:
    """Return True when a process object still represents a live child."""
    try:
        return proc.is_alive()
    except Exception:
        return False


def _terminate_processes(procs: List[mp.Process]) -> None:
    """Terminate demux worker processes and force-kill any survivors."""
    for proc in procs:
        try:
            if getattr(proc, "pid", None) is not None and proc.exitcode is None:
                proc.terminate()
        except Exception:
            pass
    deadline = time.monotonic() + 1.0
    while time.monotonic() < deadline:
        alive = [proc for proc in procs if _process_is_alive(proc)]
        if not alive:
            break
        time.sleep(0.05)
    for proc in procs:
        try:
            if _process_is_alive(proc):
                proc.kill()
        except Exception:
            pass
    for proc in procs:
        try:
            proc.join(timeout=0.2)
        except Exception:
            pass


def _close_processes(procs: Sequence[mp.Process]) -> None:
    """Release multiprocessing Process resources after they have exited."""
    for proc in procs:
        try:
            proc.close()
        except Exception:
            pass


def _abort_demux_workers(
    readers: List[mp.Process],
    writers: List[mp.Process],
    queues: Sequence[Any],
) -> List[mp.Process]:
    """Cancel queue feeder joins and terminate all reader/writer workers."""
    procs = _collect_demux_processes(readers, writers)
    for qobj in queues:
        _cancel_queue_join(qobj)
    _terminate_processes(procs)
    return procs


def _close_demux_resources(
    procs: Sequence[mp.Process],
    queues: Sequence[Any],
    *,
    cancel_join: bool,
) -> None:
    """Close queue and process objects without waiting on queue feeder threads."""
    for qobj in queues:
        if cancel_join:
            _cancel_queue_join(qobj)
        _close_queue(qobj)
    _close_processes(procs)


def _cleanup_after_interrupt(
    readers: List[mp.Process],
    writers: List[mp.Process],
    queues: Sequence[Any],
    spool_dir: Path,
) -> None:
    """Tear down demux workers while ignoring repeated Ctrl-C during cleanup."""
    previous_handler = signal.getsignal(signal.SIGINT)
    try:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    except Exception:
        previous_handler = None
    try:
        procs = _abort_demux_workers(readers, writers, queues)
        _close_demux_resources(procs, queues, cancel_join=True)
        shutil.rmtree(spool_dir, ignore_errors=True)
    finally:
        if previous_handler is not None:
            try:
                signal.signal(signal.SIGINT, previous_handler)
            except Exception:
                pass


def _compress_temp_fastq(path: Path) -> None:
    """Compress one temporary FASTQ file with pigz and remove the original."""
    outpath = path.with_suffix(path.suffix + ".gz")
    run_pipeline([["pigz", "-p", "1", "-c", str(path)]], outfile=outpath)
    path.unlink()


def _compress_demux_outputs_with_pigz(
    temp_fastqs: Sequence[Path],
    cores: int,
    log_level: str,
) -> None:
    """Compress the current run's temporary demux FASTQs in parallel with pigz."""
    temp_fastqs = sorted(path for path in temp_fastqs if path.exists())
    if not temp_fastqs:
        return
    logger.info(f"compressing {len(temp_fastqs)} demultiplexed FASTQ files with pigz")
    jobs = {
        path.name: (_compress_temp_fastq, {"path": path})
        for path in temp_fastqs
    }
    run_with_pool(
        jobs,
        log_level=log_level,
        max_workers=min(cores, len(temp_fastqs)),
        msg="Pigz FASTQs",
    )


def run_demux_pipeline(
    filenames_to_fastqs: Dict[str, Tuple[Path, Path | None]],
    config: DemuxRunConfig,
    cores: int,
) -> Dict[str, Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[bytes, int]]]:
    """Demultiplex raw FASTQs with concurrent reader and persistent writer processes."""
    if config.pigz and shutil.which("pigz") is None:
        raise IPyradError(
            "--pigz requires pigz to be installed and on PATH."
        )

    items = list(filenames_to_fastqs.items())
    reader_workers, writer_workers = _reader_writer_counts(
        len(items),
        cores,
        pigz=config.pigz,
    )
    writer_route_map = _build_writer_route_map(
        config.barcodes_to_names,
        config.merge_technical_replicates,
        writer_workers,
    )
    spool_dir = _demux_spool_dir(config.outdir)
    spool_dir.mkdir(parents=True, exist_ok=True)
    logger.info(
        f"demultiplexing on {get_demux_mode_label(config)} with "
        f"{reader_workers} reader(s) and {writer_workers} writer(s) using "
        f"{'pigz' if config.pigz else 'python'} compression"
    )
    for fname, fastq_tuple in items:
        short = tuple(path.name if path else "" for path in fastq_tuple)
        logger.info(f"processing {fname} {short}")

    ctx = mp.get_context("spawn")
    writer_queues = [ctx.Queue(maxsize=2) for _ in range(writer_workers)]
    stats_queue = ctx.SimpleQueue()
    progress_queue = ctx.SimpleQueue()
    error_queue = ctx.SimpleQueue()
    all_queues = writer_queues + [stats_queue, progress_queue, error_queue]

    reader_assignments = _partition_inputs(items, reader_workers)
    readers = [
        ctx.Process(
            name=f"demux-reader-{idx}",
            target=_reader_process,
            args=(
                idx,
                assignments,
                config,
                writer_queues,
                writer_route_map,
                spool_dir,
                stats_queue,
                progress_queue,
                error_queue,
            ),
        )
        for idx, assignments in enumerate(reader_assignments)
    ]
    writers = [
        ctx.Process(
            name=f"demux-writer-{idx}",
            target=_writer_process,
            args=(
                idx,
                writer_queues[idx],
                config.outdir,
                config.pigz,
                config.writer_flush_bytes,
                config.log_level,
                error_queue,
            ),
        )
        for idx in range(writer_workers)
    ]
    all_procs = writers + readers

    file_stats: Dict[str, Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[bytes, int]]] = {}
    progress_by_reader: Dict[int, Tuple[int, int]] = {}
    chunk_reads = max(1, config.chunksize)
    next_log_threshold = chunk_reads
    try:
        for proc in all_procs:
            proc.start()

        while any(proc.is_alive() for proc in readers):
            _raise_queued_demux_error(error_queue)
            _drain_queue_into_dict(stats_queue, file_stats)
            _drain_progress_queue(progress_queue, progress_by_reader)
            next_log_threshold = _maybe_log_chunk_progress(
                progress_by_reader,
                chunk_reads,
                next_log_threshold,
            )
            _raise_on_bad_exit(readers, "reader")
            _raise_on_bad_exit(writers, "writer")
            time.sleep(0.05)

        for proc in readers:
            proc.join()
        _raise_queued_demux_error(error_queue)
        _drain_progress_queue(progress_queue, progress_by_reader)
        next_log_threshold = _maybe_log_chunk_progress(
            progress_by_reader,
            chunk_reads,
            next_log_threshold,
        )
        deadline = time.monotonic() + 0.5
        while len(file_stats) < len(items) and time.monotonic() < deadline:
            _drain_queue_into_dict(stats_queue, file_stats)
            if len(file_stats) < len(items):
                time.sleep(0.02)
        _raise_on_bad_exit(readers, "reader")
        _raise_on_bad_exit(writers, "writer")

        if len(file_stats) != len(items):
            raise RuntimeError(
                "demux reader stage completed without returning stats for all input files"
            )

        for writer_queue in writer_queues:
            _put_with_supervision(
                writer_queue,
                None,
                error_queue,
                writers,
                timeout=config.queue_put_timeout,
            )

        while any(proc.is_alive() for proc in writers):
            _raise_queued_demux_error(error_queue)
            _drain_progress_queue(progress_queue, progress_by_reader)
            next_log_threshold = _maybe_log_chunk_progress(
                progress_by_reader,
                chunk_reads,
                next_log_threshold,
            )
            _raise_on_bad_exit(writers, "writer")
            time.sleep(0.05)

        for proc in writers:
            proc.join()
        _raise_queued_demux_error(error_queue)
        _raise_on_bad_exit(writers, "writer")

    except KeyboardInterrupt:
        _cleanup_after_interrupt(readers, writers, all_queues, spool_dir)
        raise
    except Exception:
        _abort_demux_workers(readers, writers, all_queues)
        raise
    finally:
        _close_demux_resources(all_procs, all_queues, cancel_join=False)
        shutil.rmtree(spool_dir, ignore_errors=True)

    if config.pigz:
        paired_output = any(fastq_tuple[1] is not None for _, fastq_tuple in items)
        temp_fastqs = [
            _writer_output_path(config.outdir, sample_name, 1, True)
            for sample_name in sorted(writer_route_map)
        ]
        if paired_output:
            temp_fastqs.extend(
                _writer_output_path(config.outdir, sample_name, 2, True)
                for sample_name in sorted(writer_route_map)
            )
        _compress_demux_outputs_with_pigz(temp_fastqs, cores, config.log_level)
    return file_stats
