"""Bounded ordered alignment, failure cleanup and actual native-call accounting."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import threading
import time

import pytest

from ipyrad2.denovo import align
from ipyrad2.denovo.alignment_metrics import AlignmentMetrics


def jobs(n, length=10):
    for i in range(n):
        yield i, dict(locus_id=i, locus_name=f'locus_{i}',
                      record=[('a', 'A' * length), ('b', 'C' * length)],
                      mafft_binary='test', output_spacer_len=0)


def test_slow_first_locus_bounds_all_unwritten_results(monkeypatch):
    release = threading.Event()
    later_done = threading.Event()
    generated = []
    metrics = AlignmentMetrics()

    def source():
        for idx, kw in jobs(40):
            generated.append(idx)
            yield idx, kw

    def worker(**kw):
        idx = kw['locus_id']
        if idx == 0:
            assert release.wait(5)
        if idx == 3:
            later_done.set()
        return idx, kw['locus_name'], 'ACGT', False

    monkeypatch.setattr(align, 'worker_build_consensus', worker)
    with ThreadPoolExecutor(1) as executor:
        task = executor.submit(lambda: list(align._iter_threaded_alignment_results(
            source(), max_workers=2, heartbeat_s=.05, metrics=metrics, max_outstanding=4)))
        try:
            assert later_done.wait(3)
            assert generated == list(range(4))
            assert not task.done()
            deadline = time.monotonic() + 3
            while not metrics.details['backpressure_events'] and time.monotonic() < deadline:
                time.sleep(.001)
            assert metrics.details['backpressure_events'] > 0
        finally:
            release.set()
        result = task.result(timeout=5)
    assert [i for i, _ in result] == list(range(40))
    assert metrics.details['peak_outstanding_loci'] == 4
    assert metrics.details['peak_buffered_results'] <= 4
    assert metrics.details['backpressure_events'] > 0


def test_byte_reservations_bound_work_and_isolate_oversized_locus(monkeypatch):
    lock = threading.Lock()
    active = set()
    observed = []
    metrics = AlignmentMetrics()
    lengths = [20, 500, 20, 20, 20]

    def source():
        for idx, length in enumerate(lengths):
            yield idx, dict(locus_id=idx, locus_name=f'locus_{idx}',
                           record=[('a', 'A' * length), ('b', 'C' * length)],
                           mafft_binary='test', output_spacer_len=0)

    def worker(**kw):
        idx = kw['locus_id']
        with lock:
            active.add(idx)
            observed.append(frozenset(active))
        time.sleep(.005)
        with lock:
            active.remove(idx)
        return idx, kw['locus_name'], 'A' * lengths[idx], False

    monkeypatch.setattr(align, 'worker_build_consensus', worker)
    result = list(align._iter_threaded_alignment_results(source(), max_workers=4,
        heartbeat_s=.1, metrics=metrics, output_budget_bytes=100))
    assert [i for i, _ in result] == list(range(5))
    assert all(len(s) == 1 for s in observed if 1 in s)
    assert metrics.details['oversized_loci'] == 1
    assert metrics.details['peak_reserved_payload_bytes'] == 1007
    assert metrics.details['peak_outstanding_loci'] <= 2
    assert metrics.details['peak_buffered_payload_bytes'] <= 1007


@pytest.mark.parametrize('where', ['worker', 'iterator', 'close'])
def test_abort_kills_children_and_stops_dispatch(monkeypatch, where):
    release = threading.Event()
    killed = []
    seen = []

    def kill():
        killed.append(True)
        release.set()

    def worker(**kw):
        if where == 'worker':
            raise RuntimeError('native failure')
        if kw['locus_id']:
            release.wait(3)
        return kw['locus_id'], kw['locus_name'], 'A', False

    def source():
        for idx, kw in jobs(30):
            if where == 'iterator' and idx == 1:
                raise ValueError('input failure')
            seen.append(idx)
            yield idx, kw

    monkeypatch.setattr(align.pipeline_module, '_kill_all_children', kill)
    monkeypatch.setattr(align, 'worker_build_consensus', worker)
    it = align._iter_threaded_alignment_results(source(), max_workers=2, heartbeat_s=.1)
    if where == 'close':
        assert next(it)[0] == 0
        it.close()
    else:
        with pytest.raises((ValueError, RuntimeError), match='failure'):
            list(it)
    assert killed == [True]
    assert len(seen) <= 2


def test_native_call_counts_include_multiple_calls_and_parse_failures(monkeypatch):
    @contextmanager
    def context(binary, text):
        yield binary, {}, None

    bad = False
    def pipeline(argv, **kw):
        output = '>unexpected\nACGT\n' if bad else kw['stdin_text']
        return 0, output.encode(), b''

    monkeypatch.setattr(align, 'mafft_process_context', context)
    monkeypatch.setattr(align, 'run_pipeline', pipeline)
    metrics = AlignmentMetrics()
    with metrics.worker():
        for _ in range(2):
            assert align.mafft_align_one([('a','ACGT'),('b','ACGA')], 'test') == [('a','ACGT'),('b','ACGA')]
            assert align.consensus_from_aligned([('a','ACGT'),('b','ACGA')]) == 'ACGA'
        bad = True
        with pytest.raises(RuntimeError, match='unexpected'):
            align.mafft_align_one([('a','ACGT')], 'test')
    d = metrics.details
    assert (d['mafft_invocations'], d['mafft_succeeded'], d['mafft_failed']) == (3,2,1)
    assert sum(d['mafft_duration_counts']) == 3
    assert d['loci_completed'] == 1
    assert d['mafft_execute_seconds_sum'] > 0 and d['mafft_parse_seconds_sum'] > 0
    assert d['consensus_seconds_sum'] > 0
    assert d['mafft_duration_min_seconds'] <= d['mafft_duration_max_seconds']


def test_missing_and_trivial_alignment_have_explicit_zero_calls(tmp_path):
    mapping = tmp_path/'mapping.tsv';mapping.write_text('locus\tcore\n1\ta\n')
    summary = tmp_path/'summary.tsv'
    summary.write_text('seed\tsample\trecord_type\tcluster_sequence\tarm_boundary\na\ts\tsingle\tACGT\t4\n')
    for mode in ['none','mafft']:
        d = {}
        out = tmp_path/f'{mode}.fa'
        result = align.write_ordered_consensus_stream_to_file(mapping,summary,out,'unused',alignment_mode=mode,performance=d)
        assert out.read_text() == '>locus_1\nACGT\n'
        assert result.performance is d
        assert d['mafft_invocations'] == d['loci_submitted'] == 0
        assert d['mafft_duration_min_seconds'] is None
        assert d['planning_seconds'] >= 0 and d['fasta_write_seconds'] > 0


def test_result_cannot_silently_exceed_reserved_bound(monkeypatch):
    monkeypatch.setattr(align, 'worker_build_consensus', lambda **kw: (0,'locus_0','A'*1000,False))
    monkeypatch.setattr(align.pipeline_module, '_kill_all_children', lambda:None)
    with pytest.raises(RuntimeError, match='sequence bound for locus_0'):
        list(align._iter_threaded_alignment_results(jobs(1),max_workers=1,heartbeat_s=.1))
