"""Bounded, run-local alignment measurements; no changes to scientific inputs."""
from bisect import bisect_right
from contextlib import contextmanager
from functools import wraps
import threading
import time


_LOCAL = threading.local()
_BOUNDS = (.01, .025, .05, .1, .25, .5, 1., 2., 5., 10., 30., 60., 300., 900.)
_SUMS = ('worker_compute_seconds_sum', 'mafft_setup_cleanup_seconds_sum',
         'mafft_execute_seconds_sum', 'mafft_parse_seconds_sum', 'consensus_seconds_sum')


def add_seconds(name, seconds):
    metrics = getattr(_LOCAL, 'metrics', None)
    if metrics is not None:
        metrics[name] += seconds


def timed_phase(name):
    def decorate(func):
        @wraps(func)
        def measured(*args, **kwargs):
            if getattr(_LOCAL, 'metrics', None) is None:
                return func(*args, **kwargs)
            started = time.perf_counter()
            try:
                return func(*args, **kwargs)
            finally:
                add_seconds(name, time.perf_counter() - started)
        return measured
    return decorate


def measured_mafft(func):
    @wraps(func)
    def measured(*args, **kwargs):
        metrics = getattr(_LOCAL, 'metrics', None)
        if metrics is None:
            return func(*args, **kwargs)
        started = time.perf_counter()
        metrics['mafft_invocations'] += 1
        try:
            result = func(*args, **kwargs)
        except BaseException:
            metrics['mafft_failed'] += 1
            raise
        else:
            metrics['mafft_succeeded'] += 1
            return result
        finally:
            elapsed = time.perf_counter() - started
            metrics['mafft_duration_counts'][bisect_right(_BOUNDS, elapsed)] += 1
            metrics['mafft_duration_seconds_sum'] += elapsed
            old = metrics['mafft_duration_min_seconds']
            metrics['mafft_duration_min_seconds'] = elapsed if old is None else min(old, elapsed)
            metrics['mafft_duration_max_seconds'] = max(metrics['mafft_duration_max_seconds'] or 0, elapsed)
    return measured


def _worker_metrics():
    return dict.fromkeys(_SUMS, 0.0) | dict(
        mafft_invocations=0, mafft_succeeded=0, mafft_failed=0,
        mafft_duration_counts=[0] * (len(_BOUNDS) + 1),
        mafft_duration_seconds_sum=0.0,
        mafft_duration_min_seconds=None, mafft_duration_max_seconds=None,
    )


class AlignmentMetrics:
    """Merge small per-worker accumulators; never retain per-call records."""
    def __init__(self, details=None):
        self.details = {} if details is None else details
        self.details.update(_worker_metrics(), schema_version='1.0',
                            read_summary_seconds=0.0, planning_seconds=0.0,
                            job_build_seconds=0.0, dispatch_collect_seconds=0.0,
                            fasta_write_seconds=0.0, backpressure_wait_seconds=0.0,
                            backpressure_events=0, loci_submitted=0, loci_completed=0,
                            loci_failed=0, peak_buffered_results=0,
                            peak_buffered_payload_bytes=0, peak_outstanding_loci=0,
                            peak_reserved_payload_bytes=0, oversized_loci=0,
                            outstanding_locus_limit=0, output_payload_budget_bytes=0,
                            mafft_duration_upper_bounds_seconds=list(_BOUNDS) + [None],
                            timing_definition='Worker/native times are sums across parallel loci and overlap. '
                            'MAFFT execute includes process launch and waiting; setup/cleanup includes scratch/runtime preparation. '
                            'Histogram measures complete MAFFT calls including parsing and failures; bins are lower-inclusive, upper-exclusive.',
                            payload_definition='UTF-8 bytes of result locus names and consensus strings, excluding Python object overhead, '
                            'input records, plans and native-tool memory. Admission reserves a conservative output bound. '
                            'An oversized locus runs alone.')
        self.lock = threading.Lock()

    @contextmanager
    def worker(self):
        previous = getattr(_LOCAL, 'metrics', None)
        metrics = _worker_metrics()
        _LOCAL.metrics = metrics
        started = time.perf_counter()
        failed = False
        try:
            yield
        except BaseException:
            failed = True
            raise
        finally:
            metrics['worker_compute_seconds_sum'] = time.perf_counter() - started
            _LOCAL.metrics = previous
            with self.lock:
                data = self.details
                data['loci_failed' if failed else 'loci_completed'] += 1
                for key, value in metrics.items():
                    if key == 'mafft_duration_counts':
                        data[key] = [a + b for a, b in zip(data[key], value)]
                    elif key.endswith('_min_seconds') or key.endswith('_max_seconds'):
                        if value is not None:
                            fn = min if key.endswith('_min_seconds') else max
                            data[key] = value if data[key] is None else fn(data[key], value)
                    else:
                        data[key] += value
