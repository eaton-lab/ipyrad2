"""Ordered pool delivery bounds completed results as well as running work."""
from concurrent.futures import Future
import time
import sys
from pathlib import Path

import pytest

from ipyrad2.utils.parallel import pool, run_pipeline


def identity(value):
    return value


class ControlledPool(pool._ManagedProcessPool):
    def __init__(self, failure=None):
        self.submitted = []
        self.aborted = False
        self.failure = failure

    def submit(self, func, kwargs):
        if self.failure == len(self.submitted):
            raise RuntimeError('submission failure')
        future = Future()
        self.submitted.append((future, kwargs['value']))
        return future

    def collect(self):
        pass

    def abort(self, *args, **kwargs):
        self.aborted = True


def jobs(count=30):
    # Duplicated, non-sortable keys must not define delivery order.
    for i in range(count):
        yield None if i % 2 else 'repeat', (identity, {'value': i})


def test_stalled_first_batch_bounds_running_and_completed(monkeypatch):
    controlled = ControlledPool()
    consumed = []
    waits = []
    def wait(futures, **kwargs):
        outstanding = len(controlled.submitted) - len(consumed)
        assert outstanding <= 4
        waits.append(outstanding)
        # Complete all later jobs before the first, repeatedly keeping it stalled.
        for future, value in controlled.submitted:
            if future in futures and (value != 0 or len(waits) >= 4):
                future.set_result(value)
        return {f for f in futures if f.done()}, set()
    monkeypatch.setattr(pool, 'wait', wait)
    for key, result in controlled.iter_results(jobs(), max_inflight=4, ordered=True):
        assert result == len(consumed)
        consumed.append(result)
    assert consumed == list(range(30))
    assert waits[:4] == [4, 4, 4, 4]
    assert not controlled.aborted


def test_later_worker_error_is_detected_while_first_stalls(monkeypatch):
    controlled = ControlledPool()
    def wait(futures, **kwargs):
        failed = controlled.submitted[2][0]
        failed.set_exception(ValueError('later batch failed'))
        return {failed}, set(futures) - {failed}
    monkeypatch.setattr(pool, 'wait', wait)
    with pytest.raises(pool.ParallelJobError) as error:
        next(controlled.iter_results(jobs(), max_inflight=4, ordered=True))
    assert isinstance(error.value.original, ValueError)
    assert controlled.aborted and len(controlled.submitted) == 4


def test_iterator_error_aborts_submitted_work():
    controlled = ControlledPool()
    def broken():
        yield from jobs(1)
        raise OSError('input failure')
    with pytest.raises(OSError, match='input failure'):
        next(controlled.iter_results(broken(), max_inflight=4, ordered=True))
    assert controlled.aborted and len(controlled.submitted) == 1


def test_submission_failure_aborts_submitted_work():
    controlled = ControlledPool(failure=1)
    with pytest.raises(pool.ParallelJobError) as error:
        next(controlled.iter_results(jobs(), max_inflight=4, ordered=True))
    assert isinstance(error.value.original, RuntimeError)
    assert controlled.aborted


def test_close_closes_primed_source_and_aborts(monkeypatch):
    controlled = ControlledPool()
    closed = []
    def source():
        try:
            yield from jobs()
        finally:
            closed.append(True)
    def wait(futures, **kwargs):
        f, value = controlled.submitted[0]
        f.set_result(value)
        return {f}, set(futures) - {f}
    monkeypatch.setattr(pool, 'wait', wait)
    results = controlled.iter_results(pool._prime_jobs_iter(source()), max_inflight=1, ordered=True)
    assert next(results) == ('repeat', 0)
    results.close()
    assert closed == [True] and controlled.aborted


@pytest.mark.parametrize('value', [None, 1, 'true'])
def test_ordered_requires_bool(value):
    with pytest.raises(ValueError, match='ordered'):
        list(pool.run_with_pool_iter(jobs(1), 'ERROR', ordered=value))


def native_child(marker):
    # Delay creates an observable side effect if early close misses a native child.
    command = [sys.executable, '-c',
        'import pathlib,time; '
        f'pathlib.Path({(str(marker)+".started")!r}).write_text("started"); time.sleep(1.5); '
        f'pathlib.Path({marker!r}).write_text("escaped")']
    run_pipeline([command])
    return 1


def wait_for_child(value, started):
    deadline = time.monotonic() + 5
    while not Path(started).exists():
        if time.monotonic() >= deadline:
            raise TimeoutError('native child did not start')
        time.sleep(.01)
    return value


def test_ordered_writer_failure_close_stops_native_children(tmp_path):
    marker = tmp_path/'escaped'
    work = [('first', (wait_for_child, {'value': 0, 'started': str(marker)+'.started'})),
            ('child', (native_child, {'marker': str(marker)}))]
    results = pool.run_with_pool_iter(work, 'ERROR', max_workers=2, ordered=True)
    with pytest.raises(OSError, match='writer'):
        try:
            assert next(results) == ('first', 0)
            assert Path(str(marker)+'.started').exists()
            raise OSError('writer failed')
        finally:
            results.close()
    time.sleep(1.6)
    assert not marker.exists()
