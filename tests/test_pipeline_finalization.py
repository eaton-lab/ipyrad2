import os, subprocess, sys, time
from ipyrad2.utils.parallel import run_pipeline
import pytest
from ipyrad2.utils.parallel import stream_pipeline_lines, PipelineTimeoutError
from ipyrad2.utils.parallel import pipeline as m

def emitter(delay=1.4, code=0):
    return [sys.executable,'-c',f"import os,sys,time; os.write(1,b'complete\\n'); os.close(1); time.sleep({delay}); sys.stderr.write('late status\\n'); sys.exit({code})"]

def test_delayed_success_is_not_killed():
    assert list(stream_pipeline_lines([emitter()])) == ['complete']

def test_delayed_failure_is_not_accepted():
    with pytest.raises(RuntimeError,match='rc=17'):
        list(stream_pipeline_lines([emitter(code=17)]))

def test_multistage_upstream_delayed_failure_is_not_accepted():
    copy=[sys.executable,'-c','import sys; sys.stdout.write(sys.stdin.read())']
    with pytest.raises(RuntimeError,match='rc=17'):
        list(stream_pipeline_lines([emitter(code=17),copy]))

def test_finalization_is_bounded_and_reaps_children(monkeypatch):
    monkeypatch.setattr(m,'STREAM_FINALIZE_TIMEOUT_S',.15)
    stream=stream_pipeline_lines([emitter(delay=10)])
    start=time.monotonic()
    with pytest.raises(PipelineTimeoutError): list(stream)
    assert time.monotonic()-start<2
    assert all(p.poll() is not None for p in stream.procs)
    assert not any(p in m._CHILD_PROCS for p in stream.procs)
    assert all(not os.path.exists(f.name) for f in stream.stderr_tmpfiles)

def test_cleanup_still_runs_if_reaping_raises(monkeypatch):
    stream=stream_pipeline_lines([emitter(delay=0)])
    def failed_reap(procs, **kwargs):
        for p in procs: p.wait(timeout=2)
        raise subprocess.TimeoutExpired('simulated uninterruptible exit',.2)
    monkeypatch.setattr(m,'_wait_for_pipeline',failed_reap)
    with pytest.raises(subprocess.TimeoutExpired): list(stream)
    assert not any(p in m._CHILD_PROCS for p in stream.procs)
    assert all(not os.path.exists(f.name) for f in stream.stderr_tmpfiles)

def test_early_close_releases_pipe_and_uses_separate_reap_budget(monkeypatch):
    # Simulate a layout probe ending after its first complete alignment line.
    producer = [sys.executable, '-c',
                "import os,time; os.write(1,b'first\\n'); time.sleep(60)"]
    stream = stream_pipeline_lines([producer])
    original = m._wait_for_pipeline
    budgets = []
    def reap(procs, *, timeout_s):
        budgets.append(timeout_s)
        assert stream.last.stdout.closed
        return original(procs, timeout_s=timeout_s)
    monkeypatch.setattr(m, '_wait_for_pipeline', reap)
    with stream as lines:
        assert next(lines) == 'first'
    assert budgets == [m.STREAM_EARLY_CLOSE_REAP_TIMEOUT_S]
    assert all(p.poll() is not None for p in stream.procs)
    assert not any(p in m._CHILD_PROCS for p in stream.procs)


def test_delayed_reap_requires_longer_grace(monkeypatch):
    """Model the kernel delaying SIGKILL completion on a shared filesystem."""
    class DelayedChild:
        args = ['samtools', 'view', 'bam']
        def __init__(self):
            self.exit_at = time.monotonic() + .08
        def poll(self):
            return -9 if time.monotonic() >= self.exit_at else None
    child = DelayedChild()
    with pytest.raises(PipelineTimeoutError, match='cleanup exceeded'):
        m._wait_for_pipeline([child], timeout_s=.02)
    assert m._wait_for_pipeline([child], timeout_s=.5) == [-9]


def test_cleanup_budget_is_shared_not_per_process():
    class StuckChild:
        args = ['samtools', 'view', 'bam']
        def poll(self):
            return None
    start = time.monotonic()
    with pytest.raises(PipelineTimeoutError, match='unreaped children'):
        m._wait_for_pipeline([StuckChild() for _ in range(12)], timeout_s=.1)
    assert time.monotonic() - start < .5


def test_explicit_timeout_keeps_short_cleanup_and_original_error(monkeypatch):
    original = m._wait_for_pipeline
    budgets = []
    def delayed_reap(procs, *, timeout_s):
        budgets.append(timeout_s)
        original(procs, timeout_s=timeout_s)
        raise PipelineTimeoutError('simulated slow cleanup')
    monkeypatch.setattr(m, '_wait_for_pipeline', delayed_reap)
    start = time.monotonic()
    with pytest.raises(PipelineTimeoutError, match=r'timed out after.*timeout=0.1s') as exc:
        run_pipeline([[sys.executable, '-c', 'import time; time.sleep(30)']], timeout_s=.1)
    assert time.monotonic() - start < 2
    assert isinstance(exc.value.__cause__, PipelineTimeoutError)
    assert budgets == [m.PIPELINE_FAILURE_REAP_TIMEOUT_S]
    assert m.PIPELINE_FAILURE_REAP_TIMEOUT_S < m.STREAM_EARLY_CLOSE_REAP_TIMEOUT_S


def test_primary_failure_survives_delayed_cleanup(monkeypatch):
    original = m._wait_for_pipeline
    def delayed_reap(procs, *, timeout_s):
        original(procs, timeout_s=timeout_s)
        raise PipelineTimeoutError('simulated slow cleanup')
    monkeypatch.setattr(m, '_wait_for_pipeline', delayed_reap)
    start = time.monotonic()
    with pytest.raises(RuntimeError, match='rc=17') as exc:
        run_pipeline([[sys.executable, '-c', 'import sys; print("bad", file=sys.stderr); sys.exit(17)']])
    assert time.monotonic() - start < 2
    assert 'bad' in str(exc.value)
    assert isinstance(exc.value.__cause__, PipelineTimeoutError)


def test_partial_launch_cleans_up_previous_child(tmp_path, monkeypatch):
    original = m.safe_popen
    launched = []
    def track(argv, **kwargs):
        proc = original(argv, **kwargs)
        launched.append(proc)
        return proc
    monkeypatch.setattr(m, 'safe_popen', track)
    with pytest.raises(FileNotFoundError):
        run_pipeline([[sys.executable, '-c', 'import time; time.sleep(30)'],
                      [str(tmp_path / 'nonexistent-executable')]])
    assert len(launched) == 1
    assert launched[0].poll() is not None
    assert launched[0] not in m._CHILD_PROCS


@pytest.mark.parametrize("timeout", [float("inf"), float("nan"), -1, 0, True])
def test_invalid_execution_deadline_rejected(timeout):
    with pytest.raises(ValueError, match="positive number"):
        run_pipeline([[sys.executable, '-c', 'pass']], timeout_s=timeout)



def test_early_close_does_not_hide_already_observed_failure():
    stream = stream_pipeline_lines([emitter(delay=0, code=17)])
    assert next(stream) == "complete"
    # Ensure the command has really failed before the caller cancels.
    assert stream.procs[0].wait(timeout=2) == 17
    with pytest.raises(RuntimeError, match="rc=17"):
        stream.close()
    assert stream.last.stdout.closed
    assert not any(p in m._CHILD_PROCS for p in stream.procs)


def test_consumer_exception_is_not_masked_by_cleanup_failure(monkeypatch):
    original = m._wait_for_pipeline
    budgets = []
    def fail_reap(procs, *, timeout_s):
        original(procs, timeout_s=timeout_s)
        budgets.append(timeout_s)
        raise PipelineTimeoutError("simulated slow cleanup")
    monkeypatch.setattr(m, "_wait_for_pipeline", fail_reap)
    with pytest.raises(ValueError, match="consumer error"):
        with stream_pipeline_lines([emitter(delay=10)]) as lines:
            assert next(lines) == "complete"
            raise ValueError("consumer error")
    assert budgets == [m.PIPELINE_FAILURE_REAP_TIMEOUT_S]



def test_multistage_timeout_remains_fast_and_atomic(tmp_path):
    output = tmp_path / "result"
    output.write_bytes(b"previous")
    stubborn = [sys.executable, "-c",
                "import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)"]
    start = time.monotonic()
    with pytest.raises(PipelineTimeoutError, match="timeout=0.2s"):
        run_pipeline([stubborn] * 4, outfile=output, timeout_s=.2)
    assert time.monotonic() - start < 2
    assert output.read_bytes() == b"previous"
    assert not list(tmp_path.glob(".result.*.tmp"))


def test_unreaped_child_remains_registered_until_later_poll():
    class Child:
        stdin = stdout = stderr = None
        returncode = None
        def poll(self):
            return self.returncode
    child = Child()
    with m._CHILD_PROCS_LOCK:
        m._CHILD_PROCS.add(child)
    try:
        m._cleanup_pipeline_processes([child], [])
        assert child in m._CHILD_PROCS
        child.returncode = -9
        proc = m.safe_popen([sys.executable, "-c", "pass"])
        try:
            assert child not in m._CHILD_PROCS
            proc.wait(timeout=2)
        finally:
            m._cleanup_pipeline_processes([proc], [])
    finally:
        with m._CHILD_PROCS_LOCK:
            m._CHILD_PROCS.discard(child)
