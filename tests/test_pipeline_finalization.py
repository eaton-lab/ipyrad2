import os, subprocess, sys, time
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
    def failed_reap(procs):
        for p in procs: p.wait(timeout=2)
        raise subprocess.TimeoutExpired('simulated uninterruptible exit',.2)
    monkeypatch.setattr(m,'_wait_for_pipeline',failed_reap)
    with pytest.raises(subprocess.TimeoutExpired): list(stream)
    assert not any(p in m._CHILD_PROCS for p in stream.procs)
    assert all(not os.path.exists(f.name) for f in stream.stderr_tmpfiles)
