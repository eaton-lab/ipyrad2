# tests/test_ctrl_c_parallel.py

from __future__ import annotations

import os
import signal
import time
import unittest
import multiprocessing as mp

from ipyrad2.utils.parallel3 import run_with_pool, run_with_pool_iter, run_pipeline


# ---------------------------- helpers for tests ---------------------------- #

def _sleep_worker(seconds: float) -> str:
    time.sleep(seconds)
    return f"slept:{seconds:.1f}s"


def _pipeline_sleep_worker(seconds: float) -> int:
    rc, _out, _err = run_pipeline([["bash", "-lc", f"sleep {seconds}"]])
    return rc


def _run_pool_until_interrupted(
    mode: str,
    n_jobs: int,
    sleep_each: float,
    log_level: str = "INFO",
) -> None:
    if mode == "with_pool":
        jobs = {f"job-{i}": (_sleep_worker, {"seconds": sleep_each}) for i in range(n_jobs)}
        run_with_pool(jobs, log_level, max_workers=2)

    elif mode == "with_pool_iter":
        def job_iter():
            for i in range(n_jobs):
                yield (f"job-{i}", (_sleep_worker, {"seconds": sleep_each}))

        for _key, _res in run_with_pool_iter(
            job_iter(), log_level, max_workers=2, max_inflight=2
        ):
            pass

    elif mode == "pipeline":
        jobs = {f"pipe-{i}": (_pipeline_sleep_worker, {"seconds": sleep_each}) for i in range(n_jobs)}
        run_with_pool(jobs, log_level, max_workers=2)

    else:
        raise ValueError(f"unknown mode: {mode}")


def _spawn_and_sigint(
    target,
    *args,
    grace: float = 2.0,
    start_delay: float = 0.60,
    epsilon: float = 0.15,
) -> int:
    """
    Spawn child, let it start, send SIGINT, and assert timely exit.

    grace
        Required exit time window (seconds).
    start_delay
        Time to let the child initialize the pool before signaling.
    epsilon
        Small tolerance for scheduler jitter beyond 'grace'.
    """
    ctx = mp.get_context("spawn")
    p = ctx.Process(target=target, args=args)
    p.start()

    # Allow worker pool and subprocess PGIDs to initialize.
    time.sleep(start_delay)

    # Deliver SIGINT; a second nudge reduces flakiness on some kernels.
    try:
        os.kill(p.pid, signal.SIGINT)
        time.sleep(0.05)
        os.kill(p.pid, signal.SIGINT)
    except ProcessLookupError:
        pass  # already gone

    t0 = time.time()
    p.join(timeout=grace + epsilon)
    elapsed = time.time() - t0

    if p.is_alive():
        # Hard kill to avoid hanging the test suite.
        try:
            os.kill(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        p.join(0.5)
        raise AssertionError(
            f"child did not exit within {grace:.2f}s (+{epsilon:.2f}s epsilon) after SIGINT; "
            f"elapsed={elapsed:.3f}s"
        )
    return p.exitcode


# --------------------------------- tests ---------------------------------- #

@unittest.skipUnless(os.name == "posix", "Ctrl-C/PGID tests require POSIX")
class TestCtrlCShutdown(unittest.TestCase):
    def test_fast_shutdown_with_pool(self):
        code = _spawn_and_sigint(
            _run_pool_until_interrupted, "with_pool", 50, 5.0, "WARNING", grace=2.0
        )
        self.assertGreaterEqual(code, 128, f"unexpected exit code: {code}")

    def test_fast_shutdown_with_pool_iter(self):
        code = _spawn_and_sigint(
            _run_pool_until_interrupted, "with_pool_iter", 200, 3.0, "WARNING", grace=2.0
        )
        self.assertGreaterEqual(code, 128, f"unexpected exit code: {code}")

    def test_fast_shutdown_pipeline_children(self):
        code = _spawn_and_sigint(
            _run_pool_until_interrupted, "pipeline", 20, 10.0, "WARNING", grace=2.0
        )
        self.assertGreaterEqual(code, 128, f"unexpected exit code: {code}")

    def test_normal_completion_returns_results(self):
        jobs = {f"ok-{i}": (_sleep_worker, {"seconds": 0.05}) for i in range(8)}
        res = run_with_pool(jobs, log_level="WARNING", max_workers=2)
        self.assertEqual(set(res.keys()), set(jobs.keys()))
        for v in res.values():
            self.assertTrue(v.startswith("slept:0."), f"unexpected worker result: {v}")


if __name__ == "__main__":
    unittest.main()
