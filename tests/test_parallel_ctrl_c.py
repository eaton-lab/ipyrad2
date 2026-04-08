from __future__ import annotations

import multiprocessing as mp
import os
from pathlib import Path
import signal
import time

import pytest

from ipyrad2.denovo import align as align_module
from ipyrad2.utils.parallel import run_pipeline
from ipyrad2.utils.parallel import run_with_pool
from ipyrad2.utils.parallel import run_with_pool_iter


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
        jobs = {f"job-{idx}": (_sleep_worker, {"seconds": sleep_each}) for idx in range(n_jobs)}
        run_with_pool(jobs, log_level, max_workers=2)
        return

    if mode == "with_pool_iter":
        def _job_iter():
            for idx in range(n_jobs):
                yield f"job-{idx}", (_sleep_worker, {"seconds": sleep_each})

        for _key, _res in run_with_pool_iter(
            _job_iter(),
            log_level,
            max_workers=2,
            max_inflight=2,
        ):
            pass
        return

    if mode == "pipeline":
        jobs = {
            f"pipe-{idx}": (_pipeline_sleep_worker, {"seconds": sleep_each})
            for idx in range(n_jobs)
        }
        run_with_pool(jobs, log_level, max_workers=2)
        return

    raise ValueError(f"unknown mode: {mode}")


def _run_denovo_align_until_interrupted(
    mapping_tsv: Path,
    summary_tsv: Path,
    out_fa: Path,
    mafft_binary: Path,
    pid_file: Path,
) -> None:
    os.environ["IPYRAD2_TEST_MAFFT_PID_FILE"] = str(pid_file)
    align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary=str(mafft_binary),
        cores=2,
        mafft_timeout_s=30.0,
    )


def _spawn_and_sigint(
    target,
    *args,
    grace: float = 2.0,
    start_delay: float = 0.60,
    epsilon: float = 0.15,
) -> int:
    ctx = mp.get_context("spawn")
    proc = ctx.Process(target=target, args=args)
    proc.start()

    time.sleep(start_delay)

    try:
        os.kill(proc.pid, signal.SIGINT)
        time.sleep(0.05)
        os.kill(proc.pid, signal.SIGINT)
    except ProcessLookupError:
        pass

    start = time.time()
    proc.join(timeout=grace + epsilon)
    elapsed = time.time() - start

    if proc.is_alive():
        try:
            os.kill(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.join(0.5)
        raise AssertionError(
            f"child did not exit within {grace:.2f}s (+{epsilon:.2f}s epsilon); "
            f"elapsed={elapsed:.3f}s"
        )
    return proc.exitcode


pytestmark = pytest.mark.skipif(os.name != "posix", reason="Ctrl-C tests require POSIX")


def _assert_interrupted_exit(code: int) -> None:
    assert code == -signal.SIGINT or code >= 128


def _pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def test_fast_shutdown_with_pool() -> None:
    code = _spawn_and_sigint(
        _run_pool_until_interrupted,
        "with_pool",
        50,
        5.0,
        "WARNING",
        grace=2.0,
    )
    _assert_interrupted_exit(code)


def test_fast_shutdown_with_pool_iter() -> None:
    code = _spawn_and_sigint(
        _run_pool_until_interrupted,
        "with_pool_iter",
        200,
        3.0,
        "WARNING",
        grace=2.0,
    )
    _assert_interrupted_exit(code)


def test_fast_shutdown_pipeline_children() -> None:
    code = _spawn_and_sigint(
        _run_pool_until_interrupted,
        "pipeline",
        20,
        10.0,
        "WARNING",
        grace=2.0,
    )
    _assert_interrupted_exit(code)


def test_fast_shutdown_denovo_alignment_children(tmp_path: Path) -> None:
    mapping_tsv = tmp_path / "loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"
    pid_file = tmp_path / "mafft_pids.txt"
    mafft = tmp_path / "mafft"

    mapping_tsv.write_text(
        "locus\tsample\tn_reads\tn_unique\tlength\tmerged\tcluster_id\tcore\n"
        "1\ts1\t5\t1\t10\t0\t0\ts1;J1\n"
        "1\ts2\t6\t1\t10\t0\t0\ts2;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        "sample\tcluster_id\tseed\tlength\tcluster_length\tn_unique\tn_reads\trecord_type\tcluster_sequence\tarm_boundary\n"
        "s1\t0\ts1;J1\t10\t10\t1\t5\tsingle\tAAAAAAAAAA\t10\n"
        "s2\t0\ts2;J2\t10\t10\t1\t6\tsingle\tAAAAAAAATA\t10\n",
        encoding="utf-8",
    )
    mafft.write_text(
        "#!/bin/sh\n"
        "echo $$ >> \"$IPYRAD2_TEST_MAFFT_PID_FILE\"\n"
        "sleep 10\n",
        encoding="utf-8",
    )
    mafft.chmod(0o755)

    ctx = mp.get_context("spawn")
    proc = ctx.Process(
        target=_run_denovo_align_until_interrupted,
        args=(mapping_tsv, summary_tsv, out_fa, mafft, pid_file),
    )
    proc.start()

    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline and not pid_file.exists():
        time.sleep(0.05)

    assert pid_file.exists()

    try:
        os.kill(proc.pid, signal.SIGINT)
        time.sleep(0.05)
        os.kill(proc.pid, signal.SIGINT)
    except ProcessLookupError:
        pass

    start = time.time()
    proc.join(timeout=2.15)
    elapsed = time.time() - start
    if proc.is_alive():
        try:
            os.kill(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.join(0.5)
        raise AssertionError(
            f"child did not exit within 2.00s (+0.15s epsilon); elapsed={elapsed:.3f}s"
        )
    code = proc.exitcode
    _assert_interrupted_exit(code)
    for line in pid_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        pid = int(line.strip())
        deadline = time.monotonic() + 1.0
        while time.monotonic() < deadline and _pid_is_alive(pid):
            time.sleep(0.05)
        assert not _pid_is_alive(pid)
