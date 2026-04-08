from __future__ import annotations

import gzip
import subprocess
import sys
import time
from pathlib import Path

import pytest

import ipyrad2.utils.parallel.pool as parallel_pool_mod
from ipyrad2.demuxer.match import BarMatching
from ipyrad2.utils.parallel import ParallelJobError
from ipyrad2.utils.parallel import run_pipeline
from ipyrad2.utils.parallel import run_with_pool
from ipyrad2.utils.parallel import run_with_pool_iter


REPO_ROOT = Path(__file__).resolve().parents[1]


def _add(a: int, b: int) -> int:
    return a + b


def _fail_worker(label: str) -> str:
    raise ValueError(f"boom:{label}")


def _sleep_pipeline_worker(seconds: float) -> float:
    run_pipeline([["bash", "-lc", f"sleep {seconds}"]])
    return seconds


def _write_chunk_fail(path: Path, data: list[bytes]) -> None:
    raise RuntimeError(f"intentional write failure for {Path(path).name}")


class _OneChunkMatcher(BarMatching):
    def _format_check(self) -> None:
        pass

    def _iter_matched_barcode(self):
        return iter(())

    def _iter_matched_chunks(self):
        yield {"sample": [b"@r1\nAAAA\n+\n!!!!\n"]}, {"sample": [b"@r2\nTTTT\n+\n!!!!\n"]}


class _FailingChunkMatcher(_OneChunkMatcher):
    def _build_write_jobs(self, read1s, read2s):
        return {
            "sample_R1": (
                _write_chunk_fail,
                {"path": self.outdir / "sample_R1.fastq.gz", "data": read1s["sample"]},
            ),
        }


def test_run_with_pool_empty_jobs_returns_immediately(monkeypatch: pytest.MonkeyPatch) -> None:
    class _SentinelExecutor:
        def __init__(self, *args, **kwargs) -> None:
            raise AssertionError("executor should not be created for empty jobs")

    monkeypatch.setattr(parallel_pool_mod, "ProcessPoolExecutor", _SentinelExecutor)

    assert run_with_pool({}, "WARNING") == {}


def test_run_with_pool_iter_empty_iterator_returns_immediately(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _SentinelExecutor:
        def __init__(self, *args, **kwargs) -> None:
            raise AssertionError("executor should not be created for empty iterators")

    monkeypatch.setattr(parallel_pool_mod, "ProcessPoolExecutor", _SentinelExecutor)

    assert list(run_with_pool_iter(iter(()), "WARNING")) == []


@pytest.mark.parametrize("name,value", [("max_workers", 0), ("max_workers", False), ("max_inflight", 0)])
def test_run_with_pool_validates_positive_ints(name: str, value: object) -> None:
    jobs = {"ok": (_add, {"a": 1, "b": 2})}

    with pytest.raises(ValueError):
        run_with_pool(jobs, "WARNING", **{name: value})


def test_run_with_pool_iter_validates_positive_ints() -> None:
    jobs = [("ok", (_add, {"a": 1, "b": 2}))]

    with pytest.raises(ValueError):
        list(run_with_pool_iter(jobs, "WARNING", max_workers=0))


def test_run_with_pool_iter_requires_msg_and_njobs_together() -> None:
    jobs = [("ok", (_add, {"a": 1, "b": 2}))]

    with pytest.raises(ValueError):
        list(run_with_pool_iter(jobs, "WARNING", msg="Processing"))

    with pytest.raises(ValueError):
        list(run_with_pool_iter(jobs, "WARNING", njobs=1))


def test_run_with_pool_returns_results() -> None:
    jobs = {
        "a": (_add, {"a": 1, "b": 2}),
        "b": (_add, {"a": 3, "b": 4}),
    }

    assert run_with_pool(jobs, "WARNING", max_workers=2) == {"a": 3, "b": 7}


def test_run_with_pool_iter_progress_uses_custom_increment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[tuple[str, int, str] | tuple[str, int]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, message))

        def update(self) -> None:
            events.append(("update", self.finished))

        def close(self) -> None:
            events.append(("close", self.finished))

    monkeypatch.setattr(parallel_pool_mod, "ProgressBar", _ProgressStub)

    jobs = [
        ("a", (_add, {"a": 1, "b": 2})),
        ("b", (_add, {"a": 3, "b": 4})),
    ]

    results = list(
        run_with_pool_iter(
            jobs,
            "WARNING",
            max_workers=1,
            msg="Aligning loci",
            njobs=5,
            progress_increment=lambda key, _result: 2 if key == "a" else 3,
        )
    )

    assert results == [("a", 3), ("b", 7)]
    assert ("init", 5, "Aligning loci - total jobs: 5") in events
    assert ("update", 0) in events
    assert ("update", 2) in events
    assert ("update", 5) in events
    assert events[-1] == ("close", 5)


def test_run_with_pool_iter_progress_closes_on_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[tuple[str, int, str] | tuple[str, int]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, message))

        def update(self) -> None:
            events.append(("update", self.finished))

        def close(self) -> None:
            events.append(("close", self.finished))

    monkeypatch.setattr(parallel_pool_mod, "ProgressBar", _ProgressStub)

    jobs = [
        ("ok", (_add, {"a": 1, "b": 2})),
        ("bad", (_fail_worker, {"label": "sample"})),
    ]

    with pytest.raises(ParallelJobError):
        list(
            run_with_pool_iter(
                jobs,
                "WARNING",
                max_workers=1,
                msg="Aligning loci",
                njobs=2,
            )
        )

    assert ("init", 2, "Aligning loci - total jobs: 2") in events
    assert ("update", 0) in events
    assert events[-1][0] == "close"


def test_run_with_pool_raises_parallel_job_error_with_metadata() -> None:
    jobs = {
        "ok": (_add, {"a": 1, "b": 2}),
        "bad": (_fail_worker, {"label": "sample"}),
    }

    with pytest.raises(ParallelJobError) as excinfo:
        run_with_pool(jobs, "WARNING", max_workers=2)

    exc = excinfo.value
    assert exc.key == "bad"
    assert exc.func_name == "_fail_worker"
    assert exc.kwargs == {"label": "sample"}
    assert isinstance(exc.original, ValueError)


def test_run_with_pool_kills_long_running_pipeline_worker_on_sibling_failure() -> None:
    jobs = {
        "slow": (_sleep_pipeline_worker, {"seconds": 10.0}),
        "bad": (_fail_worker, {"label": "sample"}),
    }

    start = time.monotonic()
    with pytest.raises(ParallelJobError):
        run_with_pool(jobs, "WARNING", max_workers=2)
    elapsed = time.monotonic() - start

    assert elapsed < 3.0


def test_failure_path_does_not_emit_resource_tracker_or_afterfork(tmp_path: Path) -> None:
    script = tmp_path / "parallel_failure_script.py"
    script.write_text(
        "from ipyrad2.utils.parallel import run_with_pool, run_pipeline\n"
        "\n"
        "def slow(seconds):\n"
        "    run_pipeline([[\"bash\", \"-lc\", f\"sleep {seconds}\"]])\n"
        "    return seconds\n"
        "\n"
        "def fail():\n"
        "    raise ValueError(\"boom\")\n"
        "\n"
        "if __name__ == \"__main__\":\n"
        "    jobs = {\n"
        "        \"slow\": (slow, {\"seconds\": 10.0}),\n"
        "        \"bad\": (fail, {}),\n"
        "    }\n"
        "    try:\n"
        "        run_with_pool(jobs, \"WARNING\", max_workers=2)\n"
        "    except Exception:\n"
        "        pass\n",
        encoding="utf-8",
    )

    res = subprocess.run(
        [sys.executable, str(script)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert res.returncode == 0
    assert "_afterFork" not in res.stderr
    assert "resource_tracker" not in res.stderr


def test_bar_matching_run_writes_chunk_outputs(tmp_path: Path) -> None:
    matcher = _OneChunkMatcher(
        fastqs=(tmp_path / "in_R1.fastq.gz", None),
        barcodes_to_names={},
        barcode_lengths1=(),
        barcode_lengths2=(),
        cuts1=[],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path,
        log_level="WARNING",
        workers=1,
        chunksize=10,
    )

    matcher.run()

    with gzip.open(tmp_path / "sample_R1.fastq.gz", "rb") as indata:
        assert indata.read() == b"@r1\nAAAA\n+\n!!!!\n"
    with gzip.open(tmp_path / "sample_R2.fastq.gz", "rb") as indata:
        assert indata.read() == b"@r2\nTTTT\n+\n!!!!\n"


def test_bar_matching_run_raises_parallel_job_error_on_write_failure(tmp_path: Path) -> None:
    matcher = _FailingChunkMatcher(
        fastqs=(tmp_path / "in_R1.fastq.gz", None),
        barcodes_to_names={},
        barcode_lengths1=(),
        barcode_lengths2=(),
        cuts1=[],
        cuts2=[],
        merge_technical_replicates=False,
        outdir=tmp_path,
        log_level="WARNING",
        workers=1,
        chunksize=10,
    )

    with pytest.raises(ParallelJobError) as excinfo:
        matcher.run()

    exc = excinfo.value
    assert exc.key == "sample_R1"
    assert exc.func_name == "_write_chunk_fail"
    assert "<list len=1>" in str(exc)
