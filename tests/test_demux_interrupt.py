from __future__ import annotations

import gzip
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from ipyrad2.demuxer.demux import Demux
from ipyrad2.demuxer.demux_pipeline import _demux_spool_dir


pytestmark = pytest.mark.skipif(os.name != "posix", reason="interrupt tests require POSIX")


def _write_fastq(path: Path, reads: list[str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for idx, read in enumerate(reads):
            out.write(f"@r{idx}\n{read}\n+\n{'I' * len(read)}\n")
    return path


def _write_many_inline_reads(path: Path, n_reads: int) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for idx in range(n_reads):
            barcode = "ACGT" if idx % 2 == 0 else "TGCA"
            payload = f"{idx % 10000:04d}AAAA"
            seq = f"{barcode}ATCGG{payload}"
            out.write(f"@r{idx}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


def _run_demux_until_interrupted(
    fastqs: list[Path],
    barcodes: Path,
    outdir: Path,
    pigz: bool,
    fake_pigz_dir: Path | None = None,
    pid_file: Path | None = None,
) -> None:
    if fake_pigz_dir is not None:
        os.environ["PATH"] = f"{fake_pigz_dir}{os.pathsep}{os.environ['PATH']}"
    if pid_file is not None:
        os.environ["IPYRAD2_TEST_PIGZ_PID_FILE"] = str(pid_file)

    tool = Demux(
        fastqs=fastqs,
        barcodes=barcodes,
        cutsite_1="ATCGG",
        cutsite_2=None,
        max_mismatch=0,
        cores=2,
        chunksize=1,
        merge_technical_replicates=False,
        outdir=outdir,
        i7=False,
        disable_infer_cutsite_motifs=True,
        max_reads=None,
        max_reads_kmer=100,
        log_level="WARNING",
        pigz=pigz,
    )

    try:
        tool.run()
    except KeyboardInterrupt as exc:
        raise SystemExit(130) from exc


def _write_interrupt_driver(path: Path) -> Path:
    """Write a standalone Python driver for subprocess-based demux interrupt tests."""
    path.write_text(
        "from pathlib import Path\n"
        "import os\n"
        "import sys\n"
        "from ipyrad2.demuxer.demux import Demux\n"
        "\n"
        "raw = Path(sys.argv[1])\n"
        "barcodes = Path(sys.argv[2])\n"
        "outdir = Path(sys.argv[3])\n"
        "pigz = sys.argv[4] == '1'\n"
        "fake_pigz_dir = sys.argv[5] if len(sys.argv) > 5 else ''\n"
        "pid_file = sys.argv[6] if len(sys.argv) > 6 else ''\n"
        "if fake_pigz_dir:\n"
        "    os.environ['PATH'] = fake_pigz_dir + os.pathsep + os.environ['PATH']\n"
        "if pid_file:\n"
        "    os.environ['IPYRAD2_TEST_PIGZ_PID_FILE'] = pid_file\n"
        "tool = Demux(\n"
        "    fastqs=[raw],\n"
        "    barcodes=barcodes,\n"
        "    cutsite_1='ATCGG',\n"
        "    cutsite_2=None,\n"
        "    max_mismatch=0,\n"
        "    cores=2,\n"
        "    chunksize=1,\n"
        "    merge_technical_replicates=False,\n"
        "    outdir=outdir,\n"
        "    i7=False,\n"
        "    disable_infer_cutsite_motifs=True,\n"
        "    max_reads=None,\n"
        "    max_reads_kmer=100,\n"
        "    log_level='WARNING',\n"
        "    pigz=pigz,\n"
        ")\n"
        "try:\n"
        "    tool.run()\n"
        "except KeyboardInterrupt as exc:\n"
        "    raise SystemExit(130) from exc\n",
        encoding="utf-8",
    )
    return path


def _spawn_and_sigint(
    target,
    *args,
    grace: float = 3.0,
    start_delay: float = 0.25,
    repeat_delay: float = 0.20,
    epsilon: float = 0.25,
) -> int:
    ctx = mp.get_context("spawn")
    proc = ctx.Process(target=target, args=args)
    proc.start()

    time.sleep(start_delay)

    try:
        os.kill(proc.pid, signal.SIGINT)
        time.sleep(repeat_delay)
        if proc.is_alive():
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


def _run_subprocess_and_sigint(
    cmd: list[str],
    *,
    grace: float = 3.0,
    start_delay: float = 0.25,
    repeat_delay: float = 0.20,
    epsilon: float = 0.25,
) -> int:
    """Run a Python subprocess, interrupt it, and return its exit code."""
    proc = subprocess.Popen(cmd)
    time.sleep(start_delay)

    try:
        proc.send_signal(signal.SIGINT)
        time.sleep(repeat_delay)
        if proc.poll() is None:
            proc.send_signal(signal.SIGINT)
    except ProcessLookupError:
        pass

    start = time.time()
    try:
        return proc.wait(timeout=grace + epsilon)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=0.5)
        elapsed = time.time() - start
        raise AssertionError(
            f"child did not exit within {grace:.2f}s (+{epsilon:.2f}s epsilon); "
            f"elapsed={elapsed:.3f}s"
        ) from None


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


def _tracked_mp_child_pids() -> set[int]:
    """Return current spawn/resource-tracker PIDs used by multiprocessing."""
    pids: set[int] = set()
    for pattern in (
        "from multiprocessing.resource_tracker import main;main",
        "from multiprocessing.spawn import spawn_main",
    ):
        res = subprocess.run(
            ["pgrep", "-f", pattern],
            capture_output=True,
            text=True,
            check=False,
        )
        if res.returncode not in (0, 1):
            raise AssertionError(f"pgrep failed for pattern {pattern!r}: {res.stderr.strip()}")
        for line in res.stdout.splitlines():
            line = line.strip()
            if line:
                pids.add(int(line))
    return pids


def _describe_pids(pids: set[int]) -> str:
    """Return a compact `ps` summary for leftover multiprocessing children."""
    if not pids:
        return ""
    res = subprocess.run(
        ["ps", "-o", "pid,ppid,stat,etime,cmd", "-p", ",".join(str(pid) for pid in sorted(pids))],
        capture_output=True,
        text=True,
        check=False,
    )
    return res.stdout.strip() or ",".join(str(pid) for pid in sorted(pids))


def _assert_no_new_mp_children(before: set[int], settle_seconds: float = 5.0) -> None:
    """Assert that no new multiprocessing helper processes survive after a test run."""
    deadline = time.monotonic() + settle_seconds
    while time.monotonic() < deadline:
        leaked = _tracked_mp_child_pids() - before
        if not leaked:
            return
        time.sleep(0.05)
    leaked = _tracked_mp_child_pids() - before
    assert not leaked, f"leftover multiprocessing children detected:\n{_describe_pids(leaked)}"


def test_fast_shutdown_demux_pipeline(tmp_path: Path) -> None:
    before = _tracked_mp_child_pids()
    raw = _write_many_inline_reads(tmp_path / "lane.fastq.gz", n_reads=50_000)
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    code = _spawn_and_sigint(
        _run_demux_until_interrupted,
        [raw],
        barcodes,
        tmp_path / "out",
        False,
        None,
        None,
        grace=3.0,
        start_delay=0.20,
    )
    _assert_interrupted_exit(code)
    assert not _demux_spool_dir(tmp_path / "out").exists()


def test_fast_shutdown_demux_pigz_stage(tmp_path: Path) -> None:
    raw = _write_fastq(
        tmp_path / "lane.fastq.gz",
        ["ACGTATCGGAAAA", "TGCAATCGGCCCC", "ACGTATCGGTTTT", "TGCAATCGGGGGG"],
    )
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

    fake_pigz_dir = tmp_path / "bin"
    fake_pigz_dir.mkdir()
    pid_file = tmp_path / "pigz_pids.txt"
    fake_pigz = fake_pigz_dir / "pigz"
    fake_pigz.write_text(
        "#!/usr/bin/env bash\n"
        "echo $$ >> \"$IPYRAD2_TEST_PIGZ_PID_FILE\"\n"
        "infile=\"${@: -1}\"\n"
        "sleep 10\n"
        "cat \"$infile\"\n",
        encoding="utf-8",
    )
    fake_pigz.chmod(0o755)

    code = _spawn_and_sigint(
        _run_demux_until_interrupted,
        [raw],
        barcodes,
        tmp_path / "out",
        True,
        fake_pigz_dir,
        pid_file,
        grace=3.0,
        start_delay=0.60,
    )
    _assert_interrupted_exit(code)

    if pid_file.exists():
        for line in pid_file.read_text(encoding="utf-8").splitlines():
            pid = int(line.strip())
            deadline = time.monotonic() + 1.0
            while time.monotonic() < deadline and _pid_is_alive(pid):
                time.sleep(0.05)
            assert not _pid_is_alive(pid)


def test_repeated_fast_shutdown_demux_pipeline_subprocess_does_not_grow_mp_children(tmp_path: Path) -> None:
    before = _tracked_mp_child_pids()
    driver = _write_interrupt_driver(tmp_path / "run_demux_interrupt.py")
    for idx in range(2):
        case_dir = tmp_path / f"repeat_{idx}"
        case_dir.mkdir()
        raw = _write_many_inline_reads(case_dir / "lane.fastq.gz", n_reads=20_000)
        barcodes = case_dir / "barcodes.tsv"
        barcodes.write_text("sample1 ACGT\nsample2 TGCA\n", encoding="utf-8")

        code = _run_subprocess_and_sigint(
            [
                sys.executable,
                str(driver),
                str(raw),
                str(barcodes),
                str(case_dir / "out"),
                "0",
                "",
                "",
            ],
            grace=3.0,
            start_delay=0.20,
        )
        _assert_interrupted_exit(code)
        assert not _demux_spool_dir(case_dir / "out").exists()
        _assert_no_new_mp_children(before)
