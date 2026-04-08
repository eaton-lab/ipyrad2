from __future__ import annotations

import os
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pytest

from ipyrad2.assembler.paralogs import make_indel_mask_bed
from ipyrad2.assembler.variants import BIN_BCF
from ipyrad2.utils.parallel import pipeline as pipeline_module
from ipyrad2.utils.parallel import PipelineTimeoutError
from ipyrad2.utils.parallel import run_pipeline


def _py(*lines: str) -> list[str]:
    """Build a Python -c command from one or more statements."""
    return [sys.executable, "-c", "; ".join(lines)]


def test_run_pipeline_basic_stdout_transform() -> None:
    cmd1 = _py('import sys; sys.stdout.write("hello")')
    cmd2 = _py('import sys; sys.stdout.write(sys.stdin.read().upper())')

    rc, out, err = run_pipeline([cmd1, cmd2])

    assert rc == 0
    assert out == b"HELLO"
    assert err == b""


def test_run_pipeline_stdin_text_roundtrip() -> None:
    payload = "XYZ-123\n"
    cmd = _py("import sys; sys.stdout.write(sys.stdin.read())")

    rc, out, err = run_pipeline([cmd], stdin_text=payload)

    assert rc == 0
    assert out == payload.encode("utf-8")
    assert err == b""


def test_run_pipeline_outfile_success_returns_empty_stdout(tmp_path: Path) -> None:
    outpath = tmp_path / "pipe.out"
    cmd1 = _py('import sys; sys.stdout.write("abcDEF")')
    cmd2 = _py("import sys; sys.stdout.write(sys.stdin.read())")

    rc, out, err = run_pipeline([cmd1, cmd2], outfile=outpath)

    assert rc == 0
    assert out == b""
    assert err == b""
    assert outpath.read_bytes() == b"abcDEF"


def test_run_pipeline_outfile_failure_leaves_existing_file_untouched(tmp_path: Path) -> None:
    outpath = tmp_path / "pipe.out"
    outpath.write_bytes(b"original")
    cmd1 = _py('import sys; sys.stdout.write("ABCDEF")')
    cmd2 = _py(
        "import sys",
        "data = sys.stdin.read()",
        "sys.stdout.write(data[:3])",
        "sys.stderr.write('boom\\n')",
        "sys.exit(3)",
    )

    with pytest.raises(RuntimeError, match="boom"):
        run_pipeline([cmd1, cmd2], outfile=outpath)

    assert outpath.read_bytes() == b"original"
    assert list(tmp_path.glob(".pipe.out.*.tmp")) == []


def test_run_pipeline_large_stderr_does_not_block() -> None:
    cmd1 = _py(
        "import sys",
        "sys.stderr.write('X' * 1024 * 1024)",
        "sys.stdout.write('ok')",
    )
    cmd2 = _py("import sys; sys.stdout.write(sys.stdin.read())")

    rc, out, err = run_pipeline([cmd1, cmd2])

    assert rc == 0
    assert out == b"ok"
    assert len(err) >= 1024 * 1024


def test_run_pipeline_nonzero_returncode_fails_fast() -> None:
    sleeper = _py("import time; time.sleep(10)")
    fail_fast = _py("import sys; sys.exit(3)")

    start = time.monotonic()
    with pytest.raises(RuntimeError) as excinfo:
        run_pipeline([sleeper, fail_fast])
    elapsed = time.monotonic() - start

    assert elapsed < 2.0
    assert "rc=3" in str(excinfo.value)
    assert "sys.exit(3)" in str(excinfo.value)


def test_run_pipeline_timeout_raises_quickly() -> None:
    sleeper = _py("import time; time.sleep(10)")

    start = time.monotonic()
    with pytest.raises(PipelineTimeoutError) as excinfo:
        run_pipeline([sleeper], timeout_s=0.2)
    elapsed = time.monotonic() - start

    assert elapsed < 2.0
    assert "timed out after" in str(excinfo.value)
    assert "timeout=0.2s" in str(excinfo.value)


def test_run_pipeline_disk_full_error_is_summarized() -> None:
    noisy = _py(
        "import sys",
        "sys.stderr.write('bwa noise\\n' * 5000)",
        "sys.stdout.write('ok')",
    )
    disk_full = _py(
        "import sys",
        "sys.stdin.read()",
        "sys.stderr.write('[E::bgzf_flush] File write failed (wrong size)\\n')",
        "sys.stderr.write('samtools sort: failed writing to \"/tmp/sample.tmp.bam\": No space left on device\\n')",
        "sys.stderr.write('[E::bgzf_close] File write failed\\n')",
        "sys.exit(1)",
    )

    with pytest.raises(RuntimeError) as excinfo:
        run_pipeline([noisy, disk_full])

    message = str(excinfo.value)
    assert "disk space error:" in message
    assert "No space left on device" in message
    assert "/tmp/sample.tmp.bam" in message
    assert "bwa noise" not in message


def test_run_pipeline_binary_bytes_passthrough() -> None:
    cmd = _py(
        "import sys",
        "data = bytes([0, 1, 2, 255]) + b'ABC\\x00DEF'",
        "sys.stdout.buffer.write(data)",
    )

    rc, out, err = run_pipeline([cmd])

    assert rc == 0
    assert out == bytes([0, 1, 2, 255]) + b"ABC\x00DEF"
    assert err == b""


def test_safe_popen_uses_start_new_session(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    class _DummyProc:
        pid = os.getpid()
        stdin = None
        stdout = None
        stderr = None

        def poll(self) -> int:
            return 0

    def _fake_popen(argv, **kwargs):
        captured["argv"] = argv
        captured["kwargs"] = kwargs
        return _DummyProc()

    monkeypatch.setattr(pipeline_module.sp, "Popen", _fake_popen)
    monkeypatch.setattr(pipeline_module, "_CHILD_PG_QUEUE", None)

    proc = pipeline_module.safe_popen(["echo", "x"])

    assert proc is not None
    assert captured["argv"] == ["echo", "x"]
    assert captured["kwargs"]["start_new_session"] is True
    assert "preexec_fn" not in captured["kwargs"]

    with pipeline_module._CHILD_PROCS_LOCK:
        pipeline_module._CHILD_PROCS.clear()


def _write_stress_vcf(tmp_path: Path, stem: str, nsites: int, seed: int) -> Path:
    rng = random.Random(seed)
    plain = tmp_path / f"{stem}.vcf"
    with plain.open("w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write('##INFO=<ID=INDEL,Number=0,Type=Flag,Description="indel">\n')
        out.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">\n')
        out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">\n')
        out.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
        pos = 100
        for idx in range(nsites):
            pos += rng.randint(1, 4)
            if idx % 17 == 0:
                ref = "AT"
                alt = "A"
                info = "INDEL;DP=8"
                ad = "2,6"
            else:
                ref = rng.choice("ACGT")
                alt = rng.choice([base for base in "ACGT" if base != ref])
                info = "DP=8"
                ad = "4,4"
            out.write(
                f"chr1\t{pos}\t.\t{ref}\t{alt}\t60\t.\t{info}\tGT:DP:AD:GQ\t0/1:8:{ad}:50\n"
            )
    vcf_gz = tmp_path / f"{stem}.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(vcf_gz), str(plain)]])
    return vcf_gz


def test_make_indel_mask_bed_concurrent_stress(tmp_path: Path) -> None:
    vcfs = [_write_stress_vcf(tmp_path, f"s{i}", nsites=3000, seed=i) for i in range(24)]

    def _work(vcf: Path) -> Path:
        out_bed = tmp_path / f"{vcf.name}.bed"
        make_indel_mask_bed(vcf, out_bed, pad_bp=10)
        return out_bed

    outputs: list[Path] = []
    with ThreadPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(_work, vcf): vcf for vcf in vcfs}
        for fut in as_completed(futs):
            outputs.append(fut.result())

    assert len(outputs) == len(vcfs)
    assert all(path.exists() for path in outputs)
    assert all(path.read_text(encoding="utf-8").strip() for path in outputs)
