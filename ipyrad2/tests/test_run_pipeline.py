# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
import time
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ipyrad2.utils.parallel import run_pipeline


def _py(*lines: str) -> list[str]:
    """Build a Python -c one-liner as argv."""
    code = "; ".join(lines)
    return [sys.executable, "-c", code]


class TestRunPipeline(unittest.TestCase):
    def test_stdout_basic(self) -> None:
        # Stage1 prints "hello", Stage2 uppercases it.
        s1 = _py('import sys; sys.stdout.write("hello")')
        s2 = _py('import sys; sys.stdout.write(sys.stdin.read().upper())')
        rc, out, err = run_pipeline([s1, s2])
        self.assertEqual(rc, 0)
        self.assertEqual(out, b"HELLO")
        self.assertEqual(err, b"")

    def test_stdin_text_roundtrip(self) -> None:
        payload = "XYZ-123\n"
        cat = _py("import sys; sys.stdout.write(sys.stdin.read())")
        rc, out, err = run_pipeline([cat], stdin_text=payload)
        self.assertEqual(rc, 0)
        self.assertEqual(out, payload.encode("utf-8"))
        self.assertEqual(err, b"")

    def test_outfile_write_and_empty_out(self) -> None:
        with TemporaryDirectory() as td:
            outpath = Path(td) / "pipe.out"
            s1 = _py('import sys; sys.stdout.write("abcDEF")')
            s2 = _py("import sys; sys.stdout.write(sys.stdin.read())")
            rc, out, err = run_pipeline([s1, s2], outfile=outpath)
            self.assertEqual(rc, 0)
            self.assertEqual(out, b"")  # when outfile is used, out is empty
            self.assertEqual(err, b"")
            self.assertTrue(outpath.exists())
            self.assertEqual(outpath.read_bytes(), b"abcDEF")

    def test_large_stderr_does_not_block(self) -> None:
        # Stage1 spews ~1MB to stderr, then prints "ok" on stdout.
        s1 = _py(
            "import sys",
            "sys.stderr.write('X'*1024*1024)",
            "sys.stdout.write('ok')",
        )
        # Stage2 just passes through
        s2 = _py("import sys; sys.stdout.write(sys.stdin.read())")
        rc, out, err = run_pipeline([s1, s2])
        self.assertEqual(rc, 0)
        self.assertEqual(out, b"ok")
        self.assertGreaterEqual(len(err), 1024 * 1024)

    def test_nonzero_returncode_raises_and_kills_earlier_stages_fast(self) -> None:
        # Stage1 pretends to be long-running (sleep 10s).
        sleeper = _py("import time; time.sleep(10)")
        # Stage2 fails immediately.
        fail_fast = _py("import sys; sys.exit(3)")

        t0 = time.monotonic()
        with self.assertRaises(RuntimeError) as ctx:
            run_pipeline([sleeper, fail_fast])
        elapsed = time.monotonic() - t0

        # Should bail quickly (long before 10s)
        self.assertLess(elapsed, 2.0, f"pipeline did not abort fast enough (elapsed={elapsed:.3f}s)")
        # Error message should include rc and last stage argv
        msg = str(ctx.exception)
        self.assertIn("rc=3", msg)
        self.assertIn("sys.exit(3)", msg)

    def test_large_stdout(self) -> None:
        # Emit ~5MB on stdout; ensure it’s captured fully.
        s1 = _py("import sys", "sys.stdout.write('A'* (5*1024*1024))")
        rc, out, err = run_pipeline([s1])
        self.assertEqual(rc, 0)
        self.assertEqual(len(out), 5 * 1024 * 1024)
        self.assertEqual(err, b"")

    def test_binary_bytes_passthrough(self) -> None:
        # Write raw bytes (including NULs); confirm exact bytes are captured.
        s1 = _py(
            "import sys",
            "data=bytes([0,1,2,255])+b'ABC\\x00DEF'",
            "sys.stdout.buffer.write(data)",
        )
        rc, out, err = run_pipeline([s1])
        self.assertEqual(rc, 0)
        self.assertEqual(out, bytes([0, 1, 2, 255]) + b"ABC\x00DEF")
        self.assertEqual(err, b"")

    def test_stderr_warning_only(self) -> None:
        # Only stderr output; stdout empty; rc=0 should still succeed.
        s1 = _py("import sys; sys.stderr.write('warn\\n')")
        rc, out, err = run_pipeline([s1])
        self.assertEqual(rc, 0)
        self.assertEqual(out, b"")
        self.assertIn(b"warn", err)

    def test_pipeline_transform(self) -> None:
        # Stage1 prints numbers with commas; Stage2 strips commas.
        s1 = _py('import sys; sys.stdout.write("1,234,567")')
        s2 = _py("import sys; sys.stdout.write(sys.stdin.read().replace(',', ''))")
        rc, out, err = run_pipeline([s1, s2])
        self.assertEqual(rc, 0)
        self.assertEqual(out, b"1234567")
        self.assertEqual(err, b"")


if __name__ == "__main__":
    # Allow running directly: `python -m ipyrad2.tests.test_run_pipeline -q`
    unittest.main(verbosity=2)
