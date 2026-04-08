from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

from ipyrad2.utils.parallel import stream_pipeline_lines


def test_stream_pipeline_lines_yields_output_incrementally() -> None:
    cmd = [
        sys.executable,
        "-c",
        (
            "import sys, time; "
            "print('first', flush=True); "
            "time.sleep(0.35); "
            "print('second', flush=True)"
        ),
    ]

    start = time.monotonic()
    lines = stream_pipeline_lines([cmd])
    first = next(lines)
    elapsed = time.monotonic() - start
    second = next(lines)

    assert first == "first"
    assert second == "second"
    assert elapsed < 0.25


def test_stream_pipeline_lines_raises_on_pipeline_failure() -> None:
    cmd = [
        sys.executable,
        "-c",
        "import sys; print('ok'); sys.exit(3)",
    ]

    with pytest.raises(RuntimeError, match="pipeline failed"):
        list(stream_pipeline_lines([cmd]))


def test_stream_pipeline_lines_accepts_stdin_text() -> None:
    cmd = [
        sys.executable,
        "-c",
        "import sys; sys.stdout.write(sys.stdin.read().upper())",
    ]

    assert list(stream_pipeline_lines([cmd], stdin_text="abc")) == ["ABC"]


def test_stream_pipeline_lines_returns_final_line_without_newline() -> None:
    cmd = [
        sys.executable,
        "-c",
        "import sys; sys.stdout.write('line-without-newline')",
    ]

    assert list(stream_pipeline_lines([cmd])) == ["line-without-newline"]


def test_stream_pipeline_lines_supports_multistage_transform() -> None:
    cmd1 = [
        sys.executable,
        "-c",
        "import sys; sys.stdout.write('1,234,567\\n')",
    ]
    cmd2 = [
        sys.executable,
        "-c",
        "import sys; sys.stdout.write(sys.stdin.read().replace(',', ''))",
    ]

    assert list(stream_pipeline_lines([cmd1, cmd2])) == ["1234567"]


def test_stream_pipeline_lines_failure_includes_stderr() -> None:
    cmd = [
        sys.executable,
        "-c",
        "import sys; sys.stderr.write('boom\\n'); sys.exit(3)",
    ]

    with pytest.raises(RuntimeError, match="boom"):
        list(stream_pipeline_lines([cmd]))


def test_stream_pipeline_lines_close_terminates_child_pipeline(tmp_path: Path) -> None:
    marker = tmp_path / "done.txt"
    cmd = [
        sys.executable,
        "-c",
        (
            "import pathlib, time; "
            "print('first', flush=True); "
            f"time.sleep(1.0); pathlib.Path(r'{marker}').write_text('done')"
        ),
    ]

    lines = stream_pipeline_lines([cmd])
    assert next(lines) == "first"
    lines.close()
    time.sleep(1.2)

    assert not marker.exists()


def test_stream_pipeline_lines_context_manager_closes_early(tmp_path: Path) -> None:
    marker = tmp_path / "done.txt"
    cmd = [
        sys.executable,
        "-c",
        (
            "import pathlib, time; "
            "print('first', flush=True); "
            f"time.sleep(1.0); pathlib.Path(r'{marker}').write_text('done')"
        ),
    ]

    with stream_pipeline_lines([cmd]) as lines:
        assert next(lines) == "first"

    time.sleep(1.2)
    assert not marker.exists()
