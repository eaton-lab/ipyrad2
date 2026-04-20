import re
import time

import pytest

import ipyrad2.utils.logger as logger_mod
from ipyrad2.utils.logger import set_log_level
from ipyrad2.utils.progress import ProgressBar


def test_progress_bar_uses_log_prefix_and_caller_file(capsys, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.stderr.isatty", lambda: True)
    monkeypatch.setattr(logger_mod, "color_support", lambda: False)
    set_log_level("INFO")
    prog = ProgressBar(4, start=time.time() - 5, message="Processing jobs")
    prog.finished = 1
    prog.update()
    prog.finished = 2
    prog.update()
    prog.close()

    err = capsys.readouterr().err
    assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| INFO\s+\|", err)
    assert "test_progress_logging.py" in err
    assert "Processing jobs" in err
    assert "% | Processing jobs" in err
    assert "0:00:05" not in err
    assert not re.search(r"\]\s+\d+%\s+\d+:\d{2}:\d{2}\s+\|", err)
    assert "\r" in err
    assert err.count("\n") == 1


def test_progress_bar_is_hidden_when_info_level_is_disabled(capsys, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.stderr.isatty", lambda: True)
    monkeypatch.setattr(logger_mod, "color_support", lambda: False)
    set_log_level("WARNING")
    prog = ProgressBar(2, start=time.time(), message="Hidden jobs")
    prog.finished = 1
    prog.update()
    prog.close()

    assert capsys.readouterr().err == ""


def test_progress_bar_deduplicates_identical_renders(capsys, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.stderr.isatty", lambda: True)
    monkeypatch.setattr(logger_mod, "color_support", lambda: False)
    set_log_level("INFO")
    prog = ProgressBar(200, start=time.time(), message="Busy jobs")
    prog.finished = 1
    prog.update()
    prog.finished = 1
    prog.update()
    prog.finished = 2
    prog.update()
    prog.close()

    err = capsys.readouterr().err
    assert err.count("  0% | Busy jobs") == 1
    assert err.count("  1% | Busy jobs") == 1
    assert err.count("\r") == 2
    assert err.count("\n") == 1


def test_progress_bar_writes_each_update_on_new_line_when_not_tty(
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)
    monkeypatch.setattr(logger_mod, "color_support", lambda: False)
    set_log_level("INFO")
    prog = ProgressBar(4, start=time.time() - 5, message="File jobs")
    prog.finished = 1
    prog.update()
    prog.finished = 2
    prog.update()
    prog.close()

    err = capsys.readouterr().err
    assert "\r" not in err
    assert err.count(" 25% | File jobs") == 1
    assert err.count(" 50% | File jobs") == 1
    assert err.count("\n") == 2


def test_progress_bar_close_does_not_add_blank_line_when_not_tty(
    capsys,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)
    monkeypatch.setattr(logger_mod, "color_support", lambda: False)
    set_log_level("INFO")
    prog = ProgressBar(2, start=time.time(), message="Non-tty close")
    prog.finished = 1
    prog.update()
    prog.close()

    err = capsys.readouterr().err
    assert err.count("Non-tty close") == 1
    assert err.count("\n") == 1
