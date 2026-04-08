import re
import time

from ipyrad2.utils.logger import set_log_level
from ipyrad2.utils.progress import ProgressBar


def test_progress_bar_uses_log_prefix_and_caller_file(capsys) -> None:
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


def test_progress_bar_is_hidden_when_info_level_is_disabled(capsys) -> None:
    set_log_level("WARNING")
    prog = ProgressBar(2, start=time.time(), message="Hidden jobs")
    prog.finished = 1
    prog.update()
    prog.close()

    assert capsys.readouterr().err == ""


def test_progress_bar_deduplicates_identical_renders(capsys) -> None:
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
