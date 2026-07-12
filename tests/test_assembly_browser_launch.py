from pathlib import Path
import subprocess
import sys

import pytest

from ipyrad2.apps.assembly_browser import launch
from ipyrad2.utils.exceptions import IPyradError


def test_launch_rejects_missing_streamlit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(launch.importlib.util, "find_spec", lambda name: None)

    with pytest.raises(IPyradError, match="requires streamlit"):
        launch.launch_assembly_browser(tmp_path)


def test_launch_rejects_missing_directory() -> None:
    with pytest.raises(IPyradError, match="assembly output directory does not exist"):
        launch.launch_assembly_browser(Path("definitely-missing-assembly-dir"))


def test_launch_runs_streamlit_with_assembly_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = []
    monkeypatch.setattr(launch.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(
        launch.subprocess,
        "run",
        lambda cmd, check: calls.append((cmd, check)),
    )

    launch.launch_assembly_browser(tmp_path)

    assert len(calls) == 1
    cmd, check = calls[0]
    assert check is True
    assert cmd[:4] == [sys.executable, "-m", "streamlit", "run"]
    assert cmd[-2:] == ["--assembly-dir", str(tmp_path.resolve())]
    assert cmd[-3] == "--"
    assert cmd[4].endswith("ipyrad2/apps/assembly_browser/app.py")


def test_launch_propagates_streamlit_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(_cmd, check):
        raise subprocess.CalledProcessError(1, _cmd)

    monkeypatch.setattr(launch.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(launch.subprocess, "run", fail)

    with pytest.raises(subprocess.CalledProcessError):
        launch.launch_assembly_browser(tmp_path)
