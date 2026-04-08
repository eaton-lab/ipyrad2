from __future__ import annotations

import subprocess
import sys

import ipyrad2


def test_package_exposes_version_string() -> None:
    assert isinstance(ipyrad2.__version__, str)
    assert ipyrad2.__version__


def test_cli_reports_same_version_string() -> None:
    res = subprocess.run(
        [sys.executable, "-m", "ipyrad2.cli.cli_main", "-v"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert res.returncode == 0
    assert ipyrad2.__version__ in res.stdout
