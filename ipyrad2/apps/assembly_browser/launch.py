"""Launch helpers for the Streamlit assembly browser."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys

from ipyrad2.utils.exceptions import IPyradError


def _require_streamlit() -> None:
    """Raise a clear error if Streamlit is not installed."""
    if importlib.util.find_spec("streamlit") is not None:
        return
    raise IPyradError(
        "The inspect command requires streamlit. Install it with:\n"
        "  conda install -c conda-forge streamlit plotly\n"
        "or:\n"
        "  pip install 'ipyrad2[inspect]'"
    )


def launch_assembly_browser(assembly_dir: Path) -> None:
    """Launch the Streamlit browser for one assembly output directory."""
    assembly_dir = Path(assembly_dir).expanduser().resolve()
    if not assembly_dir.exists():
        raise IPyradError(f"assembly output directory does not exist: {assembly_dir}")
    if not assembly_dir.is_dir():
        raise IPyradError(f"assembly output path is not a directory: {assembly_dir}")

    _require_streamlit()
    app_path = Path(__file__).with_name("app.py")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(app_path),
            "--",
            "--assembly-dir",
            str(assembly_dir),
        ],
        check=True,
    )
