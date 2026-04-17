import json
import subprocess
import sys
import textwrap
from typing import Set
from pathlib import Path
from types import SimpleNamespace

import ipyrad2.cli.cli_analysis as cli_analysis
import ipyrad2.cli.cli_main as cli_main


ROOT = Path(__file__).resolve().parents[1]
HEAVY_MODULES = (
    "ipyrad2.analysis.extracters.window_extracter",
    "ipyrad2.analysis.methods.bpp",
    "ipyrad2.analysis.methods.popgen.runner",
    "requests",
    "toytree",
)
UNRELATED_WEX_HELP_MODULES = (
    "ipyrad2.cli.cli_bpp",
    "ipyrad2.cli.cli_popgen",
    "ipyrad2.analysis.methods.popgen.models",
)


def _loaded_modules(script_body: str, module_names: tuple[str, ...]) -> Set[str]:
    code = "\n".join(
        [
            "import json",
            "import sys",
            textwrap.dedent(script_body).strip(),
            f"print(json.dumps(sorted(name for name in sys.modules if name in {module_names!r})))",
        ]
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    return set(json.loads(result.stdout))


def _loaded_heavy_modules(script_body: str) -> Set[str]:
    return _loaded_modules(script_body, HEAVY_MODULES)


def test_importing_cli_main_does_not_load_heavy_analysis_runtimes() -> None:
    loaded = _loaded_heavy_modules("import ipyrad2.cli.cli_main")
    assert loaded == set()


def test_building_all_parsers_keeps_heavy_analysis_runtimes_unloaded() -> None:
    loaded = _loaded_heavy_modules(
        "from ipyrad2.cli.cli_main import setup_parsers\nsetup_parsers()"
    )
    assert loaded == set()


def test_building_one_analysis_parser_keeps_other_heavy_runtimes_unloaded() -> None:
    loaded = _loaded_heavy_modules(
        "import contextlib\n"
        "import io\n"
        "from ipyrad2.cli.cli_main import command_line\n"
        "buf = io.StringIO()\n"
        "try:\n"
        "    with contextlib.redirect_stdout(buf):\n"
        "        command_line(['wex', '-h'])\n"
        "except SystemExit:\n"
        "    pass"
    )
    assert loaded == set()


def test_wex_help_imports_only_the_requested_analysis_parser_path() -> None:
    loaded = _loaded_modules(
        "import contextlib\n"
        "import io\n"
        "from ipyrad2.cli.cli_main import command_line\n"
        "buf = io.StringIO()\n"
        "try:\n"
        "    with contextlib.redirect_stdout(buf):\n"
        "        command_line(['wex', '-h'])\n"
        "except SystemExit:\n"
        "    pass",
        UNRELATED_WEX_HELP_MODULES,
    )
    assert loaded == set()


def test_run_subcommand_lazily_imports_analysis_runner(monkeypatch) -> None:
    args = cli_main.setup_parsers().parse_args(["wex", "-d", "assembly.hdf5"])
    calls: list[dict] = []
    real_import_module = cli_analysis.importlib.import_module

    def fake_import_module(name, package=None):
        if name == "..analysis.extracters.window_extracter":
            return SimpleNamespace(
                run_window_extracter=lambda **kwargs: calls.append(kwargs)
            )
        return real_import_module(name, package)

    monkeypatch.setattr(cli_analysis.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(cli_main.sys, "argv", ["ipyrad2", "wex", "-d", "assembly.hdf5"])

    cli_main.run_subcommand(args, _exit=False)

    assert len(calls) == 1
    assert calls[0]["data"] == Path("assembly.hdf5")
