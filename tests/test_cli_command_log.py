from types import SimpleNamespace

from loguru import logger

import ipyrad2.cli.cli_main as cli_main
from ipyrad2.cli.command_log import MAX_LOGGED_MATCHED_PATH_CHARS
from ipyrad2.cli.command_log import format_logged_command


def _fit_paths(paths: list[str]) -> list[str]:
    included: list[str] = []
    for path in paths:
        candidate = " ".join(included + [path])
        if len(candidate) <= MAX_LOGGED_MATCHED_PATH_CHARS or not included:
            included.append(path)
            continue
        break
    return included


def test_format_logged_command_leaves_short_fastq_args_unchanged() -> None:
    argv = ["trim", "-d", "a.fastq.gz", "b.fastq.gz", "-o", "OUT"]

    result = format_logged_command(argv)

    assert result == "ipyrad2 trim -d a.fastq.gz b.fastq.gz -o OUT"


def test_format_logged_command_truncates_long_fastq_list_and_keeps_other_args() -> None:
    paths = [f"C0_DEMUX/sample_{idx:03d}_R1.fastq.gz" for idx in range(20)]
    included = _fit_paths(paths)

    result = format_logged_command(["trim", "-d", *paths, "-o", "OUT", "-q", "20"])

    assert result.startswith(f"ipyrad2 trim -d {' '.join(included)} ")
    assert f"...[truncated; {len(paths)} total matched paths]" in result
    assert result.endswith("-o OUT -q 20")


def test_format_logged_command_truncates_assemble_rad_and_wgs_lists_independently() -> None:
    rad_paths = [f"RAD/sample_{idx:03d}.bam" for idx in range(30)]
    wgs_paths = [f"WGS/sample_{idx:03d}.bam" for idx in range(30)]

    result = format_logged_command(
        ["assemble", "-d", *rad_paths, "-w", *wgs_paths, "-r", "REF.fa", "-o", "OUT"]
    )

    assert f"[truncated; {len(rad_paths)} total matched paths]" in result
    assert f"[truncated; {len(wgs_paths)} total matched paths]" in result
    assert result.count("total matched paths]") == 2
    assert result.endswith("-r REF.fa -o OUT")


def test_format_logged_command_does_not_truncate_analysis_data_arg() -> None:
    argv = ["wex", "-d", "seqs.hdf5", "-o", "OUT", "-n", "TEST"]

    result = format_logged_command(argv)

    assert result == "ipyrad2 wex -d seqs.hdf5 -o OUT -n TEST"


def test_run_subcommand_logs_truncated_command_for_trim(monkeypatch) -> None:
    fastqs = [f"C0_DEMUX/sample_{idx:03d}_R1.fastq.gz" for idx in range(20)]
    argv = ["trim", "-d", *fastqs, "-o", "OUT"]
    args = cli_main.setup_parsers().parse_args(argv)

    monkeypatch.setattr(
        cli_main.importlib,
        "import_module",
        lambda name, package=None: SimpleNamespace(run_trimmer=lambda **kwargs: None),
    )
    monkeypatch.setattr(cli_main.sys, "argv", ["ipyrad2", *argv])

    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}")
    try:
        cli_main.run_subcommand(args, _exit=False)
    finally:
        logger.remove(sink_id)

    cmd_message = next(str(msg).strip() for msg in messages if str(msg).startswith("CMD: "))
    assert f"[truncated; {len(fastqs)} total matched paths]" in cmd_message
    assert cmd_message.endswith("-o OUT")
