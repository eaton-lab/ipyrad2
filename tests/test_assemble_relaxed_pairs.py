"""Restore the historical opt-in pair policy without relaxing read quality."""
import json
import subprocess
from types import SimpleNamespace

import pytest

from ipyrad2.assembler import read_filters as rf
from ipyrad2.cli import cli_main
from test_assemble_filter_diagnostics import fixture_bam, sam


@pytest.mark.parametrize("relaxed", [False, True])
@pytest.mark.parametrize("model", ["native", "sfs"])
def test_cli_dispatches_pair_policy_with_both_callers(monkeypatch, relaxed, model):
    argv = ["assemble", "-d", "rad.bam", "-w", "wgs.bam", "-r", "ref.fa",
            "--genotype-model", model]
    if relaxed:
        argv += ["--relaxed-pairs"]
    args = cli_main.setup_parsers().parse_args(argv)
    seen = {}
    monkeypatch.setattr(cli_main.importlib, "import_module", lambda *a, **k:
                        SimpleNamespace(run_assembler=lambda **kw: seen.update(kw)))
    monkeypatch.setattr(cli_main.sys, "argv", ["ipyrad2", *argv])
    cli_main.run_subcommand(args, _exit=False)
    assert seen["relaxed_pairs"] is relaxed
    assert seen["genotype_model"] == model


def test_relaxed_native_filter_retains_only_expected_additional_reads(tmp_path):
    bam, _ = fixture_bam(tmp_path)
    settings = dict(is_paired=True, min_map_q=30, max_tlen=1000,
                    max_softclip=5, max_nm=5, min_aligned_len=100)
    results = {}
    for relaxed in (False, True):
        result = rf.prepare_filtered_analysis_bam(
            sname="s", bam_file=bam, tmpdir=tmp_path / str(relaxed), threads=1,
            relaxed_pairs=relaxed, **settings)
        rows = subprocess.check_output([rf.BIN_SAM, "view", str(result.bam_path)], text=True)
        results[relaxed] = {l.split("\t")[0]: l for l in rows.splitlines()}
        diag = json.loads(result.bam_path.with_suffix(".filters.json").read_text())
        assert diag["counts_reconcile"]
        assert diag["settings"]["relaxed_pairs"] is relaxed
    strict = {"pass", "insertions", "secondary", "hardclip", "boundary"}
    assert set(results[False]) == strict
    assert set(results[True]) == strict | {"mate", "unmappedmate", "long"}
    assert all(results[False][s] == results[True][s] for s in strict)
    # Read flags are retained: this option does not unmate the reads.
    assert int(results[True]["unmappedmate"].split("\t")[1]) == 9


def test_unmapped_read_still_rejected_and_single_end_policy_unchanged(tmp_path):
    path = tmp_path / "unmapped.sam"
    lines = [sam("mapped", flag=9), sam("unmapped", flag=5)]
    path.write_text("@SQ\tSN:chr1\tLN:2000\n" + "".join(lines))
    filters = dict(max_tlen=None, max_softclip=None, max_nm=None, min_aligned_len=None)
    expr = rf.build_mapped_read_filter_expr(is_paired=True, relaxed_pairs=True, **filters)
    rows = subprocess.check_output([rf.BIN_SAM, "view", "-e", expr, str(path)], text=True)
    assert [l.split("\t")[0] for l in rows.splitlines()] == ["mapped"]
    diag = rf.summarize_read_filter_diagnostics(
        lines, min_map_q=0, is_paired=True, relaxed_pairs=True, **filters)
    assert diag["first_failure"]["unmapped"] == 1
    assert diag["predicted_retained_records"] == 1
    for filters in (filters, dict(max_tlen=1000, max_softclip=5, max_nm=5, min_aligned_len=100)):
        assert rf.build_mapped_read_filter_expr(is_paired=False, **filters) == \
            rf.build_mapped_read_filter_expr(is_paired=False, relaxed_pairs=True, **filters)
