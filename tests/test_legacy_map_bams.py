from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from ipyrad2.mapper.legacy_map_bams import BIN_SAMTOOLS
from ipyrad2.mapper.legacy_map_bams import CURRENT_BAM_SUFFIX
from ipyrad2.mapper.legacy_map_bams import LEGACY_PLAIN_BAM_SUFFIX
from ipyrad2.mapper.legacy_map_bams import LegacyBamMigration
from ipyrad2.mapper.legacy_map_bams import MapStatsMigration
from ipyrad2.mapper.legacy_map_bams import _rewrite_bam_rg_tags
from ipyrad2.mapper.legacy_map_bams import _rewrite_stats_sidecar
from ipyrad2.mapper.legacy_map_bams import migrate_legacy_map_outputs
from ipyrad2.mapper.legacy_map_bams import migrated_bam_name
from ipyrad2.mapper.legacy_map_bams import plan_legacy_map_migration
from ipyrad2.mapper.legacy_map_bams import read_bam_sample_name
from ipyrad2.utils.exceptions import IPyradError


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "migrate_legacy_map_bams.py"
    spec = importlib.util.spec_from_file_location("migrate_legacy_map_bams_script", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_migrated_bam_name_updates_legacy_suffix_only() -> None:
    assert migrated_bam_name("sample.filtered.bam") == "sample.trimmed.sorted.bam"
    assert migrated_bam_name("sample.trimmed.filtered.bam") == "sample.trimmed.sorted.bam"
    assert migrated_bam_name("sample.trimmed.sorted.bam") == "sample.trimmed.sorted.bam"
    with pytest.raises(IPyradError):
        migrated_bam_name("sample.bam")


def test_read_bam_sample_name_uses_samtools_samples(monkeypatch, tmp_path: Path) -> None:
    bam_path = tmp_path / "sample.trimmed.filtered.bam"
    bam_path.write_text("", encoding="utf-8")
    recorded: list[list[list[str]]] = []

    def _fake_run_pipeline(cmds, **_kwargs):
        recorded.append(cmds)
        return 0, b"sample.trimmed\n", b""

    monkeypatch.setattr("ipyrad2.mapper.legacy_map_bams.run_pipeline", _fake_run_pipeline)

    assert read_bam_sample_name(bam_path) == "sample.trimmed"
    assert recorded == [[[BIN_SAMTOOLS, "samples", str(bam_path)]]]


def test_plan_legacy_map_migration_discovers_bams_stats_and_skips(monkeypatch, tmp_path: Path) -> None:
    indir = tmp_path / "indir"
    outdir = tmp_path / "outdir"
    indir.mkdir()
    (indir / "alpha.trimmed.filtered.bam").write_text("", encoding="utf-8")
    (indir / "beta.trimmed.sorted.bam").write_text("", encoding="utf-8")
    (indir / "notes.bam").write_text("", encoding="utf-8")
    (indir / "ipyrad_map_stats_run.json").write_text("{}", encoding="utf-8")
    (indir / "ipyrad_map_stats_run.txt").write_text("old", encoding="utf-8")
    (indir / "ipyrad_map_stats_orphan.txt").write_text("orphan", encoding="utf-8")

    def _fake_read_bam_sample_name(path: Path) -> str:
        if path.name.startswith("alpha"):
            return "alpha.trimmed"
        return "beta"

    monkeypatch.setattr(
        "ipyrad2.mapper.legacy_map_bams.read_bam_sample_name",
        _fake_read_bam_sample_name,
    )

    plan = plan_legacy_map_migration(indir, outdir)

    assert [item.source_bam.name for item in plan.bam_migrations] == [
        "alpha.trimmed.filtered.bam",
        "beta.trimmed.sorted.bam",
    ]
    assert [item.output_bam.name for item in plan.bam_migrations] == [
        "alpha.trimmed.sorted.bam",
        "beta.trimmed.sorted.bam",
    ]
    assert [item.canonical_sample for item in plan.bam_migrations] == ["alpha", "beta"]
    assert [path.name for path in plan.skipped_bams] == ["notes.bam"]
    assert [item.output_json.name for item in plan.stats_migrations] == [
        "ipyrad_map_stats_run.json",
    ]
    assert [item.output_txt.name for item in plan.stats_migrations] == [
        "ipyrad_map_stats_run.txt",
    ]
    assert [path.name for path in plan.orphan_stats_txt] == ["ipyrad_map_stats_orphan.txt"]


def test_plan_legacy_map_migration_accepts_plain_filtered_bams_without_stats(
    monkeypatch,
    tmp_path: Path,
) -> None:
    indir = tmp_path / "indir"
    outdir = tmp_path / "outdir"
    indir.mkdir()
    (indir / f"alpha{LEGACY_PLAIN_BAM_SUFFIX}").write_text("", encoding="utf-8")
    (indir / "notes.bam").write_text("", encoding="utf-8")

    monkeypatch.setattr(
        "ipyrad2.mapper.legacy_map_bams.read_bam_sample_name",
        lambda _path: "alpha.trimmed",
    )

    plan = plan_legacy_map_migration(indir, outdir)

    assert [item.source_bam.name for item in plan.bam_migrations] == [
        "alpha.filtered.bam",
    ]
    assert [item.output_bam.name for item in plan.bam_migrations] == [
        "alpha.trimmed.sorted.bam",
    ]
    assert [item.canonical_sample for item in plan.bam_migrations] == ["alpha"]
    assert plan.stats_migrations == []
    assert plan.orphan_stats_txt == []
    assert [path.name for path in plan.skipped_bams] == ["notes.bam"]


def test_plan_legacy_map_migration_errors_on_canonical_collision(
    monkeypatch,
    tmp_path: Path,
) -> None:
    indir = tmp_path / "indir"
    outdir = tmp_path / "outdir"
    indir.mkdir()
    (indir / "one.trimmed.filtered.bam").write_text("", encoding="utf-8")
    (indir / "two.trimmed.filtered.bam").write_text("", encoding="utf-8")

    sample_names = {
        "one.trimmed.filtered.bam": "dup.trimmed",
        "two.trimmed.filtered.bam": "dup",
    }

    monkeypatch.setattr(
        "ipyrad2.mapper.legacy_map_bams.read_bam_sample_name",
        lambda path: sample_names[path.name],
    )

    with pytest.raises(IPyradError, match="collide after stripping"):
        plan_legacy_map_migration(indir, outdir)


def test_rewrite_bam_rg_tags_runs_addreplacerg_reheader_and_index(
    monkeypatch,
    tmp_path: Path,
) -> None:
    migration = LegacyBamMigration(
        source_bam=tmp_path / "sample.trimmed.filtered.bam",
        output_bam=tmp_path / "out" / f"sample{CURRENT_BAM_SUFFIX}",
        output_index=tmp_path / "out" / f"sample{CURRENT_BAM_SUFFIX}.csi",
        legacy_sample="sample.trimmed",
        canonical_sample="sample",
    )
    calls: list[tuple[list[list[str]], Path | None]] = []
    captured_header_text: list[str] = []

    def _fake_run_pipeline(cmds, outfile=None, **_kwargs):
        calls.append((cmds, outfile))
        cmd = cmds[0]
        if cmd[1:3] == ["view", "-H"]:
            return (
                0,
                b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:sample.trimmed\tSM:sample.trimmed\n",
                b"",
            )
        if cmd[1] == "reheader":
            captured_header_text.append(Path(cmd[3]).read_text(encoding="utf-8"))
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.mapper.legacy_map_bams.run_pipeline", _fake_run_pipeline)

    _rewrite_bam_rg_tags(migration, threads=3)

    add_rg_cmd = calls[0][0][0]
    assert add_rg_cmd[:7] == [
        add_rg_cmd[0],
        "addreplacerg",
        "-m",
        "overwrite_all",
        "-w",
        "-r",
        "@RG\tID:sample\tSM:sample",
    ]
    assert add_rg_cmd[-1] == str(migration.source_bam)
    assert calls[1][0][0][1:3] == ["view", "-H"]
    assert calls[2][0][0][1] == "reheader"
    assert calls[2][1] == migration.output_bam
    assert calls[3][0][0][1] == "index"
    assert "@RG\tID:sample\tSM:sample\n" in captured_header_text[0]
    assert "sample.trimmed" not in captured_header_text[0]


def test_rewrite_stats_sidecar_normalizes_samples_and_regenerates_text(tmp_path: Path) -> None:
    source_json = tmp_path / "ipyrad_map_stats_run.json"
    output_json = tmp_path / "out" / "ipyrad_map_stats_run.json"
    output_txt = tmp_path / "out" / "ipyrad_map_stats_run.txt"
    payload = {
        "command": "ipyrad2 map -d data -r ref.fa -o out",
        "is_paired": True,
        "applied_mapping_summary": [
            {
                "sample": "sample.trimmed",
                "input_templates": 10,
                "reads_removed_unmapped_or_nonprimary": 1,
                "reads_removed_same_scaffold_pairing": 0,
                "duplicate_records_removed": 0,
                "templates_in_final_bam": 9,
                "fraction_input_templates_retained_in_final_bam": 0.9,
            }
        ],
        "assemble_read_filter_preview": {
            "description": "These preview thresholds were not applied during mapping.",
            "flags": "-qm/--min-map-q, -ms/--max-softclip, -me/--max-nm, -mt/--max-tlen",
            "mode_note": "# Preview mode: pair-level thresholds evaluated on final BAM templates.",
            "mapq_threshold": 20,
            "soft_clipped_bases_threshold": 25,
            "nm_threshold": 50,
            "absolute_tlen_threshold": 2000,
            "filter_effects": [
                {
                    "sample": "sample.trimmed",
                    "templates_failing_min_mapq_20": 1,
                    "templates_failing_max_softclip_25": 0,
                    "templates_failing_max_nm_50": 0,
                    "templates_failing_max_abs_tlen_2000": 0,
                    "templates_passing_all_preview_filters": 8,
                    "fraction_templates_passing_all_preview_filters": 0.8,
                }
            ],
            "metric_summaries": [
                {
                    "sample": "sample.trimmed",
                    "min_mapq_mean": 30.0,
                    "min_mapq_median": 30.0,
                    "min_mapq_stdev": 0.0,
                    "max_softclip_mean": 0.0,
                    "max_softclip_median": 0.0,
                    "max_softclip_stdev": 0.0,
                    "max_nm_mean": 1.0,
                    "max_nm_median": 1.0,
                    "max_nm_stdev": 0.0,
                    "abs_tlen_mean": 100.0,
                    "abs_tlen_median": 100.0,
                    "abs_tlen_stdev": 0.0,
                }
            ],
        },
    }
    source_json.write_text(json.dumps(payload), encoding="utf-8")

    _rewrite_stats_sidecar(
        MapStatsMigration(
            source_json=source_json,
            output_json=output_json,
            output_txt=output_txt,
        )
    )

    migrated = json.loads(output_json.read_text(encoding="utf-8"))
    assert migrated["command"] == payload["command"]
    assert migrated["applied_mapping_summary"][0]["sample"] == "sample"
    assert migrated["assemble_read_filter_preview"]["filter_effects"][0]["sample"] == "sample"
    assert migrated["assemble_read_filter_preview"]["metric_summaries"][0]["sample"] == "sample"
    text = output_txt.read_text(encoding="utf-8")
    assert text.startswith(f"CMD: {payload['command']}\n\n")
    assert "sample.trimmed" not in text
    assert "sample" in text


def test_migrate_legacy_map_outputs_dry_run_does_not_write(monkeypatch, tmp_path: Path) -> None:
    indir = tmp_path / "indir"
    outdir = tmp_path / "outdir"
    indir.mkdir()
    (indir / "sample.trimmed.filtered.bam").write_text("", encoding="utf-8")

    monkeypatch.setattr(
        "ipyrad2.mapper.legacy_map_bams.read_bam_sample_name",
        lambda _path: "sample.trimmed",
    )

    plan = migrate_legacy_map_outputs(
        indir=indir,
        outdir=outdir,
        dry_run=True,
    )

    assert len(plan.bam_migrations) == 1
    assert not outdir.exists()


def test_migrate_legacy_map_outputs_dry_run_accepts_plain_filtered_bams_without_stats(
    monkeypatch,
    tmp_path: Path,
) -> None:
    indir = tmp_path / "indir"
    outdir = tmp_path / "outdir"
    indir.mkdir()
    (indir / "sample.filtered.bam").write_text("", encoding="utf-8")

    monkeypatch.setattr(
        "ipyrad2.mapper.legacy_map_bams.read_bam_sample_name",
        lambda _path: "sample.trimmed",
    )

    plan = migrate_legacy_map_outputs(
        indir=indir,
        outdir=outdir,
        dry_run=True,
    )

    assert [item.output_bam.name for item in plan.bam_migrations] == [
        "sample.trimmed.sorted.bam",
    ]
    assert plan.stats_migrations == []


def test_migration_script_main_parses_dry_run(monkeypatch, tmp_path: Path) -> None:
    module = _load_script_module()
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        module,
        "set_log_level",
        lambda level: captured.setdefault("log_level", level),
    )

    def _fake_migrate_legacy_map_outputs(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(module, "migrate_legacy_map_outputs", _fake_migrate_legacy_map_outputs)
    monkeypatch.setattr(
        "sys.argv",
        [
            "migrate_legacy_map_bams.py",
            "--indir",
            str(tmp_path / "indir"),
            "--outdir",
            str(tmp_path / "outdir"),
            "--dry-run",
        ],
    )

    assert module.main() == 0
    assert captured["log_level"] == "INFO"
    assert captured["dry_run"] is True
    assert captured["threads"] == 1
