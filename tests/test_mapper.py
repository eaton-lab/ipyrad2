from __future__ import annotations

import json
from pathlib import Path

from ipyrad2.mapper.map_samples_prep import MapperSamplePlan
from ipyrad2.mapper.map_stats import MappingJobResult
from ipyrad2.mapper import mapper as mapper_module
from ipyrad2.mapper.mapper import run_mapper


def _sample_plan(
    name: str,
    source_fastqs: tuple[tuple[Path, Path | None], ...],
) -> MapperSamplePlan:
    return MapperSamplePlan(
        output_name=name,
        source_names=tuple(f"source_{idx}" for idx in range(len(source_fastqs))),
        source_fastqs=source_fastqs,
        is_paired_input=source_fastqs[0][1] is not None,
    )


def test_run_mapper_uses_paired_stats_and_writes_report(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "mapped"

    call_log = []

    def _fake_prepare_map_samples(*, fastqs, delim_str, delim_idx, imap, tmpdir, unmate):
        del fastqs, delim_str, delim_idx, imap, tmpdir, unmate
        return {
            "sample": _sample_plan(
                "sample",
                ((Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz")),),
            )
        }, True

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level, max_workers
        call_log.append((msg, jobs))
        if msg == "Mapping":
            return {
                "sample": MappingJobResult(
                    sname="sample",
                    bam_path=outdir / "sample.trimmed.sorted.bam",
                    is_paired=True,
                    nreads_processed=10,
                    nreads_filtered_before_bam_by_unmapped_or_nonprimary=2,
                    nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold=2,
                    nreads_written_before_duplicate_removal=8,
                    duplicate_stats={"duplicate_total": 1},
                )
            }
        return {
            "sample": {
                "input_templates": 5,
                "reads_removed_unmapped_or_nonprimary": 2,
                "reads_removed_same_scaffold_pairing": 1,
                "duplicate_records_removed": 1,
                "templates_in_final_bam": 4,
                "fraction_input_templates_retained_in_final_bam": 0.8,
                "templates_failing_min_mapq_20": 1,
                "templates_failing_max_softclip_25": 0,
                "templates_failing_max_nm_50": 0,
                "templates_failing_max_abs_tlen_2000": 0,
                "templates_passing_all_preview_filters": 4,
                "fraction_templates_passing_all_preview_filters": 1.0,
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
        }

    monkeypatch.setattr("ipyrad2.mapper.mapper.prepare_map_samples", _fake_prepare_map_samples)
    monkeypatch.setattr("ipyrad2.mapper.mapper._check_mapper_dependencies", lambda: None)
    index_calls = []

    def _fake_index(_reference, force_reindex=False):
        index_calls.append((_reference, force_reindex))

    monkeypatch.setattr("ipyrad2.mapper.mapper._index_ref_with_bwa", _fake_index)
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_with_pool", _fake_run_with_pool)

    run_mapper(
        fastqs=[tmp_path / "ignored.fastq.gz"],
        outdir=outdir,
        reference=reference,
        imap=None,
        unmate=False,
        cores=6,
        threads=3,
        force=False,
        reindex_reference=False,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
        delim_str=None,
        delim_idx=1,
        log_level="WARNING",
        logged_command="ipyrad2 map -d ignored.fastq.gz -r ref.fa -o mapped",
    )

    assert [msg for msg, _jobs in call_log] == ["Mapping", "Gathering mapping stats"]
    assert index_calls == [(reference, False)]
    mapping_jobs = call_log[0][1]
    stats_jobs = call_log[1][1]
    assert list(mapping_jobs) == ["sample"]
    assert mapping_jobs["sample"][0].__name__ == "_map_sample"
    assert stats_jobs["sample"][0].__name__ == "collect_paired_bam_stats"

    stats_files = sorted(outdir.glob("ipyrad_map_stats_*.txt"))
    assert len(stats_files) == 1
    report = stats_files[0].read_text(encoding="utf-8")
    report_json = json.loads((outdir / "ipyrad_map_stats_0.json").read_text(encoding="utf-8"))
    assert report.startswith("CMD: ipyrad2 map -d ignored.fastq.gz -r ref.fa -o mapped\n\n")
    assert "# Final BAMs are coordinate sorted and indexed." in report
    assert "# Paired-end final BAMs keep only mapped mates on the same scaffold." in report
    assert "## Applied mapping summary" in report
    assert "## Assemble read-filter preview (not applied during mapping)" in report
    assert "# MAPQ threshold: 20" in report
    assert "sample" in report
    assert report_json["command"] == "ipyrad2 map -d ignored.fastq.gz -r ref.fa -o mapped"
    assert report_json["applied_mapping_summary"][0]["sample"] == "sample"
    assert report_json["assemble_read_filter_preview"]["filter_effects"][0]["sample"] == "sample"


def test_run_mapper_unmate_uses_single_end_stats_and_threads_flag(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "mapped"

    call_log = []
    prepare_calls = []

    def _fake_prepare_map_samples(*, fastqs, delim_str, delim_idx, imap, tmpdir, unmate):
        del fastqs, delim_str, delim_idx, imap, tmpdir
        prepare_calls.append(unmate)
        return {
            "sample": _sample_plan(
                "sample",
                ((Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz")),),
            )
        }, True

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level, max_workers
        call_log.append((msg, jobs))
        if msg == "Mapping":
            return {
                "sample": MappingJobResult(
                    sname="sample",
                    bam_path=outdir / "sample.trimmed.sorted.bam",
                    is_paired=False,
                    nreads_processed=10,
                    nreads_filtered_before_bam_by_unmapped_or_nonprimary=2,
                    nreads_written_before_duplicate_removal=8,
                    duplicate_stats={},
                )
            }
        return {
            "sample": {
                "input_reads": 10,
                "reads_removed_unmapped_or_nonprimary": 2,
                "reads_in_final_bam": 8,
                "fraction_input_reads_retained_in_final_bam": 0.8,
                "reads_failing_min_mapq_20": 1,
                "reads_failing_max_softclip_25": 0,
                "reads_failing_max_nm_50": 0,
                "reads_passing_all_preview_filters": 8,
                "fraction_reads_passing_all_preview_filters": 1.0,
                "mapq_mean": 30.0,
                "mapq_median": 30.0,
                "mapq_stdev": 0.0,
                "softclip_mean": 0.0,
                "softclip_median": 0.0,
                "softclip_stdev": 0.0,
                "nm_mean": 1.0,
                "nm_median": 1.0,
                "nm_stdev": 0.0,
            }
        }

    monkeypatch.setattr("ipyrad2.mapper.mapper.prepare_map_samples", _fake_prepare_map_samples)
    monkeypatch.setattr("ipyrad2.mapper.mapper._check_mapper_dependencies", lambda: None)
    index_calls = []

    def _fake_index(_reference, force_reindex=False):
        index_calls.append((_reference, force_reindex))

    monkeypatch.setattr("ipyrad2.mapper.mapper._index_ref_with_bwa", _fake_index)
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_with_pool", _fake_run_with_pool)

    run_mapper(
        fastqs=[tmp_path / "ignored_R1.fastq.gz", tmp_path / "ignored_R2.fastq.gz"],
        outdir=outdir,
        reference=reference,
        imap=None,
        unmate=True,
        cores=6,
        threads=3,
        force=False,
        reindex_reference=False,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
        delim_str=None,
        delim_idx=1,
        log_level="WARNING",
    )

    assert prepare_calls == [True]
    assert [msg for msg, _jobs in call_log] == ["Mapping", "Gathering mapping stats"]
    assert index_calls == [(reference, False)]
    mapping_jobs = call_log[0][1]
    stats_jobs = call_log[1][1]
    assert mapping_jobs["sample"][1]["is_paired"] is False
    assert mapping_jobs["sample"][1]["unmate"] is True
    assert mapping_jobs["sample"][1]["sample_plan"].source_fastqs == (
        (Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz")),
    )
    assert stats_jobs["sample"][0].__name__ == "collect_single_end_bam_stats"

    stats_files = sorted(outdir.glob("ipyrad_map_stats_*.txt"))
    assert len(stats_files) == 1
    report = stats_files[0].read_text(encoding="utf-8")
    report_json = json.loads((outdir / "ipyrad_map_stats_0.json").read_text(encoding="utf-8"))
    assert "# Paired-end final BAMs keep only mapped mates on the same scaffold." not in report
    assert "input_reads" in report
    assert "command" not in report_json


def test_map_sample_unmate_cleans_temp_fastq_after_success(monkeypatch, tmp_path: Path) -> None:
    outdir = tmp_path / "mapped"
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    sample_r1 = tmp_path / "sample_R1.fastq.gz"
    sample_r2 = tmp_path / "sample_R2.fastq.gz"
    sample_r1.write_bytes(b"r1")
    sample_r2.write_bytes(b"r2")
    unmated = mapper_module._unmated_fastq_path("sample", outdir)
    pipeline_calls = []

    def _fake_run_pipeline(cmds):
        pipeline_calls.append(cmds)
        assert unmated.exists()
        assert cmds[0][-1] == str(unmated)
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", _fake_run_pipeline)
    monkeypatch.setattr(
        "ipyrad2.mapper.mapper._load_save_counts",
        lambda _path: {
            "records_processed": 10,
            "records_filter_rejected": 2,
            "records_filter_accepted": 8,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.mapper.mapper._finalize_indexed_bam",
        lambda _tmp_bam, final_bam: final_bam.write_bytes(b"bam"),
    )

    result = mapper_module._map_sample(
        sname="sample",
        sample_plan=_sample_plan("sample", ((sample_r1, sample_r2),)),
        reference=reference,
        outdir=outdir,
        threads=2,
        is_paired=False,
        unmate=True,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
    )

    assert result.is_paired is False
    assert pipeline_calls
    assert not unmated.exists()


def test_map_sample_unmate_cleans_temp_fastq_after_failure(monkeypatch, tmp_path: Path) -> None:
    outdir = tmp_path / "mapped"
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    sample_r1 = tmp_path / "sample_R1.fastq.gz"
    sample_r2 = tmp_path / "sample_R2.fastq.gz"
    sample_r1.write_bytes(b"r1")
    sample_r2.write_bytes(b"r2")
    unmated = mapper_module._unmated_fastq_path("sample", outdir)

    def _fake_run_pipeline(cmds):
        assert unmated.exists()
        assert cmds[0][-1] == str(unmated)
        raise RuntimeError("mapping failed")

    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", _fake_run_pipeline)

    try:
        mapper_module._map_sample(
            sname="sample",
            sample_plan=_sample_plan("sample", ((sample_r1, sample_r2),)),
            reference=reference,
            outdir=outdir,
            threads=2,
            is_paired=False,
            unmate=True,
            mark_dups_by_coords=False,
            mark_dups_by_umis=False,
        )
    except RuntimeError as exc:
        assert str(exc) == "mapping failed"
    else:
        raise AssertionError("expected RuntimeError")

    assert not unmated.exists()


def test_cleanup_stale_materialized_fastqs_removes_selected_samples_only(tmp_path: Path) -> None:
    outdir = tmp_path / "mapped"
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(parents=True)
    stale_keep = mapper_module._unmated_fastq_path("keep", outdir)
    stale_drop = mapper_module._unmated_fastq_path("drop", outdir)
    stale_keep.write_bytes(b"keep")
    stale_drop.write_bytes(b"drop")
    stale_r1, stale_r2 = mapper_module._merged_fastq_paths("drop", outdir)
    stale_r1.write_bytes(b"r1")
    stale_r2.write_bytes(b"r2")

    mapper_module._cleanup_stale_materialized_fastqs(["drop"], outdir)

    assert stale_keep.exists()
    assert not stale_drop.exists()
    assert not stale_r1.exists()
    assert not stale_r2.exists()


def test_map_sample_imap_merge_cleans_temp_fastqs_after_success(monkeypatch, tmp_path: Path) -> None:
    outdir = tmp_path / "mapped"
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    sample_a_r1 = tmp_path / "sampleA_R1.fastq.gz"
    sample_a_r2 = tmp_path / "sampleA_R2.fastq.gz"
    sample_b_r1 = tmp_path / "sampleB_R1.fastq.gz"
    sample_b_r2 = tmp_path / "sampleB_R2.fastq.gz"
    for path in (sample_a_r1, sample_a_r2, sample_b_r1, sample_b_r2):
        path.write_bytes(path.name.encode("utf-8"))
    merged_r1, merged_r2 = mapper_module._merged_fastq_paths("merged", outdir)

    def _fake_run_pipeline(cmds):
        assert merged_r1.exists()
        assert merged_r2.exists()
        assert cmds[0][-2:] == [str(merged_r1), str(merged_r2)]
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", _fake_run_pipeline)
    monkeypatch.setattr(
        "ipyrad2.mapper.mapper._load_save_counts",
        lambda _path: {
            "records_processed": 10,
            "records_filter_rejected": 2,
            "records_filter_accepted": 8,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.mapper.mapper._finalize_indexed_bam",
        lambda _tmp_bam, final_bam: final_bam.write_bytes(b"bam"),
    )

    result = mapper_module._map_sample(
        sname="merged",
        sample_plan=_sample_plan(
            "merged",
            ((sample_a_r1, sample_a_r2), (sample_b_r1, sample_b_r2)),
        ),
        reference=reference,
        outdir=outdir,
        threads=2,
        is_paired=True,
        unmate=False,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
    )

    assert result.is_paired is True
    assert not merged_r1.exists()
    assert not merged_r2.exists()


def test_run_mapper_imap_skip_does_not_materialize_merged_fastqs(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "mapped"
    outdir.mkdir(parents=True)
    (outdir / "merged.trimmed.sorted.bam").write_bytes(b"bam")
    (outdir / "merged.trimmed.sorted.bam.csi").write_bytes(b"csi")

    def _fake_prepare_map_samples(*, fastqs, delim_str, delim_idx, imap, tmpdir, unmate):
        del fastqs, delim_str, delim_idx, imap, tmpdir, unmate
        return {
            "merged": _sample_plan(
                "merged",
                (
                    (Path("sampleA_R1.fastq.gz"), Path("sampleA_R2.fastq.gz")),
                    (Path("sampleB_R1.fastq.gz"), Path("sampleB_R2.fastq.gz")),
                ),
            )
        }, True

    monkeypatch.setattr("ipyrad2.mapper.mapper.prepare_map_samples", _fake_prepare_map_samples)
    monkeypatch.setattr("ipyrad2.mapper.mapper._check_mapper_dependencies", lambda: None)

    try:
        run_mapper(
            fastqs=[tmp_path / "ignored.fastq.gz"],
            outdir=outdir,
            reference=reference,
            imap=tmp_path / "imap.tsv",
            unmate=False,
            cores=6,
            threads=3,
            force=False,
            reindex_reference=False,
            mark_dups_by_coords=False,
            mark_dups_by_umis=False,
            delim_str=None,
            delim_idx=1,
            log_level="WARNING",
        )
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("expected SystemExit")

    merged_r1, merged_r2 = mapper_module._merged_fastq_paths("merged", outdir)
    assert not merged_r1.exists()
    assert not merged_r2.exists()


def test_index_ref_with_bwa_reuses_complete_existing_index(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    for path in mapper_module._reference_bwa_index_paths(reference):
        path.write_text("", encoding="utf-8")

    run_calls = []
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", lambda cmds: run_calls.append(cmds))

    mapper_module._index_ref_with_bwa(reference)

    assert run_calls == []


def test_index_ref_with_bwa_indexes_when_any_sidecar_is_missing(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    index_paths = mapper_module._reference_bwa_index_paths(reference)
    for path in index_paths[:-1]:
        path.write_text("", encoding="utf-8")

    run_calls = []
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", lambda cmds: run_calls.append(cmds))

    mapper_module._index_ref_with_bwa(reference)

    assert run_calls == [[[mapper_module.BIN_BWA, "index", str(reference)]]]


def test_index_ref_with_bwa_reindexes_when_forced(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    for path in mapper_module._reference_bwa_index_paths(reference):
        path.write_text("", encoding="utf-8")

    run_calls = []
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_pipeline", lambda cmds: run_calls.append(cmds))

    mapper_module._index_ref_with_bwa(reference, force_reindex=True)

    assert run_calls == [[[mapper_module.BIN_BWA, "index", str(reference)]]]


def test_run_mapper_force_does_not_trigger_reference_reindex(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "mapped"

    def _fake_prepare_map_samples(*, fastqs, delim_str, delim_idx, imap, tmpdir, unmate):
        del fastqs, delim_str, delim_idx, imap, tmpdir, unmate
        return {
            "sample": _sample_plan(
                "sample",
                ((Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz")),),
            )
        }, True

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level, max_workers
        if msg == "Mapping":
            return {
                "sample": MappingJobResult(
                    sname="sample",
                    bam_path=outdir / "sample.trimmed.sorted.bam",
                    is_paired=True,
                    nreads_processed=10,
                    nreads_filtered_before_bam_by_unmapped_or_nonprimary=2,
                    nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold=2,
                    nreads_written_before_duplicate_removal=8,
                    duplicate_stats={},
                )
            }
        return {
            "sample": {
                "input_templates": 5,
                "reads_removed_unmapped_or_nonprimary": 2,
                "reads_removed_same_scaffold_pairing": 1,
                "duplicate_records_removed": 0,
                "templates_in_final_bam": 4,
                "fraction_input_templates_retained_in_final_bam": 0.8,
                "templates_failing_min_mapq_20": 1,
                "templates_failing_max_softclip_25": 0,
                "templates_failing_max_nm_50": 0,
                "templates_failing_max_abs_tlen_2000": 0,
                "templates_passing_all_preview_filters": 4,
                "fraction_templates_passing_all_preview_filters": 1.0,
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
        }

    index_calls = []

    def _fake_index(_reference, force_reindex=False):
        index_calls.append((_reference, force_reindex))

    monkeypatch.setattr("ipyrad2.mapper.mapper.prepare_map_samples", _fake_prepare_map_samples)
    monkeypatch.setattr("ipyrad2.mapper.mapper._check_mapper_dependencies", lambda: None)
    monkeypatch.setattr("ipyrad2.mapper.mapper._index_ref_with_bwa", _fake_index)
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_with_pool", _fake_run_with_pool)

    run_mapper(
        fastqs=[tmp_path / "ignored.fastq.gz"],
        outdir=outdir,
        reference=reference,
        imap=None,
        unmate=False,
        cores=6,
        threads=3,
        force=True,
        reindex_reference=False,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
        delim_str=None,
        delim_idx=1,
        log_level="WARNING",
    )

    assert index_calls == [(reference, False)]
