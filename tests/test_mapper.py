from __future__ import annotations

from pathlib import Path

from ipyrad2.mapper.map_stats import MappingJobResult
from ipyrad2.mapper.mapper import run_mapper


def test_run_mapper_uses_paired_stats_and_writes_report(monkeypatch, tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "mapped"

    call_log = []

    def _fake_prepare_map_samples(*, fastqs, delim_str, delim_idx, imap, tmpdir):
        del fastqs, delim_str, delim_idx, imap, tmpdir
        return {"sample": (Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz"))}, True

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level, max_workers
        call_log.append((msg, jobs))
        if msg == "Mapping":
            return {
                "sample": MappingJobResult(
                    sname="sample",
                    bam_path=outdir / "sample.filtered.bam",
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
    monkeypatch.setattr("ipyrad2.mapper.mapper._index_ref_with_bwa", lambda _reference: None)
    monkeypatch.setattr("ipyrad2.mapper.mapper.run_with_pool", _fake_run_with_pool)

    run_mapper(
        fastqs=[tmp_path / "ignored.fastq.gz"],
        outdir=outdir,
        reference=reference,
        imap=None,
        cores=6,
        threads=3,
        force=False,
        mark_dups_by_coords=False,
        mark_dups_by_umis=False,
        delim_str=None,
        delim_idx=1,
        log_level="WARNING",
    )

    assert [msg for msg, _jobs in call_log] == ["Mapping", "Gathering mapping stats"]
    mapping_jobs = call_log[0][1]
    stats_jobs = call_log[1][1]
    assert list(mapping_jobs) == ["sample"]
    assert mapping_jobs["sample"][0].__name__ == "_map_sample"
    assert stats_jobs["sample"][0].__name__ == "collect_paired_bam_stats"

    stats_files = sorted(outdir.glob("ipyrad_map_stats_*.txt"))
    assert len(stats_files) == 1
    report = stats_files[0].read_text(encoding="utf-8")
    assert "# Final BAMs are coordinate sorted and indexed." in report
    assert "# Paired-end final BAMs keep only mapped mates on the same scaffold." in report
    assert "## Applied mapping summary" in report
    assert "## Assemble read-filter preview (not applied during mapping)" in report
    assert "# MAPQ threshold: 20" in report
    assert "sample" in report
