from __future__ import annotations

from pathlib import Path

import pytest

from ipyrad2.mapper.map_stats import MappingJobResult
from ipyrad2.mapper.map_stats import collect_paired_bam_stats
from ipyrad2.mapper.map_stats import collect_single_end_bam_stats
from ipyrad2.mapper.map_stats import parse_markdup_report
from ipyrad2.mapper.map_stats import render_map_stats_report


def _sam_line(
    qname: str,
    flag: int,
    rname: str,
    mapq: int,
    cigar: str,
    tlen: int,
    nm: int,
) -> str:
    return (
        f"{qname}\t{flag}\t{rname}\t1\t{mapq}\t{cigar}\t=\t1\t{tlen}\t"
        f"ACGT\t!!!!\tNM:i:{nm}"
    )


def test_collect_single_end_bam_stats(monkeypatch, tmp_path: Path) -> None:
    def _fake_stream_pipeline_lines(_cmds):
        yield _sam_line("r1", 0, "chr1", 30, "50M", 0, 1)
        yield _sam_line("r2", 0, "chr1", 10, "30M20S", 0, 55)

    monkeypatch.setattr(
        "ipyrad2.mapper.map_stats.stream_pipeline_lines",
        _fake_stream_pipeline_lines,
    )

    stats = collect_single_end_bam_stats(
        MappingJobResult(
            sname="sample",
            bam_path=tmp_path / "sample.trimmed.sorted.bam",
            is_paired=False,
            nreads_processed=4,
            nreads_filtered_before_bam_by_unmapped_or_nonprimary=2,
            nreads_written_before_duplicate_removal=2,
            duplicate_stats={},
        )
    )

    assert stats["reads_in_final_bam"] == 2
    assert stats["reads_failing_min_mapq_20"] == 1
    assert stats["reads_failing_max_softclip_25"] == 0
    assert stats["reads_failing_max_nm_50"] == 1
    assert stats["reads_passing_all_preview_filters"] == 1
    assert stats["mapq_mean"] == 20.0
    assert stats["mapq_median"] == 20.0


def test_collect_paired_bam_stats(monkeypatch, tmp_path: Path) -> None:
    def _fake_stream_pipeline_lines(_cmds):
        yield _sam_line("pair1", 99, "chr1", 30, "50M", 100, 1)
        yield _sam_line("pair1", 147, "chr1", 25, "48M2S", -100, 2)
        yield _sam_line("pair2", 99, "chr1", 15, "20S30M", 0, 55)
        yield _sam_line("pair2", 147, "chr2", 35, "50M", 0, 1)
        yield _sam_line("pair3", 65, "chr1", 40, "50M", 0, 0)

    monkeypatch.setattr(
        "ipyrad2.mapper.map_stats.stream_pipeline_lines",
        _fake_stream_pipeline_lines,
    )

    stats = collect_paired_bam_stats(
        MappingJobResult(
            sname="sample",
            bam_path=tmp_path / "sample.trimmed.sorted.bam",
            is_paired=True,
            nreads_processed=6,
            nreads_filtered_before_bam_by_unmapped_or_nonprimary=0,
            nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold=2,
            nreads_written_before_duplicate_removal=6,
            duplicate_stats={"duplicate_total": 2, "duplicate_pair": 1},
        )
    )

    assert stats["input_templates"] == 3
    assert stats["reads_removed_same_scaffold_pairing"] == 2
    assert stats["templates_in_final_bam"] == 3
    assert stats["duplicate_records_removed"] == 2
    assert stats["templates_failing_min_mapq_20"] == 1
    assert stats["templates_failing_max_softclip_25"] == 0
    assert stats["templates_failing_max_nm_50"] == 1
    assert stats["templates_failing_max_abs_tlen_2000"] == 0
    assert stats["templates_passing_all_preview_filters"] == 1
    assert stats["min_mapq_mean"] == pytest.approx(80.0 / 3.0)


def test_parse_markdup_report_and_render_stats_report(tmp_path: Path) -> None:
    report = tmp_path / "markdup.txt"
    report.write_text(
        "DUPLICATE TOTAL 12\n"
        "DUPLICATE PAIR 5\n"
        "DUPLICATE SINGLE 2\n",
        encoding="utf-8",
    )

    parsed = parse_markdup_report(report)
    assert parsed["duplicate_total"] == 12
    assert parsed["duplicate_pair"] == 5
    assert parsed["duplicate_single"] == 2

    rendered = render_map_stats_report(
        {
            "sample": {
                "input_templates": 10,
                "reads_removed_unmapped_or_nonprimary": 4,
                "reads_removed_same_scaffold_pairing": 2,
                "duplicate_records_removed": 1,
                "templates_in_final_bam": 7,
                "fraction_input_templates_retained_in_final_bam": 0.7,
                "templates_failing_min_mapq_20": 1,
                "templates_failing_max_softclip_25": 2,
                "templates_failing_max_nm_50": 3,
                "templates_failing_max_abs_tlen_2000": 4,
                "templates_passing_all_preview_filters": 5,
                "fraction_templates_passing_all_preview_filters": 0.5,
                "min_mapq_mean": 25.0,
                "min_mapq_median": 20.0,
                "min_mapq_stdev": 5.0,
                "max_softclip_mean": 1.0,
                "max_softclip_median": 1.0,
                "max_softclip_stdev": 0.0,
                "max_nm_mean": 2.0,
                "max_nm_median": 2.0,
                "max_nm_stdev": 0.0,
                "abs_tlen_mean": 100.0,
                "abs_tlen_median": 100.0,
                "abs_tlen_stdev": 0.0,
            }
        },
        is_paired=True,
        logged_command="ipyrad2 map -d a_R1.fastq.gz a_R2.fastq.gz -r ref.fa -o OUT",
    )
    assert rendered.startswith(
        "CMD: ipyrad2 map -d a_R1.fastq.gz a_R2.fastq.gz -r ref.fa -o OUT\n\n"
    )
    assert "## Applied mapping summary" in rendered
    assert "## Assemble read-filter preview (not applied during mapping)" in rendered
    assert "# MAPQ threshold: 20" in rendered
    assert "# These preview thresholds were not applied during mapping." in rendered
    assert "# Paired-end final BAMs keep only mapped mates on the same scaffold." in rendered
    assert "input_templates" in rendered
    assert "reads_removed_same_scaffold_pairing" in rendered
    assert "duplicate_records_removed" in rendered
    assert "templates_with_both_mates" not in rendered
    assert "duplicate_pairs_removed" not in rendered
    assert "sample" in rendered
