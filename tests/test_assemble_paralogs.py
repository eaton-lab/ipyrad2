from __future__ import annotations

from pathlib import Path

import pandas as pd

from ipyrad2.assembler.paralogs import aggregate_across_samples
from ipyrad2.assembler.paralogs import bedtools_coverage_counts
from ipyrad2.assembler.paralogs import get_sample_paralog_tables
from ipyrad2.assembler.paralogs import BIN_BED
from ipyrad2.assembler.paralogs import write_per_sample_final_good


def test_bedtools_coverage_counts_uses_sorted_sweep_with_reference_order(
    monkeypatch,
    tmp_path: Path,
) -> None:
    observed: list[list[list[str]]] = []
    bed = tmp_path / "regions.bed"
    bam = tmp_path / "sample.bam"
    out = tmp_path / "counts.tsv"
    ref_info = tmp_path / "REF_info.txt"

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del kwargs
        observed.append(cmds)
        assert outfile == out
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.paralogs.run_pipeline", _fake_run_pipeline)

    bedtools_coverage_counts(
        bed,
        bam,
        out,
        reference_sort_order=ref_info,
    )

    assert observed == [[
        [
            BIN_BED,
            "coverage",
            "-sorted",
            "-g",
            str(ref_info),
            "-a",
            str(bed),
            "-b",
            str(bam),
            "-counts",
        ]
    ]]


def test_get_sample_paralog_tables_keeps_has_data_for_loci_with_no_snps(
    monkeypatch,
    tmp_path: Path,
) -> None:
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("chr1\t0\t10\nchr1\t20\t30\n", encoding="utf-8")
    bam = tmp_path / "sample.bam"
    bam.write_text("", encoding="utf-8")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")

    monkeypatch.setattr("ipyrad2.assembler.paralogs.call_vcf_from_bam", lambda **kwargs: kwargs["out_vcf_gz"])
    monkeypatch.setattr("ipyrad2.assembler.paralogs.make_indel_mask_bed", lambda vcf_gz, out_bed, pad_bp=10: out_bed)
    monkeypatch.setattr("ipyrad2.assembler.paralogs.extract_snps_table_tsv", lambda vcf_gz, out_tsv: out_tsv)
    monkeypatch.setattr("ipyrad2.assembler.paralogs.mask_snps_table_with_bed", lambda snps_tsv, mask_bed, out_tsv: out_tsv)
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.read_snps_table",
        lambda _path: pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [1],
                "end": [2],
                "DP": [12],
                "GQ": [50],
                "AD": ["8,4"],
                "GT": ["0/1"],
            }
        ),
    )

    def _fake_all_reads_by_region(*, bam, regions_bed, out_tsv, reference_sort_order=None):
        del bam, regions_bed, reference_sort_order
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [0, 20],
                "end": [10, 30],
                "all_reads": [10, 7],
                "rid": ["chr1:0-10", "chr1:20-30"],
            }
        )
        df.to_csv(out_tsv, sep="\t", index=False)
        return df

    monkeypatch.setattr("ipyrad2.assembler.paralogs.all_reads_by_region", _fake_all_reads_by_region)

    _, df_locus = get_sample_paralog_tables(
        bam=bam,
        regions_bed=regions_bed,
        reference_fasta=reference,
        tmpdir=tmp_path,
        prefix="sample",
        min_map_q=10,
        min_base_q=13,
        softclip_len_threshold=None,
        softclip_frac_max=None,
    )

    assert list(df_locus["rid"]) == ["chr1:0-10", "chr1:20-30"]
    assert list(df_locus["has_data"]) == [True, True]
    assert int(df_locus.loc[df_locus["rid"] == "chr1:20-30", "n_snps"].iloc[0]) == 0
    assert bool(df_locus.loc[df_locus["rid"] == "chr1:20-30", "pass"].iloc[0]) is True


def test_get_sample_paralog_tables_handles_empty_masked_snp_table(
    monkeypatch,
    tmp_path: Path,
) -> None:
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("chr1\t0\t10\nchr1\t20\t30\n", encoding="utf-8")
    bam = tmp_path / "sample.bam"
    bam.write_text("", encoding="utf-8")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")

    monkeypatch.setattr("ipyrad2.assembler.paralogs.call_vcf_from_bam", lambda **kwargs: kwargs["out_vcf_gz"])
    monkeypatch.setattr("ipyrad2.assembler.paralogs.make_indel_mask_bed", lambda vcf_gz, out_bed, pad_bp=10: out_bed)
    monkeypatch.setattr("ipyrad2.assembler.paralogs.extract_snps_table_tsv", lambda vcf_gz, out_tsv: out_tsv)
    monkeypatch.setattr("ipyrad2.assembler.paralogs.mask_snps_table_with_bed", lambda snps_tsv, mask_bed, out_tsv: out_tsv)
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.read_snps_table",
        lambda _path: pd.DataFrame(
            {
                "chrom": pd.Series(dtype="string"),
                "start": pd.Series(dtype="int64"),
                "end": pd.Series(dtype="int64"),
                "DP": pd.Series(dtype="int64"),
                "GQ": pd.Series(dtype="int64"),
                "AD": pd.Series(dtype="string"),
                "GT": pd.Series(dtype="string"),
            }
        ),
    )

    def _fake_all_reads_by_region(*, bam, regions_bed, out_tsv, reference_sort_order=None):
        del bam, regions_bed, reference_sort_order
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [0, 20],
                "end": [10, 30],
                "all_reads": [10, 0],
                "rid": ["chr1:0-10", "chr1:20-30"],
            }
        )
        df.to_csv(out_tsv, sep="\t", index=False)
        return df

    monkeypatch.setattr("ipyrad2.assembler.paralogs.all_reads_by_region", _fake_all_reads_by_region)

    df_sites, df_locus = get_sample_paralog_tables(
        bam=bam,
        regions_bed=regions_bed,
        reference_fasta=reference,
        tmpdir=tmp_path,
        prefix="sample",
        min_map_q=10,
        min_base_q=13,
        softclip_len_threshold=None,
        softclip_frac_max=None,
    )

    assert df_sites.empty
    assert list(df_locus["rid"]) == ["chr1:0-10", "chr1:20-30"]
    assert list(df_locus["has_data"]) == [True, False]
    assert list(df_locus["n_snps"]) == [0, 0]


def test_aggregate_across_samples_counts_failures_only_among_samples_with_data(
    tmp_path: Path,
) -> None:
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text(
        "chr1\t0\t10\n"
        "chr1\t20\t30\n"
        "chr1\t40\t50\n",
        encoding="utf-8",
    )

    def _write_sample(prefix: str, *, has_data: dict[str, bool], good: list[str], bad: list[str]) -> None:
        df = pd.DataFrame({"rid": list(has_data), "has_data": list(has_data.values())})
        df.to_csv(tmp_path / f"{prefix}.locus_metrics.tsv", sep="\t", index=False)

        def _bed_text(rids: list[str]) -> str:
            lines = []
            for rid in rids:
                chrom, span = rid.split(":")
                start, end = span.split("-")
                lines.append(f"{chrom}\t{start}\t{end}\n")
            return "".join(lines)

        (tmp_path / f"{prefix}.good.bed").write_text(_bed_text(good), encoding="utf-8")
        (tmp_path / f"{prefix}.paralog_like.bed").write_text(_bed_text(bad), encoding="utf-8")

    _write_sample(
        "s1",
        has_data={"chr1:0-10": True, "chr1:20-30": True, "chr1:40-50": False},
        good=["chr1:0-10"],
        bad=["chr1:20-30"],
    )
    _write_sample(
        "s2",
        has_data={"chr1:0-10": True, "chr1:20-30": False, "chr1:40-50": True},
        good=["chr1:20-30", "chr1:40-50"],
        bad=["chr1:0-10"],
    )

    metrics = aggregate_across_samples(
        regions_bed=regions_bed,
        sample_prefixes=["s1", "s2"],
        in_dir=tmp_path,
        out_prefix=tmp_path / "paralogs",
        fail_frac_max=0.5,
    )

    metrics = metrics.set_index("rid")
    assert int(metrics.loc["chr1:0-10", "n_data"]) == 2
    assert int(metrics.loc["chr1:0-10", "n_fail"]) == 1
    assert float(metrics.loc["chr1:0-10", "fail_frac_among_data"]) == 0.5
    assert bool(metrics.loc["chr1:0-10", "keep_global"]) is True

    assert int(metrics.loc["chr1:20-30", "n_data"]) == 1
    assert int(metrics.loc["chr1:20-30", "n_fail"]) == 1
    assert bool(metrics.loc["chr1:20-30", "drop_global"]) is True

    assert int(metrics.loc["chr1:40-50", "n_data"]) == 1
    assert int(metrics.loc["chr1:40-50", "n_fail"]) == 0
    assert bool(metrics.loc["chr1:40-50", "keep_global"]) is True


def test_aggregate_across_samples_preserves_denovo_region_order(
    tmp_path: Path,
) -> None:
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text(
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )

    df = pd.DataFrame(
        {
            "rid": ["locus_1:0-10", "locus_2:0-10", "locus_11:0-10"],
            "has_data": [True, True, True],
        }
    )
    df.to_csv(tmp_path / "s1.locus_metrics.tsv", sep="\t", index=False)

    (tmp_path / "s1.good.bed").write_text(
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )
    (tmp_path / "s1.paralog_like.bed").write_text("", encoding="utf-8")

    metrics = aggregate_across_samples(
        regions_bed=regions_bed,
        sample_prefixes=["s1"],
        in_dir=tmp_path,
        out_prefix=tmp_path / "paralogs",
        fail_frac_max=0.5,
    )

    assert list(metrics["rid"]) == ["locus_1:0-10", "locus_2:0-10", "locus_11:0-10"]
    assert (tmp_path / "paralogs.shared_good.final.bed").read_text(encoding="utf-8") == (
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_11\t0\t10\n"
    )


def test_write_per_sample_final_good_accepts_subset_of_denovo_contigs(
    tmp_path: Path,
) -> None:
    phase_dir = tmp_path / "phase"
    bed_dir = tmp_path / "beds"
    phase_dir.mkdir()
    bed_dir.mkdir()
    (tmp_path / "REF_info.txt").write_text(
        "locus_1\t10\n"
        "locus_2\t10\n"
        "locus_3\t10\n"
        "locus_11\t10\n",
        encoding="utf-8",
    )

    (phase_dir / "s1.good.bed").write_text(
        "locus_1\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )
    shared_good_bed = phase_dir / "paralogs.shared_good.final.bed"
    shared_good_bed.write_text(
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )

    write_per_sample_final_good(
        sample_prefixes=["s1"],
        in_dir=phase_dir,
        shared_good_bed=shared_good_bed,
        out_dir=bed_dir,
    )

    assert (bed_dir / "s1.final.good.bed").read_text(encoding="utf-8") == (
        "locus_1\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n"
    )
