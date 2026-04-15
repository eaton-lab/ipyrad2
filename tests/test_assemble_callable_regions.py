from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from ipyrad2.assembler import assemble as assemble_module
from ipyrad2.assembler import beds as beds_module
from ipyrad2.assembler.hdf5_utils import get_fai_values
from ipyrad2.assembler.paralogs import get_sample_paralog_tables
from ipyrad2.assembler.paralogs import read_snps_table
from ipyrad2.utils.exceptions import IPyradError


def test_write_callable_regions_bed_splits_non_acgt_and_preserves_order(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(
        ">locus_2\nACGTRNNAC\n>locus_1\nNNACGT\n",
        encoding="utf-8",
    )
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text(
        "locus_2\t0\t9\nlocus_1\t0\t6\n",
        encoding="utf-8",
    )

    out_bed = beds_module.write_callable_regions_bed(
        regions_bed,
        reference,
        tmp_path / "callable.bed",
    )

    assert out_bed.read_text(encoding="utf-8") == (
        "locus_2\t0\t4\nlocus_2\t7\t9\nlocus_1\t2\t6\n"
    )


def test_write_callable_regions_bed_can_be_empty(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">locus_1\nNNRY-\n", encoding="utf-8")
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("locus_1\t0\t5\n", encoding="utf-8")

    out_bed = beds_module.write_callable_regions_bed(
        regions_bed,
        reference,
        tmp_path / "callable.empty.bed",
    )

    assert out_bed.read_text(encoding="utf-8") == ""


def test_get_reference_sort_order_refreshes_stale_fai(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">locus_1\nACGT\n", encoding="utf-8")
    fai = reference.with_suffix(reference.suffix + ".fai")
    fai.write_text("locus_1\t10\t9\t10\t11\n", encoding="utf-8")
    observed: dict[str, int] = {"faidx_calls": 0}

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del kwargs
        cmd = cmds[0]
        if cmd[:2] == [beds_module.BIN_SAM, "faidx"]:
            observed["faidx_calls"] += 1
            fai.write_text("locus_1\t4\t9\t4\t5\n", encoding="utf-8")
        elif cmd == ["cut", "-f", "1,2", str(fai)]:
            assert outfile is not None
            outfile.write_text("locus_1\t4\n", encoding="utf-8")
        else:
            raise AssertionError(cmd)
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.beds.run_pipeline", _fake_run_pipeline)

    out_path = beds_module.get_reference_sort_order(reference, tmp_path)

    assert observed["faidx_calls"] == 1
    assert fai.read_text(encoding="utf-8") == "locus_1\t4\t9\t4\t5\n"
    assert out_path.read_text(encoding="utf-8") == "locus_1\t4\n"


def test_write_callable_regions_bed_refreshes_reference_fai(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">locus_1\nACGT\n", encoding="utf-8")
    fai = reference.with_suffix(reference.suffix + ".fai")
    fai.write_text("locus_1\t10\t9\t10\t11\n", encoding="utf-8")
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("locus_1\t0\t4\n", encoding="utf-8")
    observed: dict[str, int] = {"faidx_calls": 0}

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del outfile, kwargs
        cmd = cmds[0]
        if cmd[:2] != [beds_module.BIN_SAM, "faidx"]:
            raise AssertionError(cmd)
        observed["faidx_calls"] += 1
        fai.write_text("locus_1\t4\t9\t4\t5\n", encoding="utf-8")
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.beds.run_pipeline", _fake_run_pipeline)

    out_bed = beds_module.write_callable_regions_bed(
        regions_bed,
        reference,
        tmp_path / "callable.refresh.bed",
    )

    assert observed["faidx_calls"] == 1
    assert out_bed.read_text(encoding="utf-8") == "locus_1\t0\t4\n"


def test_get_fai_values_observes_rewritten_fai(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")
    fai = reference.with_suffix(reference.suffix + ".fai")
    fai.write_text("chr1\t10\t6\t10\t11\n", encoding="utf-8")

    assert list(get_fai_values(reference, "length")) == [10]

    fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    stat = fai.stat()
    os.utime(fai, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

    assert list(get_fai_values(reference, "length")) == [4]


def test_read_snps_table_handles_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.snps.tsv"
    path.write_text("", encoding="utf-8")

    df = read_snps_table(path)

    assert df.empty
    assert list(df.columns) == ["chrom", "start", "end", "DP", "GQ", "AD", "GT"]
    assert str(df.dtypes["chrom"]) == "string"
    assert str(df.dtypes["start"]) == "int64"
    assert str(df.dtypes["AD"]) == "string"


def test_normalize_user_loci_bed_rejects_interval_beyond_reference_length(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t4\n", encoding="utf-8")
    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text("chr1\t0\t5\n", encoding="utf-8")

    with pytest.raises(
        IPyradError,
        match="--loci-bed line 1 exceeds reference length for chr1: 5 > 4",
    ):
        assemble_module._normalize_user_loci_bed(loci_bed, tmpdir)


def test_get_sample_paralog_tables_uses_callable_bed_only_for_variant_call(
    monkeypatch,
    tmp_path: Path,
) -> None:
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("chr1\t0\t10\nchr1\t20\t30\n", encoding="utf-8")
    callable_bed = tmp_path / "regions.callable.bed"
    callable_bed.write_text("chr1\t0\t4\nchr1\t6\t10\n", encoding="utf-8")
    bam = tmp_path / "sample.bam"
    bam.write_text("", encoding="utf-8")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGTNNACGT\n", encoding="utf-8")
    calls: dict[str, Path] = {}

    def _fake_call_vcf_from_bam(**kwargs):
        calls["variant_regions_bed"] = kwargs["regions_bed"]
        return kwargs["out_vcf_gz"]

    def _fake_make_indel_mask_bed(vcf_gz, out_bed, pad_bp=10):
        del vcf_gz, pad_bp
        out_bed.write_text("", encoding="utf-8")
        return out_bed

    def _fake_extract_snps_table_tsv(vcf_gz, out_tsv):
        del vcf_gz
        out_tsv.write_text("", encoding="utf-8")
        return out_tsv

    def _fake_mask_snps_table_with_bed(snps_tsv, mask_bed, out_tsv):
        del snps_tsv, mask_bed
        out_tsv.write_text("", encoding="utf-8")
        return out_tsv

    def _fake_all_reads_by_region(*, bam, regions_bed, out_tsv):
        del bam
        calls["coverage_regions_bed"] = regions_bed
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [0, 20],
                "end": [10, 30],
                "all_reads": [8, 5],
                "rid": ["chr1:0-10", "chr1:20-30"],
            }
        )
        df.to_csv(out_tsv, sep="\t", index=False)
        return df

    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.call_vcf_from_bam", _fake_call_vcf_from_bam
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.make_indel_mask_bed", _fake_make_indel_mask_bed
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.extract_snps_table_tsv",
        _fake_extract_snps_table_tsv,
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.mask_snps_table_with_bed",
        _fake_mask_snps_table_with_bed,
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.paralogs.all_reads_by_region", _fake_all_reads_by_region
    )

    df_sites, df_locus = get_sample_paralog_tables(
        bam=bam,
        regions_bed=regions_bed,
        callable_regions_bed=callable_bed,
        reference_fasta=reference,
        tmpdir=tmp_path,
        prefix="sample",
        min_map_q=10,
        min_base_q=13,
        softclip_len_threshold=None,
        softclip_frac_max=None,
    )

    assert df_sites.empty
    assert calls["variant_regions_bed"] == callable_bed
    assert calls["coverage_regions_bed"] == regions_bed
    assert list(df_locus["rid"]) == ["chr1:0-10", "chr1:20-30"]
    assert list(df_locus["has_data"]) == [True, True]


def test_run_paralog_stage_passes_callable_bed_to_sample_jobs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    bed_dir = tmp_path / "beds"
    phase_dir = tmp_path / "phase"
    bed_dir.mkdir()
    phase_dir.mkdir()
    sample_bams = {"sample": tmp_path / "sample.bam"}
    sample_bams["sample"].write_text("", encoding="utf-8")
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("chr1\t0\t10\n", encoding="utf-8")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTNNACGT\n", encoding="utf-8")
    captured: dict[str, object] = {}

    def _fake_write_callable_regions_bed(regions_bed, reference_fasta, out_bed):
        del regions_bed, reference_fasta
        out_bed.write_text("chr1\t0\t4\nchr1\t6\t10\n", encoding="utf-8")
        return out_bed

    def _fake_run_with_pool(jobs, log_level, workers, msg="Processing"):
        del log_level, workers, msg
        captured["job_kwargs"] = jobs["sample"][1]
        return {}

    def _fake_aggregate_across_samples(**kwargs):
        out_bed = Path(f"{kwargs['out_prefix']}.shared_good.final.bed")
        out_bed.write_text("chr1\t0\t10\n", encoding="utf-8")
        return pd.DataFrame({"rid": ["chr1:0-10"], "keep_global": [True]})

    def _fake_sort_bed_by_reference_order(in_bed, out_bed, ref_info):
        del ref_info
        out_bed.write_text(Path(in_bed).read_text(encoding="utf-8"), encoding="utf-8")
        return out_bed

    monkeypatch.setattr(
        assemble_module, "write_callable_regions_bed", _fake_write_callable_regions_bed
    )
    monkeypatch.setattr(assemble_module, "run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr(
        assemble_module, "aggregate_across_samples", _fake_aggregate_across_samples
    )
    monkeypatch.setattr(
        assemble_module,
        "sort_bed_by_reference_order",
        _fake_sort_bed_by_reference_order,
    )
    monkeypatch.setattr(
        assemble_module,
        "write_per_sample_final_good",
        lambda **kwargs: {"sample": kwargs["out_dir"] / "sample.final.good.bed"},
    )

    outputs = assemble_module._run_paralog_stage(
        sample_bams=sample_bams,
        regions_bed=regions_bed,
        reference=reference,
        bed_dir=bed_dir,
        phase_dir=phase_dir,
        min_map_q=20,
        min_base_q=13,
        softclip_len_threshold=20,
        softclip_frac_max=0.5,
        depth_z_max=7.0,
        third_frac_cut=0.1,
        min_3allele_sites=2,
        maf_threshold=0.2,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.5,
        workers=1,
        log_level="WARNING",
    )

    job_kwargs = captured["job_kwargs"]
    assert job_kwargs["regions_bed"] == regions_bed
    assert job_kwargs["callable_regions_bed"] == phase_dir / "loci.callable.paralog.bed"
    assert outputs.shared_loci_bed.read_text(encoding="utf-8") == "chr1\t0\t10\n"


def test_run_variant_stage_chunks_callable_loci_bed(
    monkeypatch,
    tmp_path: Path,
) -> None:
    (tmp_path / "beds").mkdir()
    canonical_loci_bed = tmp_path / "beds" / "loci.bed"
    canonical_loci_bed.write_text("chr1\t0\t10\n", encoding="utf-8")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTNNACGT\n", encoding="utf-8")
    bam = tmp_path / "sample.bam"
    bam.write_text("", encoding="utf-8")
    captured: dict[str, object] = {}
    chunk_bed = tmp_path / "beds" / "chunk-0.bed"

    def _fake_write_callable_regions_bed(regions_bed, reference_fasta, out_bed):
        captured["callable_input"] = regions_bed
        del reference_fasta
        out_bed.write_text("chr1\t0\t4\nchr1\t6\t10\n", encoding="utf-8")
        return out_bed

    def _fake_get_chunked_loci_beds(tmpdir, nchunks, source_bed=None):
        del tmpdir, nchunks
        captured["chunk_source_bed"] = source_bed
        chunk_bed.write_text("chr1\t0\t4\n", encoding="utf-8")
        return [chunk_bed]

    def _fake_run_with_pool(jobs, log_level, workers, msg="Processing"):
        del log_level, workers, msg
        captured["variant_job_kwargs"] = next(iter(jobs.values()))[1]
        return {}

    monkeypatch.setattr(
        assemble_module, "write_callable_regions_bed", _fake_write_callable_regions_bed
    )
    monkeypatch.setattr(
        assemble_module, "get_chunked_loci_beds", _fake_get_chunked_loci_beds
    )
    monkeypatch.setattr(assemble_module, "run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr(
        assemble_module,
        "get_concat_chunk_vcfs",
        lambda tmpdir, threads: tmpdir / "vcfs" / "loci.raw.vcf.gz",
    )
    monkeypatch.setattr(
        assemble_module,
        "get_filtered_vcf",
        lambda tmpdir, min_sample_depth, min_geno_q, min_site_q, threads: (
            tmpdir / "vcfs" / "loci.filtered.vcf.gz"
        ),
    )
    monkeypatch.setattr(
        assemble_module,
        "get_vcf_with_indels_resolved",
        lambda tmpdir, reference, threads: tmpdir / "vcfs" / "variants.resolved.vcf.gz",
    )

    result = assemble_module._run_variant_stage(
        tmpdir=tmp_path,
        reference=reference,
        bam_dict={"sample": bam},
        group_samples_file=None,
        min_map_q=20,
        min_base_q=13,
        min_sample_depth=6,
        min_geno_q=20,
        min_site_q=20,
        cores=1,
        threads=1,
        log_level="WARNING",
        wgs_samples=None,
    )

    assert captured["callable_input"] == canonical_loci_bed
    assert (
        captured["chunk_source_bed"] == tmp_path / "beds" / "loci.callable.variant.bed"
    )
    assert captured["variant_job_kwargs"]["locus_chunk"] == chunk_bed
    assert canonical_loci_bed.read_text(encoding="utf-8") == "chr1\t0\t10\n"
    assert result == tmp_path / "vcfs" / "variants.resolved.vcf.gz"
