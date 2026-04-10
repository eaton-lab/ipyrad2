from __future__ import annotations

import gzip
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from ipyrad2.assembler import run_assembler as exported_run_assembler
from ipyrad2.assembler.assemble import _run_variant_stage
from ipyrad2.assembler.assemble import _run_paralog_stage
from ipyrad2.assembler.assemble import _write_consensus_and_outputs
from ipyrad2.assembler.assemble import _normalize_bam_rename_file
from ipyrad2.assembler.assemble import _normalize_user_loci_bed
from ipyrad2.assembler.assemble import _normalize_populations_file
from ipyrad2.assembler.assemble import run_assembler
from ipyrad2.assembler.hdf5_utils import choose_hdf5_cache_settings
from ipyrad2.assembler.hdf5_utils import choose_unsigned_int_dtype
from ipyrad2.assembler.hdf5_utils import get_fai_values
from ipyrad2.assembler.beds import BIN_BED
from ipyrad2.assembler.beds import get_across_sample_loci_bed
from ipyrad2.assembler.beds import get_coverage_bed_graphs
from ipyrad2.assembler.beds import get_sample_depth_stats_in_final_loci
from ipyrad2.assembler.read_filters import BIN_SAM
from ipyrad2.assembler.read_filters import bam_appears_paired
from ipyrad2.assembler.read_filters import classify_bam_layout
from ipyrad2.assembler.loci import filter_trim_locus
from ipyrad2.assembler.loci import get_consensus
from ipyrad2.assembler.loci import get_indel_overlap_mask_path
from ipyrad2.assembler.loci import get_lowdepth_mask_path
from ipyrad2.assembler.loci import get_paralog_mask_path
from ipyrad2.assembler.loci import make_lowdepth_mask
from ipyrad2.assembler.loci import make_paralog_mask
from ipyrad2.assembler.loci import merge_sample_mask_beds
from ipyrad2.assembler.loci import write_final_outputs
from ipyrad2.assembler.loci import write_assemble_stats_report
from ipyrad2.assembler.loci import write_loci_and_stats_files
from ipyrad2.assembler.read_filters import build_mapped_read_filter_expr
from ipyrad2.assembler.variants import apply_sample_region_masks_to_resolved_vcf
from ipyrad2.assembler.variants import apply_wgs_het_allele_balance_mask
from ipyrad2.assembler.variants import _write_overlapping_indel_cluster_masks
from ipyrad2.assembler.variants import BIN_BCF
from ipyrad2.assembler.variants import compact_resolved_vcf_to_final_loci_contigs
from ipyrad2.assembler.variants import get_group_called_variants_in_vcf_chunks
from ipyrad2.assembler.variants import get_concat_chunk_vcfs
from ipyrad2.assembler.variants import get_indel_overlap_clusters_bed_path
from ipyrad2.assembler.variants import get_vcf_with_indels_resolved
from ipyrad2.assembler.variants import summarize_variant_support_by_sample_type
from ipyrad2.assembler.variants import load_variant_resolution_stats
from ipyrad2.assembler.variants import write_vcf
from ipyrad2.assembler.write_seqs import write_seqs_hdf5
from ipyrad2.assembler.write_snps import write_snps_hdf5
from ipyrad2.utils.parallel import run_pipeline
from ipyrad2.utils.exceptions import IPyradError


def test_assembler_package_exports_active_entrypoint() -> None:
    assert exported_run_assembler is run_assembler


def test_choose_hdf5_cache_settings_is_bounded() -> None:
    small = choose_hdf5_cache_settings(total_ram_bytes=8 * 1024**3)
    medium = choose_hdf5_cache_settings(total_ram_bytes=64 * 1024**3)
    large = choose_hdf5_cache_settings(total_ram_bytes=256 * 1024**3)

    assert small["rdcc_nbytes"] == 128 * 1024**2
    assert small["rdcc_nslots"] == 524_287
    assert medium["rdcc_nbytes"] == 512 * 1024**2
    assert medium["rdcc_nslots"] == 1_000_003
    assert large["rdcc_nbytes"] == 1024**3
    assert large["rdcc_nslots"] == 2_000_003


def test_choose_unsigned_int_dtype_falls_back_to_uint64() -> None:
    assert choose_unsigned_int_dtype(100) == np.dtype(np.uint32)
    assert choose_unsigned_int_dtype(2**32) == np.dtype(np.uint64)


def test_get_fai_values_reads_reference_index_columns(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(
        "chr1\t10\t6\t10\t11\n"
        "chr2\t20\t23\t10\t11\n",
        encoding="utf-8",
    )

    assert list(get_fai_values(reference, "scaffold")) == ["chr1", "chr2"]
    np.testing.assert_array_equal(
        get_fai_values(reference, "length"),
        np.array([10, 20], dtype=np.int64),
    )


def test_build_mapped_read_filter_expr_combines_requested_filters() -> None:
    expr = build_mapped_read_filter_expr(
        is_paired=True,
        max_tlen=2000,
        max_softclip=25,
        max_nm=12,
    )

    assert expr is not None
    assert 'rnext=="=" || rnext==rname' in expr
    assert "tlen>=-2000 && tlen<=2000" in expr
    assert "sclen <= 25" in expr
    assert "[NM] <= 12" in expr


def test_build_mapped_read_filter_expr_ignores_pair_filters_for_single_end() -> None:
    expr = build_mapped_read_filter_expr(
        is_paired=False,
        max_tlen=2000,
        max_softclip=None,
        max_nm=None,
    )

    assert expr is None


def test_classify_bam_layout_uses_primary_mapped_counts(monkeypatch, tmp_path: Path) -> None:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("", encoding="utf-8")
    observed_cmds: list[list[list[str]]] = []

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del outfile, kwargs
        observed_cmds.append(cmds)
        cmd = cmds[0]
        if cmd == [BIN_SAM, "view", "-c", "-F", "0x904", "-f", "0x1", str(bam_file)]:
            return 0, b"12\n", b""
        if cmd == [BIN_SAM, "view", "-c", "-F", "0x905", str(bam_file)]:
            return 0, b"0\n", b""
        raise AssertionError(f"unexpected command: {cmd}")

    monkeypatch.setattr("ipyrad2.assembler.read_filters.run_pipeline", _fake_run_pipeline)

    assert classify_bam_layout(bam_file) == "paired"
    assert observed_cmds == [
        [[BIN_SAM, "view", "-c", "-F", "0x904", "-f", "0x1", str(bam_file)]],
        [[BIN_SAM, "view", "-c", "-F", "0x905", str(bam_file)]],
    ]


def test_classify_bam_layout_rejects_hybrid_primary_mapped_layout(
    monkeypatch,
    tmp_path: Path,
) -> None:
    bam_file = tmp_path / "hybrid.bam"
    bam_file.write_text("", encoding="utf-8")

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del outfile, kwargs
        cmd = cmds[0]
        if cmd == [BIN_SAM, "view", "-c", "-F", "0x904", "-f", "0x1", str(bam_file)]:
            return 0, b"4\n", b""
        if cmd == [BIN_SAM, "view", "-c", "-F", "0x905", str(bam_file)]:
            return 0, b"3\n", b""
        raise AssertionError(f"unexpected command: {cmd}")

    monkeypatch.setattr("ipyrad2.assembler.read_filters.run_pipeline", _fake_run_pipeline)

    with pytest.raises(
        IPyradError,
        match="mixed single-end and paired-end primary mapped reads",
    ):
        classify_bam_layout(bam_file)


def test_bam_appears_paired_reflects_layout_classifier(monkeypatch, tmp_path: Path) -> None:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("", encoding="utf-8")

    monkeypatch.setattr(
        "ipyrad2.assembler.read_filters.classify_bam_layout",
        lambda _bam_file: "paired",
    )

    assert bam_appears_paired(bam_file) is True


@pytest.mark.parametrize(
    ("is_paired", "expected_cmd2"),
    [
        (False, [BIN_BED, "bamtobed", "-i", "-"]),
        (True, [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]),
    ],
)
def test_get_coverage_bed_graphs_uses_layout_specific_bamtobed_command(
    monkeypatch,
    tmp_path: Path,
    is_paired: bool,
    expected_cmd2: list[str],
) -> None:
    bam_file = tmp_path / "sample.bam"
    bam_file.write_text("", encoding="utf-8")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference.with_suffix(reference.suffix + ".fai").write_text(
        "chr1\t4\t6\t4\t5\n",
        encoding="utf-8",
    )
    observed: dict[str, object] = {}

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del kwargs
        observed["cmds"] = cmds
        if outfile is not None:
            outfile.parent.mkdir(parents=True, exist_ok=True)
            outfile.write_text("", encoding="utf-8")
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.beds.run_pipeline", _fake_run_pipeline)

    out_bed = get_coverage_bed_graphs(
        sname="sample",
        bam_file=bam_file,
        is_paired=is_paired,
        reference=reference,
        tmpdir=tmp_path / "TMP",
        min_map_q=10,
        min_sample_depth=1,
        min_merge_distance=50,
        threads=2,
    )

    assert out_bed == tmp_path / "TMP" / "beds" / "sample.fragments.merged.bed"
    assert observed["cmds"][1] == expected_cmd2


def test_get_consensus_uses_shared_reference_fasta(monkeypatch, tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "vcfs").mkdir(parents=True)
    (tmpdir / "beds").mkdir(parents=True)

    reference_fasta = tmpdir / "consensus_seqs" / "assembly_reference_sequence.consensus.fa"
    reference_fasta.parent.mkdir(parents=True)
    reference_fasta.write_text(">chr1:1-4\nAAAA\n", encoding="utf-8")

    vcf_gz = tmpdir / "vcfs" / "variants.resolved.vcf.gz"
    vcf_gz.write_text("", encoding="utf-8")
    mask_bed = tmpdir / "beds" / "s1.mask.bed"
    mask_bed.write_text("", encoding="utf-8")

    observed_cmds: list[list[list[str]]] = []

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del kwargs
        observed_cmds.append(cmds)
        if outfile is not None:
            outfile.parent.mkdir(parents=True, exist_ok=True)
            outfile.write_text(">chr1:1-4\nAAAA\n", encoding="utf-8")
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.loci.run_pipeline", _fake_run_pipeline)

    out_fasta = get_consensus(
        sname="s1",
        reference_fasta=reference_fasta,
        tmpdir=tmpdir,
        keep_insertions=False,
    )

    assert out_fasta == tmpdir / "consensus_seqs" / "s1.consensus.fa"
    assert observed_cmds == [[
        [
            BIN_BCF,
            "consensus",
            "-f",
            str(reference_fasta),
            "-s",
            "s1",
            "-M",
            "N",
            "--mask",
            str(mask_bed),
            "--mask-with",
            "N",
            "--mark-del",
            "-",
            "--mark-ins",
            "+",
            "--regions-overlap",
            "1",
            str(vcf_gz),
        ],
        ["tr", "-d", "'+"],
    ]]


def test_write_consensus_and_outputs_uses_one_pool_with_stage_specific_consensus_workers(
    monkeypatch,
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")

    observed_calls: list[tuple[list[str], int | None, str]] = []

    monkeypatch.setattr("ipyrad2.assembler.assemble.write_sam_faidx", lambda _tmpdir: _tmpdir / "loci.faidx.txt")
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_reference_in_loci_beds",
        lambda _tmpdir, _reference: _tmpdir / "consensus_seqs" / "assembly_reference_sequence.consensus.fa",
    )
    monkeypatch.setattr("ipyrad2.assembler.assemble.build_locus_fasta_database", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_final_outputs",
        lambda **kwargs: {
            "nloci_before_filtering": 5,
            "nloci_after_filtering": 0,
            "nsites_after_filtering": 0,
            "filter_counts": {},
            "site_totals": {},
            "sample_locus_counts": {},
            "samples_per_locus_counts": {},
            "locus_length_counts": {},
            "alignment_nonmissing_sample_bases": 0,
        },
    )

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level
        observed_calls.append((list(jobs), max_workers, msg))
        return {key: Path(f"/tmp/{key}.fa") for key in jobs}

    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_vcf",
        lambda *args, **kwargs: pytest.fail("write_vcf should not run when no loci survive"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_snps_hdf5",
        lambda *args, **kwargs: pytest.fail("write_snps_hdf5 should not run when no loci survive"),
    )

    with pytest.raises(IPyradError, match="No loci passed final trimming/filtering"):
        _write_consensus_and_outputs(
            name="assembly",
            outdir=tmp_path,
            tmpdir=tmpdir,
            snames=["s1", "s2", "s3", "s4", "s5"],
            reference=reference,
            masks=None,
            sample_masks=None,
            shared_loci_after_delimiting=5,
            shared_loci_after_paralog_filtering=5,
            min_locus_sample_coverage=1,
            min_locus_trim_sample_coverage=1,
            min_locus_length=25,
            max_locus_hetero_frequency=0.3,
            max_locus_variant_frequency=1.0,
            max_sample_hetero_frequency=0.10,
            consensus_workers=3,
            workers=1,
            threads=1,
            log_level="WARNING",
        )

    assert observed_calls == [
        (["s1", "s2", "s3", "s4", "s5"], 3, "Extracting consensus sequences"),
    ]


def test_run_variant_stage_caps_inflight_jobs_to_assemble_worker_budget(
    monkeypatch,
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "vcfs").mkdir(parents=True)
    (tmpdir / "beds").mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")
    (tmpdir / "beds" / "loci.bed").write_text("chr1\t0\t4\n", encoding="utf-8")

    observed: dict[str, object] = {}

    def _fake_get_chunked_loci_beds(_tmpdir, nchunks, source_bed=None):
        observed["chunk_count"] = nchunks
        observed["source_bed"] = source_bed
        chunk0 = tmpdir / "beds" / "chunk-0.bed"
        chunk1 = tmpdir / "beds" / "chunk-1.bed"
        chunk0.write_text("chr1\t0\t10\n", encoding="utf-8")
        chunk1.write_text("chr1\t10\t20\n", encoding="utf-8")
        return [chunk0, chunk1]

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level
        observed["msg"] = msg
        observed["max_workers"] = max_workers
        observed["jobs"] = jobs
        return {key: None for key in jobs}

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_chunked_loci_beds", _fake_get_chunked_loci_beds)
    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr("ipyrad2.assembler.assemble.get_concat_chunk_vcfs", lambda *args, **kwargs: None)
    monkeypatch.setattr("ipyrad2.assembler.assemble.get_filtered_vcf", lambda *args, **kwargs: None)
    monkeypatch.setattr("ipyrad2.assembler.assemble.get_vcf_with_indels_resolved", lambda *args, **kwargs: tmpdir / "vcfs" / "variants.resolved.vcf.gz")

    _run_variant_stage(
        tmpdir=tmpdir,
        reference=reference,
        bam_dict={"s1": tmp_path / "s1.bam", "s2": tmp_path / "s2.bam"},
        group_samples_file=None,
        min_map_q=20,
        min_base_q=20,
        min_sample_depth=4,
        min_geno_q=20,
        min_site_q=20,
        cores=6,
        threads=3,
        log_level="WARNING",
    )

    assert observed["msg"] == "Calling variants"
    assert observed["max_workers"] == 2
    assert observed["chunk_count"] == 8
    assert observed["source_bed"] == tmpdir / "beds" / "loci.callable.variant.bed"
    first_job = next(iter(observed["jobs"].values()))
    assert first_job[1]["threads"] == 2


@pytest.mark.parametrize("use_group_samples_file", [False, True])
def test_get_group_called_variants_in_vcf_chunks_sets_bcftools_group_arg(
    monkeypatch,
    tmp_path: Path,
    use_group_samples_file: bool,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "vcfs").mkdir(parents=True)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")
    locus_chunk = tmp_path / "chunk-0.bed"
    locus_chunk.write_text("chr1\t0\t4\n", encoding="utf-8")

    group_samples_file = None
    expected_group_arg = "-"
    if use_group_samples_file:
        group_samples_file = tmp_path / "groups.tsv"
        group_samples_file.write_text("s1\tpop1\ns2\tpop2\n", encoding="utf-8")
        expected_group_arg = str(group_samples_file)

    observed_cmds: list[list[list[str]]] = []

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del outfile, kwargs
        observed_cmds.append(cmds)
        return 0, b"", b""

    monkeypatch.setattr("ipyrad2.assembler.variants.run_pipeline", _fake_run_pipeline)

    out_vcf = get_group_called_variants_in_vcf_chunks(
        tmpdir=tmpdir,
        reference=reference,
        bam_files=[tmp_path / "s1.bam", tmp_path / "s2.bam"],
        locus_chunk=locus_chunk,
        min_map_q=20,
        min_base_q=20,
        threads=2,
        group_samples_file=group_samples_file,
    )

    assert out_vcf == tmpdir / "vcfs" / "chunk-0.vcf.gz"
    assert len(observed_cmds) == 1
    cmd2 = observed_cmds[0][1]
    assert cmd2[cmd2.index("-G") + 1] == expected_group_arg


def test_get_concat_chunk_vcfs_raises_clear_error_when_no_chunks_exist(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "vcfs").mkdir(parents=True)

    with pytest.raises(IPyradError, match="No chunk VCFs found"):
        get_concat_chunk_vcfs(tmpdir, threads=1)


def test_apply_wgs_het_allele_balance_mask_masks_out_of_balance_hets(
    monkeypatch,
    tmp_path: Path,
) -> None:
    vcf_gz = tmp_path / "loci.filtered.vcf.gz"
    with gzip.open(vcf_gz, "wt", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trad\twgs\n")
        out.write("chr1\t1\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:18,2\n")
        out.write("chr1\t2\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:9,11\n")

    def _fake_run_pipeline(cmds, outfile=None, **kwargs):
        del outfile, kwargs
        first = cmds[0]
        if len(cmds) == 2 and first[1] == "+fill-tags":
            plain_vcf = Path(first[2])
            out_gz = Path(cmds[1][cmds[1].index("-o") + 1])
            with plain_vcf.open("r", encoding="utf-8") as src, gzip.open(out_gz, "wt", encoding="utf-8") as dst:
                dst.write(src.read())
            return 0, b"", b""
        if len(cmds) == 1 and first[1] == "index":
            Path(f"{first[-1]}.csi").write_text("", encoding="utf-8")
            return 0, b"", b""
        raise AssertionError(f"unexpected run_pipeline call: {cmds}")

    monkeypatch.setattr("ipyrad2.assembler.variants.run_pipeline", _fake_run_pipeline)

    stats = apply_wgs_het_allele_balance_mask(vcf_gz, ["wgs"], low=0.20, high=0.80)

    assert stats == {
        "wgs_het_genotypes_masked_by_allele_balance": 1,
        "wgs_het_genotypes_examined_for_allele_balance": 2,
    }
    with gzip.open(vcf_gz, "rt", encoding="utf-8") as handle:
        records = [line.rstrip("\n").split("\t") for line in handle if line and not line.startswith("#")]
    assert records[0][9] == "0/1:10,10"
    assert records[0][10] == "./.:18,2"
    assert records[1][10] == "0/1:9,11"


def test_summarize_variant_support_by_sample_type_counts_support_categories(
    tmp_path: Path,
) -> None:
    vcf_gz = tmp_path / "assembly.vcf.gz"
    with gzip.open(vcf_gz, "wt", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\trad\twgs\n")
        out.write("chr1\t1\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\t0/0\n")
        out.write("chr1\t2\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\t0/1\n")
        out.write("chr1\t3\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\t0/1\n")
        out.write("chr1\t4\t.\tA\tC\t50\tPASS\t.\tGT\t./.\t./.\n")

    stats = summarize_variant_support_by_sample_type(vcf_gz, ["rad"], ["wgs"])

    assert stats == {
        "sites_supported_rad_only": 1,
        "sites_supported_wgs_only": 1,
        "sites_supported_both": 1,
        "sites_supported_neither": 1,
    }


def test_write_overlapping_indel_cluster_masks_prunes_clusters_and_writes_sample_masks(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir(parents=True)

    plain_vcf = vcf_dir / "variants.resolved.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=20>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t5\t.\tAAAG\tA\t50\tPASS\t.\tGT\t1/1\t0/0\n"
        "chr1\t7\t.\tAGAA\tA\t50\tPASS\t.\tGT\t1/1\t0/0\n"
        "chr1\t15\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0\n",
            encoding="utf-8",
    )
    resolved_vcf = vcf_dir / "variants.resolved.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(resolved_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(resolved_vcf)]])

    sample_masks = _write_overlapping_indel_cluster_masks(tmpdir)

    with gzip.open(resolved_vcf, "rt", encoding="utf-8") as handle:
        rows = [
            line.rstrip("\n").split("\t")
            for line in handle
            if line and not line.startswith("#")
        ]
    assert [(row[0], row[1], row[3], row[4]) for row in rows] == [
        ("chr1", "15", "C", "T"),
    ]

    overlap_bed = get_indel_overlap_clusters_bed_path(tmpdir)
    assert overlap_bed.read_text(encoding="utf-8") == "chr1\t4\t10\n"

    s1_mask = get_indel_overlap_mask_path("s1", tmpdir)
    s2_mask = get_indel_overlap_mask_path("s2", tmpdir)
    assert sample_masks == {"s1": s1_mask, "s2": s2_mask}
    assert s1_mask.read_text(encoding="utf-8") == "chr1\t4\t10\n"
    assert s2_mask.read_text(encoding="utf-8") == ""

    pre_overlap_vcf = vcf_dir / "variants.resolved.pre_overlap_clusters.vcf.gz"
    assert pre_overlap_vcf.exists()


def test_write_overlapping_indel_cluster_masks_handles_header_only_vcf(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir(parents=True)

    plain_vcf = vcf_dir / "variants.resolved.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=20>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
        encoding="utf-8",
    )
    resolved_vcf = vcf_dir / "variants.resolved.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(resolved_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(resolved_vcf)]])

    sample_masks = _write_overlapping_indel_cluster_masks(tmpdir)

    assert sample_masks == {
        "s1": get_indel_overlap_mask_path("s1", tmpdir),
        "s2": get_indel_overlap_mask_path("s2", tmpdir),
    }
    assert get_indel_overlap_clusters_bed_path(tmpdir).read_text(encoding="utf-8") == ""
    assert get_indel_overlap_mask_path("s1", tmpdir).read_text(encoding="utf-8") == ""
    assert get_indel_overlap_mask_path("s2", tmpdir).read_text(encoding="utf-8") == ""
    assert load_variant_resolution_stats(tmpdir) == {
        "overlapping_indel_clusters_masked": 0,
        "overlapping_indel_records_removed": 0,
        "overlapping_indel_bp_masked": 0,
        "indel_records_inspected": 0,
    }


def test_merge_sample_mask_beds_includes_overlap_cluster_masks(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8")

    (bed_dir / "s1.lowdepth.mask.bed").write_text("chr1\t0\t5\n", encoding="utf-8")
    (bed_dir / "s1.paralog.mask.bed").write_text("chr1\t20\t25\n", encoding="utf-8")
    get_indel_overlap_mask_path("s1", tmpdir).write_text("chr1\t4\t10\n", encoding="utf-8")

    out_bed = merge_sample_mask_beds("s1", tmpdir)

    assert out_bed.read_text(encoding="utf-8") == "chr1\t0\t10\nchr1\t20\t25\n"


def test_merge_sample_mask_beds_handles_denovo_nested_locus_ids(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_3_8\t100\n"
        "locus_3_16\t100\n",
        encoding="utf-8",
    )

    (bed_dir / "s1.lowdepth.mask.bed").write_text(
        "locus_3_8\t0\t5\n"
        "locus_3_16\t20\t25\n",
        encoding="utf-8",
    )
    (bed_dir / "s1.paralog.mask.bed").write_text(
        "locus_3_8\t4\t10\n"
        "locus_3_16\t24\t30\n",
        encoding="utf-8",
    )

    out_bed = merge_sample_mask_beds("s1", tmpdir)

    assert out_bed.read_text(encoding="utf-8") == (
        "locus_3_8\t0\t10\n"
        "locus_3_16\t20\t30\n"
    )


def test_make_lowdepth_mask_accepts_subset_of_denovo_contigs(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_1\t10\n"
        "locus_2\t10\n"
        "locus_3\t10\n"
        "locus_11\t10\n"
        "locus_10010\t10\n",
        encoding="utf-8",
    )
    (bed_dir / "loci.bed").write_text(
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )
    (bed_dir / "s1.fragments.bedgraph").write_text(
        "locus_1\t0\t10\t5\n"
        "locus_3\t0\t10\t5\n"
        "locus_11\t0\t10\t5\n",
        encoding="utf-8",
    )

    out_bed = make_lowdepth_mask("s1", 1, tmpdir)

    assert out_bed == get_lowdepth_mask_path("s1", tmpdir)
    assert out_bed.read_text(encoding="utf-8") == "locus_2\t0\t10\n"


def test_make_lowdepth_mask_handles_denovo_nested_locus_ids(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_3_8\t100\n"
        "locus_3_16\t100\n",
        encoding="utf-8",
    )
    (bed_dir / "loci.bed").write_text(
        "locus_3_8\t0\t10\n"
        "locus_3_16\t20\t30\n",
        encoding="utf-8",
    )
    (bed_dir / "s1.fragments.bedgraph").write_text(
        "locus_3_8\t0\t10\t5\n",
        encoding="utf-8",
    )

    out_bed = make_lowdepth_mask("s1", 1, tmpdir)

    assert out_bed == get_lowdepth_mask_path("s1", tmpdir)
    assert out_bed.read_text(encoding="utf-8") == "locus_3_16\t20\t30\n"


def test_get_across_sample_loci_bed_handles_denovo_nested_locus_ids(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_3_8\t100\n"
        "locus_3_16\t100\n",
        encoding="utf-8",
    )
    for sample in ("s1", "s2"):
        (bed_dir / f"{sample}.fragments.merged.bed").write_text(
            "locus_3_8\t0\t10\n"
            "locus_3_16\t20\t30\n",
            encoding="utf-8",
        )

    out_bed = get_across_sample_loci_bed(
        ["s1", "s2"],
        min_sample_coverage=2,
        min_merge_distance=0,
        min_locus_length=1,
        suffix=".fragments.merged.bed",
        tmpdir=tmpdir,
    )

    assert out_bed == bed_dir / "loci.bed"
    assert out_bed.read_text(encoding="utf-8") == (
        "locus_3_8\t0\t10\t2\n"
        "locus_3_16\t20\t30\t2\n"
    )


def test_make_paralog_mask_accepts_subset_of_denovo_contigs(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_1\t10\n"
        "locus_2\t10\n"
        "locus_3\t10\n"
        "locus_11\t10\n"
        "locus_10010\t10\n",
        encoding="utf-8",
    )
    (bed_dir / "loci.bed").write_text(
        "locus_1\t0\t10\n"
        "locus_2\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )
    (bed_dir / "s1.final.good.bed").write_text(
        "locus_1\t0\t10\n"
        "locus_3\t0\t10\n"
        "locus_11\t0\t10\n",
        encoding="utf-8",
    )

    out_bed = make_paralog_mask("s1", tmpdir)

    assert out_bed == get_paralog_mask_path("s1", tmpdir)
    assert out_bed.read_text(encoding="utf-8") == "locus_2\t0\t10\n"


def test_normalize_user_loci_bed_sorts_by_reference_and_keeps_bed3(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr2\t100\nchr1\t100\n", encoding="utf-8")

    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text(
        "chr1\t5\t10\textra\n"
        "chr2\t0\t5\tignored\n",
        encoding="utf-8",
    )

    out_bed, nloci = _normalize_user_loci_bed(loci_bed, tmpdir)

    assert out_bed == tmpdir / "beds" / "loci.raw.bed"
    assert nloci == 2
    assert out_bed.read_text(encoding="utf-8") == "chr2\t0\t5\nchr1\t5\t10\n"


def test_run_paralog_stage_normalizes_shared_bed_to_reference_order(
    monkeypatch,
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    phase_dir = tmpdir / "paralogs"
    bed_dir.mkdir(parents=True)
    phase_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text(
        "locus_1\t10\n"
        "locus_2\t10\n"
        "locus_11\t10\n",
        encoding="utf-8",
    )
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">locus_1\nACGTACGTAC\n>locus_2\nACGTACGTAC\n>locus_11\nACGTACGTAC\n",
        encoding="utf-8",
    )
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text(
        "locus_1\t0\t10\nlocus_2\t0\t10\nlocus_11\t0\t10\n",
        encoding="utf-8",
    )

    def _fake_run_with_pool(jobs, log_level, workers, msg=None):
        del jobs, log_level, workers, msg
        return {}

    def _fake_aggregate_across_samples(**kwargs):
        del kwargs
        shared_good = phase_dir / "paralogs.shared_good.final.bed"
        shared_good.write_text(
            "locus_1\t0\t10\n"
            "locus_11\t0\t10\n"
            "locus_2\t0\t10\n",
            encoding="utf-8",
        )
        return pd.DataFrame({"keep_global": [True, True, True]})

    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr("ipyrad2.assembler.assemble.aggregate_across_samples", _fake_aggregate_across_samples)
    monkeypatch.setattr("ipyrad2.assembler.assemble.write_per_sample_final_good", lambda **_kwargs: None)

    final_bed = _run_paralog_stage(
        sample_bams={"sample": tmp_path / "sample.bam"},
        regions_bed=regions_bed,
        reference=reference,
        bed_dir=bed_dir,
        phase_dir=phase_dir,
        min_map_q=40,
        min_base_q=30,
        softclip_len_threshold=20,
        softclip_frac_max=0.25,
        depth_z_max=5.0,
        third_frac_cut=0.10,
        min_3allele_sites=2,
        maf_threshold=0.20,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.10,
        workers=1,
        log_level="INFO",
    )

    expected = "locus_1\t0\t10\nlocus_2\t0\t10\nlocus_11\t0\t10\n"
    assert final_bed.read_text(encoding="utf-8") == expected
    assert (bed_dir / "loci.bed").read_text(encoding="utf-8") == expected


def test_run_paralog_stage_uses_rad_aggregate_for_mixed_shared_bed(
    monkeypatch,
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    phase_dir = tmpdir / "paralogs"
    bed_dir.mkdir(parents=True)
    phase_dir.mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8")
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + ("A" * 100) + "\n", encoding="utf-8")
    regions_bed = tmp_path / "regions.bed"
    regions_bed.write_text("chr1\t0\t30\n", encoding="utf-8")

    def _fake_run_with_pool(jobs, log_level, workers, msg=None):
        del jobs, log_level, workers, msg
        return {}

    def _fake_aggregate_across_samples(*, sample_prefixes, out_prefix, **kwargs):
        del kwargs
        names = tuple(sample_prefixes)
        if names == ("rad",):
            Path(f"{out_prefix}.shared_good.final.bed").write_text("chr1\t0\t10\n", encoding="utf-8")
            return pd.DataFrame(
                {
                    "chrom": ["chr1", "chr1"],
                    "start": [0, 20],
                    "end": [10, 30],
                    "rid": ["chr1:0-10", "chr1:20-30"],
                    "n_data": [1, 1],
                    "n_good": [1, 0],
                    "n_fail": [0, 1],
                    "fail_frac_among_data": [0.0, 1.0],
                    "good_frac_among_data": [1.0, 0.0],
                    "drop_global": [False, True],
                    "keep_global": [True, False],
                }
            )
        Path(f"{out_prefix}.shared_good.final.bed").write_text("chr1\t20\t30\n", encoding="utf-8")
        return pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [0, 20],
                "end": [10, 30],
                "rid": ["chr1:0-10", "chr1:20-30"],
                "n_data": [1, 1],
                "n_good": [0, 1],
                "n_fail": [1, 0],
                "fail_frac_among_data": [1.0, 0.0],
                "good_frac_among_data": [0.0, 1.0],
                "drop_global": [True, False],
                "keep_global": [False, True],
            }
        )

    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr("ipyrad2.assembler.assemble.aggregate_across_samples", _fake_aggregate_across_samples)
    monkeypatch.setattr("ipyrad2.assembler.assemble.write_per_sample_final_good", lambda **_kwargs: None)

    final_bed = _run_paralog_stage(
        sample_bams={"rad": tmp_path / "rad.bam", "wgs": tmp_path / "wgs.bam"},
        regions_bed=regions_bed,
        reference=reference,
        bed_dir=bed_dir,
        phase_dir=phase_dir,
        min_map_q=40,
        min_base_q=30,
        softclip_len_threshold=20,
        softclip_frac_max=0.25,
        depth_z_max=5.0,
        third_frac_cut=0.10,
        min_3allele_sites=2,
        maf_threshold=0.20,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.10,
        workers=1,
        log_level="INFO",
        rad_sample_names=["rad"],
        wgs_sample_names=["wgs"],
    )

    assert final_bed.read_text(encoding="utf-8") == "chr1\t0\t10\n"
    assert (phase_dir / "paralogs.mixed_summary.tsv").exists()
    counts = (phase_dir / "paralogs.mixed.counts.tsv").read_text(encoding="utf-8")
    assert "loci_fail_paralog_rad\t1" in counts
    assert "loci_fail_paralog_wgs\t1" in counts
    assert "loci_fail_paralog_both\t0" in counts
    assert "loci_pass_paralog_rad_fail_paralog_wgs\t1" in counts


def test_normalize_user_loci_bed_rejects_unknown_scaffolds(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8")

    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text("chr2\t0\t5\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="scaffold not present in reference: chr2"):
        _normalize_user_loci_bed(loci_bed, tmpdir)


def test_normalize_user_loci_bed_rejects_overlapping_intervals(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8")

    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text("chr1\t0\t10\nchr1\t5\t12\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="overlapping intervals on chr1"):
        _normalize_user_loci_bed(loci_bed, tmpdir)


def test_normalize_user_loci_bed_rejects_empty_files(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    (tmpdir / "beds").mkdir(parents=True)
    (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8")

    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text("", encoding="utf-8")

    with pytest.raises(IPyradError, match="contains no loci"):
        _normalize_user_loci_bed(loci_bed, tmpdir)


def test_normalize_populations_file_writes_sample_ordered_group_table(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    tmpdir.mkdir(parents=True)
    populations = tmp_path / "groups.tsv"
    populations.write_text(
        "s1\tpopA\n"
        "s2\tpopB\n"
        "s3\tpopA\n",
        encoding="utf-8",
    )

    out_path, imap, minmap = _normalize_populations_file(
        populations=populations,
        tmpdir=tmpdir,
        sample_names=["s2", "s1", "s3"],
    )

    assert out_path == tmpdir / "populations.normalized.tsv"
    assert out_path.read_text(encoding="utf-8") == (
        "s2\tpopB\n"
        "s1\tpopA\n"
        "s3\tpopA\n"
    )
    assert imap == {"popA": ["s1", "s3"], "popB": ["s2"]}
    assert minmap is None


def test_normalize_populations_file_accepts_classic_pop_assign_format(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    tmpdir.mkdir(parents=True)
    populations = tmp_path / "groups.txt"
    populations.write_text(
        "s1 pop1\n"
        "s2 pop2\n"
        "# pop1:1 pop2:2\n",
        encoding="utf-8",
    )

    out_path, imap, minmap = _normalize_populations_file(
        populations=populations,
        tmpdir=tmpdir,
        sample_names=["s2", "s1"],
    )

    assert out_path.read_text(encoding="utf-8") == "s2\tpop2\ns1\tpop1\n"
    assert imap == {"pop1": ["s1"], "pop2": ["s2"]}
    assert minmap == {"pop1": 1, "pop2": 2}


def test_normalize_populations_file_rejects_duplicate_sample_assignments(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    tmpdir.mkdir(parents=True)
    populations = tmp_path / "groups.tsv"
    populations.write_text("s1\tpopA\ns1\tpopB\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="assigns sample\\(s\\) multiple times: s1"):
        _normalize_populations_file(
            populations=populations,
            tmpdir=tmpdir,
            sample_names=["s1"],
        )


def test_normalize_populations_file_rejects_missing_and_extra_samples(tmp_path: Path) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    tmpdir.mkdir(parents=True)

    missing = tmp_path / "missing.tsv"
    missing.write_text("s1\tpopA\n", encoding="utf-8")
    with pytest.raises(IPyradError, match="missing assembled sample\\(s\\): s2"):
        _normalize_populations_file(
            populations=missing,
            tmpdir=tmpdir,
            sample_names=["s1", "s2"],
        )

    extra = tmp_path / "extra.tsv"
    extra.write_text("s1\tpopA\ns2\tpopB\n", encoding="utf-8")
    with pytest.raises(IPyradError, match="not present in this assemble run: s2"):
        _normalize_populations_file(
            populations=extra,
            tmpdir=tmpdir,
            sample_names=["s1"],
        )


def test_normalize_bam_rename_file_parses_partial_basename_map(tmp_path: Path) -> None:
    bam1 = tmp_path / "rad.bam"
    bam2 = tmp_path / "wgs.bam"
    bam1.write_text("", encoding="utf-8")
    bam2.write_text("", encoding="utf-8")
    rename_bams = tmp_path / "rename.tsv"
    rename_bams.write_text(
        "# bam_basename sample_name\n"
        "rad.bam renamed_rad\n",
        encoding="utf-8",
    )

    rename_map = _normalize_bam_rename_file(rename_bams, [bam1, bam2])

    assert rename_map == {"rad.bam": "renamed_rad"}


def test_normalize_bam_rename_file_rejects_unknown_and_duplicate_inputs(tmp_path: Path) -> None:
    bam1 = tmp_path / "rad.bam"
    bam2 = tmp_path / "wgs.bam"
    bam1.write_text("", encoding="utf-8")
    bam2.write_text("", encoding="utf-8")

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text("missing.bam renamed\n", encoding="utf-8")
    with pytest.raises(IPyradError, match="not present in this assemble run: missing.bam"):
        _normalize_bam_rename_file(unknown, [bam1, bam2])

    dup_dir = tmp_path / "nested"
    dup_dir.mkdir()
    dup_bam = dup_dir / "rad.bam"
    dup_bam.write_text("", encoding="utf-8")
    rename_bams = tmp_path / "rename.tsv"
    rename_bams.write_text("rad.bam renamed\n", encoding="utf-8")
    with pytest.raises(IPyradError, match="input BAM basenames are duplicated: rad.bam"):
        _normalize_bam_rename_file(rename_bams, [bam1, dup_bam])


def test_get_vcf_with_indels_resolved_writes_stable_outputs_when_no_indels_exist(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir(parents=True)

    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + ("A" * 40) + "\n", encoding="utf-8")

    plain_vcf = vcf_dir / "loci.filtered.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=40>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t20\t.\tA\tT\t50\tPASS\t.\tGT\t1/1\t0/1\n",
        encoding="utf-8",
    )
    filtered_vcf = vcf_dir / "loci.filtered.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(filtered_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(filtered_vcf)]])

    resolved_vcf = get_vcf_with_indels_resolved(tmpdir, reference, threads=1)

    assert resolved_vcf == vcf_dir / "variants.resolved.vcf.gz"
    assert resolved_vcf.exists()
    assert resolved_vcf.with_suffix(resolved_vcf.suffix + ".csi").exists()
    assert get_indel_overlap_clusters_bed_path(tmpdir).read_text(encoding="utf-8") == ""
    assert get_indel_overlap_mask_path("s1", tmpdir).read_text(encoding="utf-8") == ""
    assert get_indel_overlap_mask_path("s2", tmpdir).read_text(encoding="utf-8") == ""
    assert load_variant_resolution_stats(tmpdir) == {
        "overlapping_indel_clusters_masked": 0,
        "overlapping_indel_records_removed": 0,
        "overlapping_indel_bp_masked": 0,
        "indel_records_inspected": 0,
    }


def test_run_assembler_uses_filtered_analysis_bams_downstream(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGTACGTAC\n", encoding="utf-8")
    rad_bam = tmp_path / "rad.bam"
    wgs_bam = tmp_path / "wgs.bam"
    rad_bam.write_text("", encoding="utf-8")
    wgs_bam.write_text("", encoding="utf-8")

    def _fake_get_name_from_bam(path: Path) -> str:
        return path.stem

    pool_calls: list[tuple[str, dict, int | None]] = []
    applied_paralog_masks: dict[str, Path] = {}
    masked_vcf: dict[str, Path] = {}
    built_database: dict[str, object] = {}
    seqs_hdf5_call: dict[str, object] = {}
    compacted_vcf: dict[str, Path] = {}
    variant_chunk_counts: list[int] = []
    variant_postfilter_stats: dict[str, int] = {}

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level
        pool_calls.append((msg, jobs, max_workers))
        if msg == "Filtering mapped reads":
            return {
                sname: tmp_path / "OUT" / "assembly_tmpdir" / "analysis_bams" / f"{sname}.analysis.filtered.bam"
                for sname in jobs
            }
        if msg == "Building sample-specific paralog masks":
            return {
                sname: tmp_path / "OUT" / "assembly_tmpdir" / "beds" / f"{sname}.paralog.mask.bed"
                for sname in jobs
            }
        if msg == "Summarizing final sample depths":
            return {
                sname: {
                    "shared_loci_with_nonzero_depth": 1,
                    "mean_depth_shared_loci": 3.0,
                    "median_depth_shared_loci": 3.0,
                    "mean_depth_nonzero_shared_loci": 3.0,
                    "median_depth_nonzero_shared_loci": 3.0,
                }
                for sname in jobs
            }
        if msg == "Merging final VCF masks":
            return {
                sname: tmp_path / "OUT" / "assembly_tmpdir" / "beds" / f"{sname}.final.vcf.mask.bed"
                for sname in jobs
            }
        return {sname: None for sname in jobs}

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", _fake_get_name_from_bam)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.classify_bam_layout",
        lambda path: "paired" if path.stem == "rad" else "single",
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_reference_sort_order",
        lambda _reference, tmpdir: (tmpdir / "REF_info.txt").write_text("chr1\t100\n", encoding="utf-8"),
    )

    def _fake_get_across_sample_loci_bed(*_args, **_kwargs):
        loci_bed = tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "loci.bed"
        loci_bed.write_text("chr1\t0\t10\t1\n", encoding="utf-8")
        return loci_bed

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_across_sample_loci_bed", _fake_get_across_sample_loci_bed)

    def _fake_aggregate_across_samples(*, regions_bed, sample_prefixes, in_dir, out_prefix, fail_frac_max, min_data_samples):
        del regions_bed, in_dir, fail_frac_max, min_data_samples
        Path(f"{out_prefix}.shared_good.final.bed").write_text("chr1\t0\t10\n", encoding="utf-8")
        Path(f"{out_prefix}.shared_good.strict_all_samples.bed").write_text("chr1\t0\t10\n", encoding="utf-8")
        sample_count = len(tuple(sample_prefixes))
        Path(f"{out_prefix}.shared_metrics.tsv").write_text(
            "chrom\tstart\tend\trid\tn_data\tn_good\tn_fail\tfail_frac_among_data\tgood_frac_among_data\tdrop_global\tkeep_global\n"
            f"chr1\t0\t10\tchr1:0-10\t{sample_count}\t{sample_count}\t0\t0.0\t1.0\tFalse\tTrue\n",
            encoding="utf-8",
        )
        return pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [0],
                "end": [10],
                "rid": ["chr1:0-10"],
                "n_data": [sample_count],
                "n_good": [sample_count],
                "n_fail": [0],
                "fail_frac_among_data": [0.0],
                "good_frac_among_data": [1.0],
                "drop_global": [False],
                "keep_global": [True],
            }
        )

    def _fake_write_per_sample_final_good(*, sample_prefixes, in_dir, shared_good_bed, out_dir, out_suffix=".final.good.bed"):
        del in_dir, shared_good_bed
        out_dir.mkdir(parents=True, exist_ok=True)
        for prefix in sample_prefixes:
            (out_dir / f"{prefix}{out_suffix}").write_text("chr1\t0\t10\n", encoding="utf-8")

    monkeypatch.setattr("ipyrad2.assembler.assemble.aggregate_across_samples", _fake_aggregate_across_samples)
    monkeypatch.setattr("ipyrad2.assembler.assemble.write_per_sample_final_good", _fake_write_per_sample_final_good)
    def _fake_get_chunked_loci_beds(tmpdir, nchunks, source_bed=None):
        variant_chunk_counts.append(nchunks)
        assert source_bed == tmpdir / "beds" / "loci.callable.variant.bed"
        chunk_bed = tmpdir / "beds" / "chunk-0.bed"
        chunk_bed.write_text("chr1\t0\t10\n", encoding="utf-8")
        return [chunk_bed]

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_chunked_loci_beds", _fake_get_chunked_loci_beds)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_concat_chunk_vcfs",
        lambda tmpdir, threads: (tmpdir / "vcfs" / "loci.raw.vcf.gz").write_text("", encoding="utf-8"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_filtered_vcf",
        lambda tmpdir, min_sample_depth, min_geno_q, min_site_q, threads: (
            tmpdir / "vcfs" / "loci.filtered.vcf.gz"
        ).write_text("", encoding="utf-8"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.apply_wgs_het_allele_balance_mask",
        lambda *args, **kwargs: {
            "wgs_het_genotypes_masked_by_allele_balance": 2,
            "wgs_het_genotypes_examined_for_allele_balance": 3,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_variant_postfilter_stats",
        lambda tmpdir, **stats: variant_postfilter_stats.update(stats),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.load_variant_postfilter_stats",
        lambda tmpdir: dict(variant_postfilter_stats),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble._load_mixed_paralog_counts",
        lambda tmpdir: {
            "loci_fail_paralog_rad": 0,
            "loci_fail_paralog_wgs": 0,
            "loci_fail_paralog_both": 0,
            "loci_pass_paralog_rad_fail_paralog_wgs": 0,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.summarize_variant_support_by_sample_type",
        lambda vcf_gz, rad_samples, wgs_samples: {
            "sites_supported_rad_only": 1,
            "sites_supported_wgs_only": 0,
            "sites_supported_both": 0,
            "sites_supported_neither": 0,
        },
    )
    def _fake_get_vcf_with_indels_resolved(tmpdir, reference, threads):
        del reference, threads
        path = tmpdir / "vcfs" / "variants.resolved.vcf.gz"
        path.write_text("", encoding="utf-8")
        return path

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_vcf_with_indels_resolved", _fake_get_vcf_with_indels_resolved)

    def _fake_write_sam_faidx(tmpdir):
        path = tmpdir / "loci.faidx.txt"
        path.write_text("chr1:1-10\n", encoding="utf-8")
        return path

    monkeypatch.setattr("ipyrad2.assembler.assemble.write_sam_faidx", _fake_write_sam_faidx)

    def _fake_get_reference_in_loci_beds(tmpdir, reference):
        del reference
        consensus_dir = tmpdir / "consensus_seqs"
        consensus_dir.mkdir(parents=True, exist_ok=True)
        path = consensus_dir / "assembly_reference_sequence.consensus.fa"
        path.write_text(">chr1:1-10\nACGT\n", encoding="utf-8")
        return path

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_reference_in_loci_beds", _fake_get_reference_in_loci_beds)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.apply_sample_region_masks_to_resolved_vcf",
        lambda tmpdir, sample_masks, vcf_gz=None: applied_paralog_masks.update(sample_masks) or masked_vcf.update({"path": vcf_gz}) or (vcf_gz if vcf_gz is not None else tmpdir / "vcfs" / "variants.resolved.vcf.gz"),
    )
    def _fake_write_final_outputs(
        *,
        snames,
        name,
        outdir,
        tmpdir,
        reference,
        min_locus_sample_coverage,
        min_locus_trim_sample_coverage,
        min_locus_length,
        max_locus_hetero_frequency,
        max_locus_variant_frequency,
        max_sample_hetero_frequency,
        cores,
        log_level,
    ):
        del reference
        del min_locus_sample_coverage
        del min_locus_trim_sample_coverage
        del min_locus_length
        del max_locus_hetero_frequency
        del max_locus_variant_frequency
        del max_sample_hetero_frequency
        del cores
        del log_level
        (outdir / f"{name}.bed").write_text("chr1\t0\t10\t1\n", encoding="utf-8")
        with gzip.open(outdir / f"{name}.loci.gz", "wt", encoding="utf-8") as out:
            out.write("// test\n")
        (outdir / f"{name}.hdf5").write_text("", encoding="utf-8")
        (outdir / "assembly_tmpdir" / "beds" / "rad.consensus_hetero.mask.bed").write_text("", encoding="utf-8")
        (outdir / "assembly_tmpdir" / "beds" / "wgs.consensus_hetero.mask.bed").write_text("", encoding="utf-8")
        summary = {
            "nloci_before_filtering": 1,
            "nloci_after_filtering": 1,
            "nsites_after_filtering": 4,
            "filter_counts": {
                "min_length": 0,
                "min_samples": 0,
                "max_variant_frequency": 0,
                "max_shared_hetero_frequency": 0,
                "max_depth_outlier": 0,
            },
            "site_totals": {
                "variant_sites": 1,
                "variant_phylo_informative_sites": 0,
                "nsites": 4,
                "nsites_sample_cov_greater_than_1": 4,
                "nsites_sample_cov_greater_than_2": 4,
                "nsites_sample_cov_greater_than_3": 0,
                "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 4,
            },
            "sample_locus_counts": {"rad": 1, "wgs": 1},
            "masked_by_max_hetero_frequency_counts": {"rad": 0, "wgs": 0},
            "loci_with_samples_masked_by_max_hetero_frequency": 0,
            "total_masked_sample_occurrences_by_max_hetero_frequency": 0,
            "samples_per_locus_counts": {2: 1},
            "locus_length_counts": {4: 1},
            "alignment_nonmissing_sample_bases": 8,
        }
        seqs_hdf5_call.update(
            {
                "snames": list(snames),
                "name": name,
                "outdir": outdir,
                "tmpdir": tmpdir,
                "loci_bed": outdir / f"{name}.bed",
                "nsites_after_filtering": summary["nsites_after_filtering"],
                "nloci_after_filtering": summary["nloci_after_filtering"],
            }
        )
        return summary

    monkeypatch.setattr("ipyrad2.assembler.assemble.write_final_outputs", _fake_write_final_outputs)
    def _fake_build_locus_fasta_database(name, snames, reference, tmpdir, masks):
        del reference
        built_database.update({"name": name, "snames": list(snames), "masks": masks})
        path = tmpdir / f"{name}.database.fa"
        path.write_text(">chr1:1-10 assembly_reference_sequence\nACGT\n", encoding="utf-8")
        return path

    monkeypatch.setattr("ipyrad2.assembler.assemble.build_locus_fasta_database", _fake_build_locus_fasta_database)
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.compact_resolved_vcf_to_final_loci_contigs",
        lambda tmpdir, reference, loci_bed: compacted_vcf.update({"loci_bed": loci_bed}) or (tmpdir / "vcfs" / "variants.resolved.vcf.gz"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_vcf",
        lambda name, outdir, tmpdir, threads: (
            (outdir / f"{name}.vcf.gz").write_text("", encoding="utf-8"),
            outdir / f"{name}.vcf.gz",
        )[1],
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.load_variant_resolution_stats",
        lambda tmpdir: {
            "overlapping_indel_clusters_masked": 0,
            "overlapping_indel_records_removed": 0,
            "overlapping_indel_bp_masked": 0,
            "indel_records_inspected": 0,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_snps_hdf5",
        lambda name, outdir, snames, reference: (
            (outdir / f"{name}.hdf5").write_text("", encoding="utf-8"),
            1,
        )[1],
    )
    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)

    run_assembler(
        rad_bams=[rad_bam],
        wgs_bams=[wgs_bam],
        reference=reference,
        outdir=tmp_path / "OUT",
        name="assembly",
        loci_bed=None,
        min_map_q=15,
        max_tlen=1500,
        max_softclip=20,
        max_nm=8,
        min_site_q=13,
        min_geno_q=13,
        min_base_q=13,
        min_sample_depth=1,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=25,
        min_locus_merge_distance=300,
        max_locus_hetero_frequency=0.3,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
        softclip_len_threshold=20,
        softclip_frac_max=0.5,
        depth_z_max=7.0,
        third_frac_cut=0.10,
        min_3allele_sites=2,
        maf_threshold=0.20,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.10,
        populations=None,
        rename_bams=None,
        masks=None,
        cores=4,
        threads=3,
        force=False,
        log_level="WARNING",
    )

    assert [msg for msg, _jobs, _max_workers in pool_calls] == [
        "Filtering mapped reads",
        "Delimiting sample coverage beds",
        "Filtering paralogs within samples",
        "Calling variants",
        "Building low-depth masks",
        "Building sample-specific paralog masks",
        "Merging sample masks",
        "Extracting consensus sequences",
        "Merging final VCF masks",
        "Summarizing final sample depths",
    ]

    assert pool_calls[0][2] == 1
    assert pool_calls[1][2] == 1
    assert pool_calls[2][2] == 1
    assert pool_calls[7][2] == 2
    assert pool_calls[8][2] == 1
    assert pool_calls[9][2] == 1

    filter_jobs = pool_calls[0][1]
    assert filter_jobs["rad"][0].__name__ == "prepare_filtered_analysis_bam"
    assert filter_jobs["wgs"][0].__name__ == "prepare_filtered_analysis_bam"
    assert filter_jobs["rad"][1]["is_paired"] is True
    assert filter_jobs["wgs"][1]["is_paired"] is False

    coverage_jobs = pool_calls[1][1]
    assert coverage_jobs["rad"][1]["bam_file"].name == "rad.analysis.filtered.bam"
    assert coverage_jobs["wgs"][1]["bam_file"].name == "wgs.analysis.filtered.bam"
    assert coverage_jobs["rad"][1]["is_paired"] is True
    assert coverage_jobs["wgs"][1]["is_paired"] is False

    paralog_jobs = pool_calls[2][1]
    assert paralog_jobs["rad"][1]["bam"].name == "rad.analysis.filtered.bam"
    assert paralog_jobs["wgs"][1]["bam"].name == "wgs.analysis.filtered.bam"

    variant_jobs = pool_calls[3][1]
    chunk_job = next(iter(variant_jobs.values()))
    assert [path.name for path in chunk_job[1]["bam_files"]] == [
        "rad.analysis.filtered.bam",
        "wgs.analysis.filtered.bam",
    ]
    assert variant_chunk_counts == [8]

    consensus_jobs = pool_calls[7][1]
    assert consensus_jobs["rad"][1]["reference_fasta"].name == "assembly_reference_sequence.consensus.fa"
    assert consensus_jobs["wgs"][1]["reference_fasta"].name == "assembly_reference_sequence.consensus.fa"
    final_vcf_mask_jobs = pool_calls[8][1]
    assert final_vcf_mask_jobs["rad"][0].__name__ == "merge_final_vcf_mask_beds"
    assert final_vcf_mask_jobs["wgs"][0].__name__ == "merge_final_vcf_mask_beds"
    assert seqs_hdf5_call == {
        "snames": ["rad", "wgs"],
        "name": "assembly",
        "outdir": tmp_path / "OUT",
        "tmpdir": tmp_path / "OUT" / "assembly_tmpdir",
        "loci_bed": tmp_path / "OUT" / "assembly.bed",
        "nsites_after_filtering": 4,
        "nloci_after_filtering": 1,
    }
    assert compacted_vcf == {"loci_bed": tmp_path / "OUT" / "assembly.bed"}

    assert applied_paralog_masks == {
        "rad": tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "rad.final.vcf.mask.bed",
        "wgs": tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "wgs.final.vcf.mask.bed",
    }
    assert masked_vcf == {"path": tmp_path / "OUT" / "assembly.vcf.gz"}
    assert built_database["snames"] == ["rad", "wgs"]
    assert variant_postfilter_stats == {
        "wgs_het_genotypes_masked_by_allele_balance": 2,
        "wgs_het_genotypes_examined_for_allele_balance": 3,
    }
    assert (tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "loci.paralog_filtered.bed").exists()
    assert (tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "loci.raw.bed").exists()
    assert (tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "loci.bed").read_text(encoding="utf-8") == "chr1\t0\t10\n"
    assert (tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "rad.final.good.bed").exists()
    assert (tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "wgs.final.good.bed").exists()
    assert (tmp_path / "OUT" / "assembly.loci.gz").exists()
    assert (tmp_path / "OUT" / "assembly.stats.txt").exists()
    assert "Mixed RAD/WGS Summary" in (tmp_path / "OUT" / "assembly.stats.txt").read_text(encoding="utf-8")
    assert (tmp_path / "OUT" / "assembly.vcf.gz").exists()
    assert (tmp_path / "OUT" / "assembly.hdf5").exists()


def test_run_assembler_rejects_duplicate_sample_names_across_rad_and_wgs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    rad_bam = tmp_path / "rad.bam"
    wgs_bam = tmp_path / "wgs.bam"
    rad_bam.write_text("", encoding="utf-8")
    wgs_bam.write_text("", encoding="utf-8")

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", lambda _path: "dup")

    with pytest.raises(IPyradError, match="duplicate sample names: dup"):
        run_assembler(
            rad_bams=[rad_bam],
            wgs_bams=[wgs_bam],
            reference=reference,
            outdir=tmp_path / "OUT",
            name="assembly",
            loci_bed=None,
            min_map_q=10,
            max_tlen=None,
            max_softclip=None,
            max_nm=None,
            min_site_q=13,
            min_geno_q=13,
            min_base_q=13,
            min_sample_depth=1,
            min_locus_sample_coverage=1,
            min_locus_trim_sample_coverage=1,
            min_locus_length=25,
            min_locus_merge_distance=300,
            max_locus_hetero_frequency=0.3,
            max_locus_variant_frequency=1.0,
            max_sample_hetero_frequency=0.10,
            softclip_len_threshold=20,
            softclip_frac_max=0.5,
            depth_z_max=7.0,
            third_frac_cut=0.10,
            min_3allele_sites=2,
            maf_threshold=0.20,
            max_sites_above_maf=8,
            paralog_fail_frac_max=0.10,
            populations=None,
            rename_bams=None,
            masks=None,
            cores=2,
            threads=1,
            force=False,
            log_level="WARNING",
        )


def test_run_assembler_rename_bams_overrides_header_names_for_populations_and_outputs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    rad_bam = tmp_path / "rad.bam"
    rad_bam.write_text("", encoding="utf-8")
    populations = tmp_path / "groups.tsv"
    populations.write_text("renamed_rad\tpop1\n", encoding="utf-8")
    rename_bams = tmp_path / "rename.tsv"
    rename_bams.write_text("rad.bam renamed_rad\n", encoding="utf-8")

    observed: dict[str, object] = {}

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", lambda _path: "header_name")
    monkeypatch.setattr("ipyrad2.assembler.assemble.classify_bam_layout", lambda _path: "paired")
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_reference_sort_order",
        lambda _reference, tmpdir: (tmpdir / "REF_info.txt").write_text("chr1\t4\n", encoding="utf-8"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble._prepare_analysis_bams",
        lambda **kwargs: {
            sname: kwargs["tmpdir"] / "analysis_bams" / f"{sname}.analysis.filtered.bam"
            for sname in kwargs["bam_dict"]
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.run_with_pool",
        lambda jobs, log_level, max_workers=None, msg="Processing": {sname: None for sname in jobs},
    )

    def _fake_get_across_sample_loci_bed(_snames, _mincov, _merge, _minlen, _suffix, tmpdir):
        loci_bed = tmpdir / "beds" / "rad.raw.bed"
        loci_bed.write_text("chr1\t0\t4\n", encoding="utf-8")
        return loci_bed

    def _fake_run_paralog_stage(**kwargs):
        final_bed = kwargs["bed_dir"] / "loci.bed"
        final_bed.write_text("chr1\t0\t4\n", encoding="utf-8")
        return final_bed

    def _fake_run_variant_stage(**kwargs):
        observed["group_samples_file"] = kwargs["group_samples_file"]
        observed["variant_bams"] = sorted(kwargs["bam_dict"])
        return kwargs["tmpdir"] / "vcfs" / "variants.resolved.vcf.gz"

    def _fake_write_consensus_and_outputs(**kwargs):
        observed["snames"] = kwargs["snames"]

    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_across_sample_loci_bed",
        _fake_get_across_sample_loci_bed,
    )
    monkeypatch.setattr("ipyrad2.assembler.assemble._run_paralog_stage", _fake_run_paralog_stage)
    monkeypatch.setattr("ipyrad2.assembler.assemble._run_variant_stage", _fake_run_variant_stage)
    monkeypatch.setattr("ipyrad2.assembler.assemble._build_sample_masks", lambda **_kwargs: {})
    monkeypatch.setattr("ipyrad2.assembler.assemble._write_consensus_and_outputs", _fake_write_consensus_and_outputs)

    run_assembler(
        rad_bams=[rad_bam],
        wgs_bams=None,
        reference=reference,
        outdir=tmp_path / "OUT",
        name="assembly",
        loci_bed=None,
        min_map_q=10,
        max_tlen=None,
        max_softclip=None,
        max_nm=None,
        min_site_q=13,
        min_geno_q=13,
        min_base_q=13,
        min_sample_depth=1,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=25,
        min_locus_merge_distance=300,
        max_locus_hetero_frequency=0.3,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
        softclip_len_threshold=20,
        softclip_frac_max=0.5,
        depth_z_max=7.0,
        third_frac_cut=0.10,
        min_3allele_sites=2,
        maf_threshold=0.20,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.10,
        populations=populations,
        rename_bams=rename_bams,
        masks=None,
        cores=2,
        threads=1,
        force=False,
        log_level="WARNING",
    )

    group_samples_file = observed["group_samples_file"]
    assert group_samples_file == tmp_path / "OUT" / "assembly_tmpdir" / "populations.normalized.tsv"
    assert group_samples_file.read_text(encoding="utf-8") == "renamed_rad\tpop1\n"
    assert observed["variant_bams"] == ["renamed_rad"]
    assert observed["snames"] == ["renamed_rad"]


def test_run_assembler_requires_rad_bams_when_no_loci_bed(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    wgs_bam = tmp_path / "wgs.bam"
    wgs_bam.write_text("", encoding="utf-8")

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", lambda path: path.stem)

    with pytest.raises(IPyradError, match="No RAD bam files found. These are required unless --loci-bed is provided."):
        run_assembler(
            rad_bams=None,
            wgs_bams=[wgs_bam],
            reference=reference,
            outdir=tmp_path / "OUT",
            name="assembly",
            loci_bed=None,
            min_map_q=10,
            max_tlen=None,
            max_softclip=None,
            max_nm=None,
            min_site_q=13,
            min_geno_q=13,
            min_base_q=13,
            min_sample_depth=1,
            min_locus_sample_coverage=1,
            min_locus_trim_sample_coverage=1,
            min_locus_length=25,
            min_locus_merge_distance=300,
            max_locus_hetero_frequency=0.3,
            max_locus_variant_frequency=1.0,
            max_sample_hetero_frequency=0.10,
            softclip_len_threshold=20,
            softclip_frac_max=0.5,
            depth_z_max=7.0,
            third_frac_cut=0.10,
            min_3allele_sites=2,
            maf_threshold=0.20,
            max_sites_above_maf=8,
            paralog_fail_frac_max=0.10,
            populations=None,
            rename_bams=None,
            masks=None,
            cores=2,
            threads=1,
            force=False,
            log_level="WARNING",
        )


def test_run_assembler_requires_at_least_one_bam_with_loci_bed(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    loci_bed = tmp_path / "loci.bed"
    loci_bed.write_text("chr1\t0\t4\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="No input BAM files found. Provide --rad-bams and/or --wgs-bams."):
        run_assembler(
            rad_bams=None,
            wgs_bams=None,
            reference=reference,
            outdir=tmp_path / "OUT",
            name="assembly",
            loci_bed=loci_bed,
            min_map_q=10,
            max_tlen=None,
            max_softclip=None,
            max_nm=None,
            min_site_q=13,
            min_geno_q=13,
            min_base_q=13,
            min_sample_depth=1,
            min_locus_sample_coverage=1,
            min_locus_trim_sample_coverage=1,
            min_locus_length=25,
            min_locus_merge_distance=300,
            max_locus_hetero_frequency=0.3,
            max_locus_variant_frequency=1.0,
            max_sample_hetero_frequency=0.10,
            softclip_len_threshold=20,
            softclip_frac_max=0.5,
            depth_z_max=7.0,
            third_frac_cut=0.10,
            min_3allele_sites=2,
            maf_threshold=0.20,
            max_sites_above_maf=8,
            paralog_fail_frac_max=0.10,
            populations=None,
            rename_bams=None,
            masks=None,
            cores=2,
            threads=1,
            force=False,
            log_level="WARNING",
        )


def test_run_assembler_accepts_loci_bed_without_rad_samples(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n>chr2\nACGT\n", encoding="utf-8")
    wgs_bam = tmp_path / "wgs.bam"
    wgs_bam.write_text("", encoding="utf-8")
    loci_bed = tmp_path / "input.bed"
    loci_bed.write_text("chr1\t5\t10\tx\nchr2\t0\t4\ty\n", encoding="utf-8")

    observed: dict[str, object] = {}
    coverage_job_names: list[str] = []

    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", lambda path: path.stem)
    monkeypatch.setattr("ipyrad2.assembler.assemble.classify_bam_layout", lambda _path: "single")
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_reference_sort_order",
        lambda _reference, tmpdir: (tmpdir / "REF_info.txt").write_text("chr2\t4\nchr1\t4\n", encoding="utf-8"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble._prepare_analysis_bams",
        lambda **kwargs: {
            sname: kwargs["tmpdir"] / "analysis_bams" / f"{sname}.analysis.filtered.bam"
            for sname in kwargs["bam_dict"]
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_across_sample_loci_bed",
        lambda *_args, **_kwargs: pytest.fail("RAD-based locus delimiting should be skipped when --loci-bed is provided"),
    )

    def _fake_run_with_pool(jobs, log_level, max_workers=None, msg="Processing"):
        del log_level, max_workers
        if msg == "Delimiting sample coverage beds":
            coverage_job_names.extend(sorted(jobs))
        return {sname: None for sname in jobs}

    def _fake_run_paralog_stage(**kwargs):
        observed["regions_bed"] = kwargs["regions_bed"]
        final_bed = kwargs["bed_dir"] / "loci.bed"
        final_bed.write_text("chr2\t0\t4\nchr1\t5\t10\n", encoding="utf-8")
        return final_bed

    def _fake_run_variant_stage(**kwargs):
        observed["variant_bams"] = sorted(kwargs["bam_dict"])
        return kwargs["tmpdir"] / "vcfs" / "variants.resolved.vcf.gz"

    def _fake_write_consensus_and_outputs(**kwargs):
        observed["snames"] = kwargs["snames"]
        observed["shared_loci_after_delimiting"] = kwargs["shared_loci_after_delimiting"]
        observed["shared_loci_after_paralog_filtering"] = kwargs["shared_loci_after_paralog_filtering"]

    monkeypatch.setattr("ipyrad2.assembler.assemble.run_with_pool", _fake_run_with_pool)
    monkeypatch.setattr("ipyrad2.assembler.assemble._run_paralog_stage", _fake_run_paralog_stage)
    monkeypatch.setattr("ipyrad2.assembler.assemble._run_variant_stage", _fake_run_variant_stage)
    monkeypatch.setattr("ipyrad2.assembler.assemble._build_sample_masks", lambda **_kwargs: {})
    monkeypatch.setattr("ipyrad2.assembler.assemble._write_consensus_and_outputs", _fake_write_consensus_and_outputs)

    run_assembler(
        rad_bams=None,
        wgs_bams=[wgs_bam],
        reference=reference,
        outdir=tmp_path / "OUT",
        name="assembly",
        loci_bed=loci_bed,
        min_map_q=10,
        max_tlen=None,
        max_softclip=None,
        max_nm=None,
        min_site_q=13,
        min_geno_q=13,
        min_base_q=13,
        min_sample_depth=1,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=1,
        min_locus_length=25,
        min_locus_merge_distance=300,
        max_locus_hetero_frequency=0.3,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
        softclip_len_threshold=20,
        softclip_frac_max=0.5,
        depth_z_max=7.0,
        third_frac_cut=0.10,
        min_3allele_sites=2,
        maf_threshold=0.20,
        max_sites_above_maf=8,
        paralog_fail_frac_max=0.10,
        populations=None,
        rename_bams=None,
        masks=None,
        cores=2,
        threads=1,
        force=False,
        log_level="WARNING",
    )

    assert coverage_job_names == ["wgs"]
    assert observed["variant_bams"] == ["wgs"]
    assert observed["snames"] == ["wgs"]
    assert observed["shared_loci_after_delimiting"] == 2
    assert observed["shared_loci_after_paralog_filtering"] == 2
    assert observed["regions_bed"] == tmp_path / "OUT" / "assembly_tmpdir" / "beds" / "loci.raw.bed"
    assert observed["regions_bed"].read_text(encoding="utf-8") == "chr2\t0\t4\nchr1\t5\t10\n"


def test_filter_trim_locus_respects_min_locus_length() -> None:
    header = "chr1:1-20"
    locus_dict = {
        "assembly_reference_sequence": "A" * 20,
        "s1": "A" * 20,
        "s2": "A" * 20,
        "s3": "A" * 20,
    }

    result_short = filter_trim_locus(
        header,
        locus_dict,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=4,
        min_locus_length=20,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
    )
    result_long = filter_trim_locus(
        header,
        locus_dict,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=4,
        min_locus_length=25,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
    )

    assert result_short[4]["min_length"] is False
    assert result_long[4]["min_length"] is True


def test_write_loci_and_stats_files_counts_max_variant_frequency_filter(
    tmp_path: Path,
) -> None:
    database = tmp_path / "assembly.database.fa"
    database.write_text(
        ">chr1:1-4 assembly_reference_sequence\nAAAA\n"
        ">chr1:1-4 s1\nAAAA\n"
        ">chr1:1-4 s2\nAATA\n"
        ">chr1:1-4 s3\nAATA\n\n",
        encoding="utf-8",
    )
    summary = write_loci_and_stats_files(
        snames=["s1", "s2", "s3"],
        name="assembly",
        outdir=tmp_path,
        tmpdir=tmp_path,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=4,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=0.2,
    )

    assert summary["nloci_after_filtering"] == 0
    assert summary["nsites_after_filtering"] == 0
    assert summary["filter_counts"]["max_variant_frequency"] == 1
    assert not (tmp_path / "assembly.stats_counts.tsv").exists()


def test_write_loci_and_stats_files_writes_gzipped_loci_and_streamed_bed(
    tmp_path: Path,
) -> None:
    database = tmp_path / "assembly.database.fa"
    database.write_text(
        ">chr1:1-4 assembly_reference_sequence\nAAAA\n"
        ">chr1:1-4 s1\nAAAA\n"
        ">chr1:1-4 s2\nAATA\n"
        ">chr1:1-4 s3\nAATA\n\n",
        encoding="utf-8",
    )

    summary = write_loci_and_stats_files(
        snames=["s1", "s2", "s3"],
        name="assembly",
        outdir=tmp_path,
        tmpdir=tmp_path,
        min_locus_sample_coverage=4,
        min_locus_trim_sample_coverage=4,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
    )

    with gzip.open(tmp_path / "assembly.loci.gz", "rt", encoding="utf-8") as handle:
        loci_lines = handle.read().splitlines()

    assert summary["nloci_after_filtering"] == 1
    assert summary["nsites_after_filtering"] == 4
    assert all(not line.startswith(">") for line in loci_lines)
    assert loci_lines[0].startswith("assembly_reference_sequence")
    assert loci_lines[0].endswith("AAAA")
    assert loci_lines[1].startswith("s1")
    assert loci_lines[1].endswith("AAAA")
    assert loci_lines[2].startswith("s2")
    assert loci_lines[2].endswith("AATA")
    assert loci_lines[3].startswith("s3")
    assert loci_lines[3].endswith("AATA")
    assert loci_lines[4].startswith("//")
    assert loci_lines[4].endswith("|0:chr1:1-4")
    assert (tmp_path / "assembly.bed").read_text(encoding="utf-8") == "chr1\t0\t4\t4\n"
    assert summary["sample_locus_counts"] == {"s1": 1, "s2": 1, "s3": 1}
    assert summary["samples_per_locus_counts"] == {3: 1}
    assert summary["locus_length_counts"] == {4: 1}
    assert not (tmp_path / "assembly.stats_counts.tsv").exists()
    assert not (tmp_path / "assembly.stats_sample_cov.txt").exists()
    assert not (tmp_path / "assembly.stats_locus_coverage.txt").exists()


def test_write_loci_and_stats_files_masks_samples_above_max_sample_hetero_frequency(
    tmp_path: Path,
) -> None:
    database = tmp_path / "assembly.database.fa"
    database.write_text(
        ">chr1:1-4 assembly_reference_sequence\nAAAA\n"
        ">chr1:1-4 s1\nARAA\n"
        ">chr1:1-4 s2\nAAAA\n"
        ">chr1:1-4 s3\nAAAA\n\n",
        encoding="utf-8",
    )

    summary = write_loci_and_stats_files(
        snames=["s1", "s2", "s3"],
        name="assembly",
        outdir=tmp_path,
        tmpdir=tmp_path,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.20,
    )

    with gzip.open(tmp_path / "assembly.loci.gz", "rt", encoding="utf-8") as handle:
        loci_lines = handle.read().splitlines()

    assert summary["nloci_after_filtering"] == 1
    assert summary["sample_locus_counts"] == {"s1": 0, "s2": 1, "s3": 1}
    assert summary["masked_by_max_hetero_frequency_counts"] == {"s1": 1, "s2": 0, "s3": 0}
    assert summary["loci_with_samples_masked_by_max_hetero_frequency"] == 1
    assert summary["total_masked_sample_occurrences_by_max_hetero_frequency"] == 1
    assert all(not line.startswith("s1") for line in loci_lines)
    assert loci_lines[1].startswith("s2")
    assert loci_lines[2].startswith("s3")
    assert loci_lines[3].startswith("//")
    assert (tmp_path / "assembly.bed").read_text(encoding="utf-8") == "chr1\t0\t4\t3\n"
    assert (tmp_path / "beds" / "s1.consensus_hetero.mask.bed").read_text(encoding="utf-8") == "chr1\t0\t4\n"
    manifest = (tmp_path / "assembly.retained_loci.tsv").read_text(encoding="utf-8")
    assert "chr1:1-4\tchr1:1-4\ts1" in manifest


def test_filter_trim_locus_counts_masked_n_at_variable_sites_for_sample_hetero_filter() -> None:
    header = "chr1:1-4"
    locus_dict = {
        "assembly_reference_sequence": "AAAA",
        "s1": "ANAA",
        "s2": "ATAA",
        "s3": "AAAA",
    }

    _header, _names, tseqs, _snps, filters, stats = filter_trim_locus(
        header,
        locus_dict,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.20,
    )

    assert not any(filters.values())
    assert stats["masked_samples_by_max_sample_hetero_frequency"] == ("s1",)
    assert bytes(tseqs[1]).decode() == "NNNN"


def test_get_sample_depth_stats_in_final_loci_uses_full_shared_locus_length(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True)
    loci_bed = tmp_path / "assembly.bed"
    loci_bed.write_text(
        "chr1\t0\t10\t2\n"
        "chr1\t20\t30\t2\n"
        "chr1\t40\t50\t2\n",
        encoding="utf-8",
    )
    (bed_dir / "s1.fragments.bedgraph").write_text(
        "chr1\t0\t5\t2\n"
        "chr1\t5\t10\t4\n"
        "chr1\t20\t25\t1\n",
        encoding="utf-8",
    )

    stats = get_sample_depth_stats_in_final_loci("s1", loci_bed, tmpdir)

    assert stats["shared_loci_with_nonzero_depth"] == 2
    assert stats["mean_depth_shared_loci"] == pytest.approx((3.0 + 0.5 + 0.0) / 3.0)
    assert stats["median_depth_shared_loci"] == pytest.approx(0.5)
    assert stats["mean_depth_nonzero_shared_loci"] == pytest.approx(1.75)
    assert stats["median_depth_nonzero_shared_loci"] == pytest.approx(1.75)


def test_write_assemble_stats_report_writes_single_text_report(tmp_path: Path) -> None:
    outpath = write_assemble_stats_report(
        name="assembly",
        outdir=tmp_path,
        snames=["s1", "s2"],
        shared_loci_after_delimiting=10,
        shared_loci_after_paralog_filtering=8,
        loci_summary={
            "nloci_after_filtering": 6,
            "nsites_after_filtering": 24,
            "filter_counts": {
                "min_length": 1,
                "min_samples": 0,
                "max_variant_frequency": 2,
                "max_shared_hetero_frequency": 0,
                "max_depth_outlier": 0,
            },
            "site_totals": {
                "variant_sites": 5,
                "variant_phylo_informative_sites": 2,
                "nsites": 24,
                "nsites_sample_cov_greater_than_1": 20,
                "nsites_sample_cov_greater_than_2": 18,
                "nsites_sample_cov_greater_than_3": 0,
                "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 18,
            },
            "sample_locus_counts": {"s1": 5, "s2": 4},
            "masked_by_max_hetero_frequency_counts": {"s1": 2, "s2": 0},
            "loci_with_samples_masked_by_max_hetero_frequency": 2,
            "total_masked_sample_occurrences_by_max_hetero_frequency": 2,
            "samples_per_locus_counts": {1: 3, 2: 3},
            "locus_length_counts": {4: 2, 5: 4},
            "alignment_nonmissing_sample_bases": 36,
        },
        sample_depth_stats={
            "s1": {
                "shared_loci_with_nonzero_depth": 6,
                "mean_depth_shared_loci": 8.0,
                "median_depth_shared_loci": 7.5,
                "mean_depth_nonzero_shared_loci": 8.0,
                "median_depth_nonzero_shared_loci": 7.5,
            },
            "s2": {
                "shared_loci_with_nonzero_depth": 4,
                "mean_depth_shared_loci": 4.0,
                "median_depth_shared_loci": 3.5,
                "mean_depth_nonzero_shared_loci": 6.0,
                "median_depth_nonzero_shared_loci": 6.0,
            },
        },
        nsnps_written=5,
        overlap_stats={
            "overlapping_indel_clusters_masked": 2,
            "overlapping_indel_records_removed": 4,
            "overlapping_indel_bp_masked": 9,
            "indel_records_inspected": 12,
        },
    )

    report = outpath.read_text(encoding="utf-8")
    assert outpath == tmp_path / "assembly.stats.txt"
    assert "# Assemble Summary" in report
    assert "# Locus Filtering" in report
    assert "# Sample Masking" in report
    assert "# Alignment Summary" in report
    assert "# Sample Summary" in report
    assert "# Locus Occupancy" in report
    assert "shared_loci_after_delimiting" in report
    assert "final_snp_sites_written" in report
    assert "masked_by_max_hetero_frequency" in report
    assert "loci_with_samples_masked_by_max_hetero_frequency" in report
    assert "total_masked_sample_occurrences_by_max_hetero_frequency" in report
    assert "loci_with_sample_masks_max_sample_heterozygosity" not in report
    assert "sample_locus_masks_max_sample_heterozygosity" not in report
    assert "loci_masked_max_sample_hetero_frequency" not in report
    assert "s1" in report
    assert "assembly_reference_sequence" not in report


def test_write_assemble_stats_report_includes_mixed_run_summary(tmp_path: Path) -> None:
    outpath = write_assemble_stats_report(
        name="assembly",
        outdir=tmp_path,
        snames=["rad", "wgs"],
        shared_loci_after_delimiting=10,
        shared_loci_after_paralog_filtering=8,
        loci_summary={
            "nloci_after_filtering": 6,
            "nsites_after_filtering": 24,
            "filter_counts": {
                "min_length": 1,
                "min_samples": 0,
                "max_variant_frequency": 2,
                "max_shared_hetero_frequency": 0,
                "max_depth_outlier": 0,
            },
            "site_totals": {
                "variant_sites": 5,
                "variant_phylo_informative_sites": 2,
                "nsites": 24,
                "nsites_sample_cov_greater_than_1": 20,
                "nsites_sample_cov_greater_than_2": 18,
                "nsites_sample_cov_greater_than_3": 0,
                "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 18,
            },
            "sample_locus_counts": {"rad": 5, "wgs": 4},
            "samples_per_locus_counts": {1: 3, 2: 3},
            "locus_length_counts": {4: 2, 5: 4},
            "alignment_nonmissing_sample_bases": 36,
        },
        sample_depth_stats={
            "rad": {
                "shared_loci_with_nonzero_depth": 6,
                "mean_depth_shared_loci": 8.0,
                "median_depth_shared_loci": 7.5,
                "mean_depth_nonzero_shared_loci": 8.0,
                "median_depth_nonzero_shared_loci": 7.5,
            },
            "wgs": {
                "shared_loci_with_nonzero_depth": 4,
                "mean_depth_shared_loci": 4.0,
                "median_depth_shared_loci": 3.5,
                "mean_depth_nonzero_shared_loci": 6.0,
                "median_depth_nonzero_shared_loci": 6.0,
            },
        },
        nsnps_written=5,
        overlap_stats={
            "overlapping_indel_clusters_masked": 2,
            "overlapping_indel_records_removed": 4,
            "overlapping_indel_bp_masked": 9,
            "indel_records_inspected": 12,
        },
        mixed_run_summary={
            "rad_samples": 1,
            "wgs_samples": 1,
            "loci_fail_paralog_rad": 1,
            "loci_fail_paralog_wgs": 2,
            "loci_fail_paralog_both": 1,
            "loci_pass_paralog_rad_fail_paralog_wgs": 1,
            "sites_supported_rad_only": 2,
            "sites_supported_wgs_only": 1,
            "sites_supported_both": 2,
            "wgs_het_genotypes_masked_by_allele_balance": 3,
        },
    )

    report = outpath.read_text(encoding="utf-8")
    assert "# Mixed RAD/WGS Summary" in report
    assert "sites_supported_wgs_only" in report
    assert "wgs_het_genotypes_masked_by_allele_balance" in report


def test_write_consensus_and_outputs_fails_cleanly_when_no_loci_survive(
    monkeypatch,
    tmp_path: Path,
) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()

    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_sam_faidx",
        lambda _tmpdir: _tmpdir / "loci.faidx.txt",
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.get_reference_in_loci_beds",
        lambda _tmpdir, _reference: _tmpdir / "consensus_seqs" / "assembly_reference_sequence.consensus.fa",
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.build_locus_fasta_database",
        lambda *args, **kwargs: tmpdir / "assembly.database.fa",
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_final_outputs",
        lambda **kwargs: {
            "nloci_before_filtering": 5,
            "nloci_after_filtering": 0,
            "nsites_after_filtering": 0,
            "filter_counts": {},
            "site_totals": {},
            "sample_locus_counts": {},
            "samples_per_locus_counts": {},
            "locus_length_counts": {},
            "alignment_nonmissing_sample_bases": 0,
        },
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.run_with_pool",
        lambda jobs, log_level, max_workers=None, msg="Processing": {name: None for name in jobs},
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_vcf",
        lambda *args, **kwargs: pytest.fail("write_vcf should not run when no loci survive"),
    )
    monkeypatch.setattr(
        "ipyrad2.assembler.assemble.write_snps_hdf5",
        lambda *args, **kwargs: pytest.fail("write_snps_hdf5 should not run when no loci survive"),
    )

    with pytest.raises(IPyradError, match="No loci passed final trimming/filtering"):
        _write_consensus_and_outputs(
            name="assembly",
            outdir=tmp_path,
            tmpdir=tmpdir,
            snames=["s1"],
            reference=reference,
            masks=None,
            sample_masks=None,
            shared_loci_after_delimiting=5,
            shared_loci_after_paralog_filtering=5,
            min_locus_sample_coverage=1,
            min_locus_trim_sample_coverage=1,
            min_locus_length=25,
            max_locus_hetero_frequency=0.3,
            max_locus_variant_frequency=1.0,
            max_sample_hetero_frequency=0.10,
            consensus_workers=1,
            workers=1,
            threads=1,
            log_level="WARNING",
        )


def test_write_snps_hdf5_supports_zero_snp_outputs(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference.with_suffix(".fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    out_h5 = tmp_path / "assembly.hdf5"
    out_h5.touch()
    (tmp_path / "assembly.bed").write_text("chr1\t0\t4\n", encoding="utf-8")
    with gzip.open(tmp_path / "assembly.vcf.gz", "wt", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
    nsnps = write_snps_hdf5("assembly", tmp_path, ["s1"], reference)
    assert nsnps == 0

    with h5py.File(out_h5, "r") as io5:
        assert int(io5.attrs["nsnps"]) == 0
        assert io5["snpsmap"].shape == (0, 5)
        assert io5["snpsmap"].dtype == np.uint32
        assert io5["genos"].shape == (1, 0, 3)
        assert io5["genos"].chunks[1] <= 131_072
        assert io5["snpsmap"].chunks[0] <= 131_072
        assert io5["reference"].shape == (0,)


def test_write_seqs_hdf5_completes_without_tail_debug_code(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")
    reference.with_suffix(".fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    loci_bed = tmp_path / "assembly.bed"
    loci_bed.write_text("chr1\t0\t4\n", encoding="utf-8")

    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()
    database = tmpdir / "assembly.database.fa"
    database.write_text(
        ">chr1:1-4 assembly_reference_sequence\nAAAA\n"
        ">chr1:1-4 s1\nAAAA\n"
        ">chr1:1-4 s2\nAATA\n"
        ">chr1:1-4 s3\nAATA\n\n",
        encoding="utf-8",
    )
    (tmpdir / "assembly.retained_loci.tsv").write_text(
        "raw_header\tfinal_header\tmasked_samples\n"
        "chr1:1-4\tchr1:1-4\t\n",
        encoding="utf-8",
    )

    write_seqs_hdf5(
        name="assembly",
        outdir=tmp_path,
        tmpdir=tmpdir,
        snames=["s1", "s2", "s3"],
        reference=reference,
        loci_bed=loci_bed,
        nsites_after_filtering=4,
        nloci_after_filtering=1,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
    )

    with h5py.File(tmp_path / "assembly.hdf5", "r") as io5:
        assert io5["phy"].shape == (4, 4)
        assert io5["phymap"].shape == (1, 5)
        assert io5["phy"].chunks[1] <= 262_144
        assert io5["phymap"].dtype == np.uint32


def test_write_seqs_hdf5_compacts_scaffold_metadata_to_retained_bed_subset(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr1\nAAAA\n"
        ">chr2\nCCCC\n"
        ">chr3\nGGGG\n",
        encoding="utf-8",
    )
    reference.with_suffix(".fa.fai").write_text(
        "chr1\t4\t6\t4\t5\n"
        "chr2\t4\t17\t4\t5\n"
        "chr3\t4\t28\t4\t5\n",
        encoding="utf-8",
    )
    loci_bed = tmp_path / "assembly.bed"
    loci_bed.write_text(
        "chr3\t4\t8\n"
        "chr2\t0\t4\n",
        encoding="utf-8",
    )

    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()
    database = tmpdir / "assembly.database.fa"
    database.write_text(
        ">chr2:1-4 assembly_reference_sequence\nCCCC\n"
        ">chr2:1-4 s1\nCCCC\n"
        ">chr2:1-4 s2\nCCCC\n"
        ">chr2:1-4 s3\nCCCC\n\n"
        ">chr3:5-8 assembly_reference_sequence\nGGGG\n"
        ">chr3:5-8 s1\nGGGG\n"
        ">chr3:5-8 s2\nGGGG\n"
        ">chr3:5-8 s3\nGGGG\n\n",
        encoding="utf-8",
    )
    (tmpdir / "assembly.retained_loci.tsv").write_text(
        "raw_header\tfinal_header\tmasked_samples\n"
        "chr2:1-4\tchr2:1-4\t\n"
        "chr3:5-8\tchr3:5-8\t\n",
        encoding="utf-8",
    )

    write_seqs_hdf5(
        name="assembly",
        outdir=tmp_path,
        tmpdir=tmpdir,
        snames=["s1", "s2", "s3"],
        reference=reference,
        loci_bed=loci_bed,
        nsites_after_filtering=8,
        nloci_after_filtering=2,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
    )

    with h5py.File(tmp_path / "assembly.hdf5", "r") as io5:
        assert list(io5.attrs["scaffold_names"]) == ["chr2", "chr3"]
        assert list(io5.attrs["scaffold_lengths"]) == [4, 4]
        np.testing.assert_array_equal(
            io5["phymap"][:],
            np.array(
                [
                    [0, 0, 4, 1, 4],
                    [1, 4, 8, 5, 8],
                ],
                dtype=np.uint32,
            ),
        )


def test_write_final_outputs_finalizes_hdf5_after_bed_is_complete(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr1\nAAAA\n"
        ">chr2\nCCCC\n"
        ">chr3\nGGGG\n",
        encoding="utf-8",
    )
    reference.with_suffix(".fa.fai").write_text(
        "chr1\t4\t6\t4\t5\n"
        "chr2\t4\t17\t4\t5\n"
        "chr3\t4\t28\t4\t5\n",
        encoding="utf-8",
    )

    outdir = tmp_path / "out"
    outdir.mkdir()
    tmpdir = outdir / "assembly_tmpdir"
    tmpdir.mkdir()
    (tmpdir / "beds").mkdir()
    database = tmpdir / "assembly.database.fa"
    database.write_text(
        ">chr2:1-4 assembly_reference_sequence\nCCCC\n"
        ">chr2:1-4 s1\nCCCC\n"
        ">chr2:1-4 s2\nCCCC\n\n"
        ">chr3:5-8 assembly_reference_sequence\nGGGG\n"
        ">chr3:5-8 s1\nGGGG\n"
        ">chr3:5-8 s2\nGGGG\n\n",
        encoding="utf-8",
    )

    summary = write_final_outputs(
        snames=["s1", "s2"],
        name="assembly",
        outdir=outdir,
        tmpdir=tmpdir,
        reference=reference,
        min_locus_sample_coverage=1,
        min_locus_trim_sample_coverage=1,
        min_locus_length=1,
        max_locus_hetero_frequency=1.0,
        max_locus_variant_frequency=1.0,
        max_sample_hetero_frequency=0.10,
        cores=2,
        log_level="WARNING",
        batch_size=1,
    )

    assert summary["nloci_after_filtering"] == 2
    assert (outdir / "assembly.bed").read_text(encoding="utf-8") == (
        "chr2\t0\t4\t3\n"
        "chr3\t4\t8\t3\n"
    )
    with h5py.File(outdir / "assembly.hdf5", "r") as io5:
        assert list(io5.attrs["scaffold_names"]) == ["chr2", "chr3"]
        assert list(io5.attrs["scaffold_lengths"]) == [4, 4]
        np.testing.assert_array_equal(
            io5["phymap"][:],
            np.array(
                [
                    [0, 0, 4, 1, 4],
                    [1, 4, 8, 5, 8],
                ],
                dtype=np.uint32,
            ),
        )


def test_compact_resolved_vcf_to_final_loci_contigs_trims_contig_headers(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr1\nAAAA\n"
        ">chr2\nCCCC\n"
        ">chr3\nGGGG\n",
        encoding="utf-8",
    )
    reference.with_suffix(".fa.fai").write_text(
        "chr1\t4\t6\t4\t5\n"
        "chr2\t4\t17\t4\t5\n"
        "chr3\t4\t28\t4\t5\n",
        encoding="utf-8",
    )

    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    vcf_dir.mkdir(parents=True)
    plain_vcf = vcf_dir / "variants.resolved.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=4>\n"
        "##contig=<ID=chr2,length=4>\n"
        "##contig=<ID=chr3,length=4>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr2\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n"
        "chr3\t2\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\n",
        encoding="utf-8",
    )
    resolved_vcf = vcf_dir / "variants.resolved.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(resolved_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(resolved_vcf)]])

    outdir = tmp_path
    loci_bed = outdir / "assembly.bed"
    loci_bed.write_text(
        "chr3\t0\t4\n"
        "chr2\t0\t4\n",
        encoding="utf-8",
    )

    compact_resolved_vcf_to_final_loci_contigs(tmpdir, reference, loci_bed)
    write_vcf("assembly", outdir, tmpdir, threads=1)

    with gzip.open(resolved_vcf, "rt", encoding="utf-8") as handle:
        resolved_lines = handle.readlines()
    assert [line.rstrip("\n") for line in resolved_lines if line.startswith("##contig=")] == [
        "##contig=<ID=chr2,length=4>",
        "##contig=<ID=chr3,length=4>",
    ]

    final_vcf = outdir / "assembly.vcf.gz"
    with gzip.open(final_vcf, "rt", encoding="utf-8") as handle:
        final_lines = handle.readlines()
    assert [line.rstrip("\n") for line in final_lines if line.startswith("##contig=")] == [
        "##contig=<ID=chr2,length=4>",
        "##contig=<ID=chr3,length=4>",
    ]
    assert [
        line.rstrip("\n")
        for line in final_lines
        if line and not line.startswith("#")
    ] == [
        "chr2\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1",
        "chr3\t2\t.\tG\tA\t50\tPASS\t.\tGT\t1/1",
    ]


def test_apply_sample_region_masks_to_resolved_vcf_masks_only_targeted_sample(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir(parents=True)

    plain_vcf = vcf_dir / "variants.resolved.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\t0/1\n",
        encoding="utf-8",
    )
    final_vcf = tmp_path / "assembly.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(final_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(final_vcf)]])

    mask_bed = bed_dir / "s1.paralog.mask.bed"
    mask_bed.write_text("chr1\t0\t15\n", encoding="utf-8")

    apply_sample_region_masks_to_resolved_vcf(tmpdir, {"s1": mask_bed}, vcf_gz=final_vcf)

    with gzip.open(final_vcf, "rt", encoding="utf-8") as handle:
        rows = [
            line.rstrip("\n").split("\t")
            for line in handle
            if line and not line.startswith("#")
        ]

    assert rows[0][9] == "./."
    assert rows[0][10] == "0/0"
    assert rows[1][9] == "1/1"
    assert rows[1][10] == "0/1"


def test_apply_sample_region_masks_to_resolved_vcf_rejects_unknown_samples(
    tmp_path: Path,
) -> None:
    tmpdir = tmp_path / "assembly_tmpdir"
    vcf_dir = tmpdir / "vcfs"
    bed_dir = tmpdir / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir(parents=True)

    plain_vcf = vcf_dir / "variants.resolved.vcf"
    plain_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
        encoding="utf-8",
    )
    final_vcf = tmp_path / "assembly.vcf.gz"
    run_pipeline([[BIN_BCF, "view", "-Oz", "-o", str(final_vcf), str(plain_vcf)]])
    run_pipeline([[BIN_BCF, "index", "-f", "-c", str(final_vcf)]])

    mask_bed = bed_dir / "unknown.mask.bed"
    mask_bed.write_text("chr1\t0\t15\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="these samples are not present"):
        apply_sample_region_masks_to_resolved_vcf(tmpdir, {"unknown": mask_bed}, vcf_gz=final_vcf)
