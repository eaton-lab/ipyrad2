"""Native-tool parity tests for diagnostic-only changes and GQ retention."""

import gzip
import json
from pathlib import Path
import subprocess

import pytest

from ipyrad2.assembler import read_filters as rf
from ipyrad2.assembler import variants


def sam(name, cigar="120M", flag=1, mq=60, mate="=", tlen=300, nm=1):
    import re

    length = sum(
        int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MIS=X"
    )
    return (
        "\t".join(
            map(
                str,
                [
                    name,
                    flag,
                    "chr1",
                    10,
                    mq,
                    cigar,
                    mate,
                    20,
                    tlen,
                    "A" * length,
                    "I" * length,
                ],
            )
        )
        + (f"\tNM:i:{nm}" if nm is not None else "")
        + "\n"
    )


def fixture_bam(tmp_path):
    lines = [
        sam("pass"),
        sam("multi", "30S70M", mq=10, nm=12),
        sam("clip", "25S100M"),
        sam("short", "94M"),
        sam("nm", nm=6),
        sam("missing", nm=None),
        sam("mate", mate="chr2"),
        sam("unmappedmate", flag=9),
        sam("long", tlen=-1100),
        sam("insertions", "70M40I"),
        sam("secondary", flag=257),
        sam("hardclip", "30H100M"),
        sam("boundary", "5S100M", mq=30, nm=5, tlen=-1000),
    ]
    path = tmp_path / "input.sam"
    path.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:2000\n"
        + "".join(lines)
    )
    bam = tmp_path / "input.bam"
    subprocess.run([rf.BIN_SAM, "view", "-b", "-o", str(bam), str(path)], check=True)
    return bam, lines


@pytest.mark.parametrize("paired", [False, True])
def test_diagnostics_match_native_filters_and_overlap_accounting(tmp_path, paired):
    bam, lines = fixture_bam(tmp_path)
    settings = dict(
        is_paired=paired,
        min_map_q=30,
        max_tlen=1000,
        max_softclip=5,
        max_nm=5,
        min_aligned_len=100,
    )
    diag = rf.summarize_read_filter_diagnostics(lines, **settings)
    expr = rf.build_mapped_read_filter_expr(
        **{k: v for k, v in settings.items() if k != "min_map_q"}
    )
    count = int(
        subprocess.check_output(
            [rf.BIN_SAM, "view", "-c", "-q", "30", "-e", expr, str(bam)]
        )
    )
    assert diag["predicted_retained_records"] == count
    assert sum(diag["first_failure"].values()) == diag["rejected_records"]
    assert sum(diag["independent_failures"].values()) > diag["rejected_records"]
    assert diag["records_without_nm"] == 1
    assert diag["input_flag_counts"]["secondary"] == 1  # reported, not newly excluded
    assert diag["independent_failures"]["mismatch_count"] == 3
    assert (
        "template_length" in diag["filter_order"]
        if paired
        else "template_length" not in diag["filter_order"]
    )
    # Native qlen-sclen includes insertions and excludes hard clipping.
    assert diag["independent_failures"]["aligned_query_length"] == 2
    # Each reason separately agrees with HTSlib, including overlapping failures
    # and missing NM; test first-failure counts against cumulative native gates.
    predicates = {
        "mapping_quality": "mapq >= 30",
        "pair_geometry": '((flag&4)==0) && ((flag&8)==0) && (rnext=="=" || rnext==rname)',
        "template_length": "tlen>=-1000 && tlen<=1000",
        "softclip": "sclen<=5",
        "mismatch_count": "[NM]<=5",
        "aligned_query_length": "qlen-sclen>=100",
    }
    previous = len(lines)
    cumulative = []
    for key in diag["filter_order"]:
        passing = int(
            subprocess.check_output(
                [rf.BIN_SAM, "view", "-c", "-e", predicates[key], str(bam)]
            )
        )
        assert diag["independent_failures"][key] == len(lines) - passing
        cumulative.append("(" + predicates[key] + ")")
        passing = int(
            subprocess.check_output(
                [rf.BIN_SAM, "view", "-c", "-e", " && ".join(cumulative), str(bam)]
            )
        )
        assert diag["first_failure"][key] == previous - passing
        previous = passing


@pytest.mark.parametrize("paired", [False, True])
def test_prepared_bam_unchanged_and_diagnostics_persist(tmp_path, paired):
    bam, _ = fixture_bam(tmp_path)
    settings = dict(
        is_paired=paired,
        min_map_q=30,
        max_tlen=1000,
        max_softclip=5,
        max_nm=5,
        min_aligned_len=100,
    )
    result = rf.prepare_filtered_analysis_bam(
        sname="s", bam_file=bam, tmpdir=tmp_path / "TMP", threads=1, **settings
    )
    expr = rf.build_mapped_read_filter_expr(
        **{k: v for k, v in settings.items() if k != "min_map_q"}
    )
    expected = subprocess.check_output(
        [rf.BIN_SAM, "view", "-q", "30", "-e", expr, str(bam)]
    ).splitlines()
    observed = subprocess.check_output(
        [rf.BIN_SAM, "view", str(result.bam_path)]
    ).splitlines()
    assert sorted(observed) == sorted(
        expected
    )  # fixture has no overlapping reciprocal mates
    diag = json.loads(result.bam_path.with_suffix(".filters.json").read_text())
    assert diag == result.filter_diagnostics
    assert diag["counts_reconcile"]
    assert result.reads_before_filtering == 13
    assert result.reads_after_filtering == len(expected)


def test_disabled_filters_missing_nm_empty_and_unmapped():
    settings = dict(
        is_paired=False,
        min_map_q=0,
        max_tlen=None,
        max_softclip=None,
        max_nm=None,
        min_aligned_len=None,
    )
    diag = rf.summarize_read_filter_diagnostics(
        [sam("missing", nm=None), sam("unmapped", flag=4)], **settings
    )
    assert diag["predicted_retained_records"] == 2
    assert diag["filter_order"] == ["mapping_quality"]
    assert rf.summarize_read_filter_diagnostics([], **settings)["input_records"] == 0


def test_gq_retention_preserves_calls_masks_and_final_output(tmp_path, monkeypatch):
    directory = tmp_path / "vcfs"
    directory.mkdir()
    raw = tmp_path / "raw.vcf"
    raw.write_text("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=2000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
chr1\t10\t.\tA\tC\t60\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:20:40:10,10:40,0,40\t1/1:20:5:0,20:40,20,0
chr1\t20\t.\tA\tG\t60\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:3:30:1,2:30,0,30\t0/0:25:60:25,0:0,60,90
chr1\t30\t.\tA\tT\t10\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:20:40:10,10:40,0,40\t0/0:20:40:20,0:0,40,90
chr1\t40\t.\tA\tC\t60\tPASS\t.\tGT:DP:GQ:AD:PL\t./.:.:.:.,.:.,.,.\t0/1:20:30:10,10:30,0,30
""")
    subprocess.run(
        [
            variants.BIN_BCF,
            "view",
            "-Oz",
            "-o",
            str(directory / "loci.raw.vcf.gz"),
            str(raw),
        ],
        check=True,
    )
    filtered = variants.get_filtered_vcf(tmp_path, 5, 13, 30, 1)
    expected = subprocess.check_output(
        [variants.BIN_BCF, "query", "-f", "%POS[\t%GT:%GQ]\n", str(filtered)]
    ).decode()
    assert expected.splitlines() == [
        "10\t0/1:40\t./.:5",
        "20\t./.:30\t0/0:60",
        "30\t./.:40\t./.:40",
        "40\t./.:.\t0/1:30",
    ]
    # Run the previous pipeline against the same raw input, changing only its
    # output-tag removal, to prove all non-GQ data and genotype masks agree.
    original_run = variants.run_pipeline

    def legacy(cmds, *args, **kwargs):
        for cmd in cmds:
            if "annotate" in cmd and "-x" in cmd:
                i = cmd.index("-x") + 1
                cmd[i] = "FORMAT/GQ," + cmd[i]
        return original_run(cmds, *args, **kwargs)

    monkeypatch.setattr(variants, "run_pipeline", legacy)
    baseline = variants.get_filtered_vcf(tmp_path, 5, 13, 30, 1)
    legacy_rows = subprocess.check_output(
        [variants.BIN_BCF, "view", "-H", str(baseline)]
    )
    # The final writer retains FORMAT fields; exercise it with real tools.
    monkeypatch.setattr(variants, "run_pipeline", original_run)
    filtered = variants.get_filtered_vcf(tmp_path, 5, 13, 30, 1)
    without_gq = subprocess.check_output(
        [variants.BIN_BCF, "annotate", "-x", "FORMAT/GQ", str(filtered)]
    )
    assert (
        b"\n".join(x for x in without_gq.splitlines() if not x.startswith(b"#")) + b"\n"
        == legacy_rows
    )
    bed = tmp_path / "loci.bed"
    bed.write_text("chr1\t0\t100\n")
    final = tmp_path / "final.vcf.gz"
    variants._write_unmasked_final_vcf(
        loci_bed=bed, in_vcf_gz=filtered, out_vcf_gz=final, threads=1
    )
    assert (
        "##FORMAT=<ID=GQ,"
        in subprocess.check_output(
            [variants.BIN_BCF, "view", "-h", str(final)]
        ).decode()
    )
    assert subprocess.check_output(
        [variants.BIN_BCF, "query", "-f", "%POS[\t%GT:%GQ]\n", str(final)]
    ).decode().splitlines() == [expected.splitlines()[i] for i in (0, 1, 3)]

    # Exercise final per-sample masking too: GQ must survive when GT is masked.
    import shutil

    shutil.copyfile(filtered, directory / "variants.resolved.vcf.gz")
    subprocess.run(
        [
            variants.BIN_BCF,
            "index",
            "-f",
            "-c",
            str(directory / "variants.resolved.vcf.gz"),
        ],
        check=True,
    )
    output = tmp_path / "OUT"
    output.mkdir()
    (output / "assembly.bed").write_text(bed.read_text())
    (tmp_path / "beds").mkdir()
    (tmp_path / "REF_info.txt").write_text("chr1\t2000\n")
    mask = tmp_path / "sample-mask.bed"
    mask.write_text("chr1\t9\t10\n")
    masked = variants.write_vcf(
        "assembly", output, tmp_path, 1, sample_masks={"s1": mask}, cores=1
    )
    masked_rows = (
        subprocess.check_output(
            [variants.BIN_BCF, "query", "-f", "%POS[\t%GT:%GQ]\n", str(masked)]
        )
        .decode()
        .splitlines()
    )
    assert masked_rows[0] == "10\t./.:40\t./.:5"
    assert masked_rows[-1] == "40\t./.:.\t0/1:30"
