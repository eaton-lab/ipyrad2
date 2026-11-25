"""Explicit BAM groups retain one sample, original data, and read provenance."""

import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

from ipyrad2.assembler.assemble import (
    _collect_named_bams, _normalize_bam_subsample_file, _prepare_grouped_input,
)
from ipyrad2.assembler.bam_grouping import merge_grouped_bams
from ipyrad2.assembler.beds import get_reference_sort_order
from ipyrad2.assembler.read_filters import prepare_filtered_analysis_bam
from ipyrad2.utils.exceptions import IPyradError

BIN = Path(sys.prefix) / "bin"
SAM = BIN / "samtools"
BCF = BIN / "bcftools"


def command(*args):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result.stdout


def make_bam(path, seq, *, rg="library", sample="original", n=3, start=31,
             ref_length=300, paired=False, low_mapq=False, loci=None):
    sam = path.with_suffix(".sam")
    header = f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr\tLN:{ref_length}\n"
    if rg is not None:
        header += f"@RG\tID:{rg}\tSM:{sample}\tLB:original-library\n"
    tags = "\tNM:i:0" + (f"\tRG:Z:{rg}" if rg is not None else "")
    records = []
    for locus, (start, seq) in enumerate(loci if loci is not None else [(start, seq)]):
        for i in range(n):
            q = 10 if low_mapq and i == 0 else 60
            name = f"read{locus}.{i}"
            if paired:
                records.append(f"{name}\t99\tchr\t{start}\t{q}\t50M\t=\t{start+30}\t80\t{seq[:50]}\t{'I'*50}{tags}\n")
                records.append(f"{name}\t147\tchr\t{start+30}\t{q}\t50M\t=\t{start}\t-80\t{seq[30:80]}\t{'I'*50}{tags}\n")
            else:
                records.append(f"{name}\t0\tchr\t{start}\t{q}\t{len(seq)}M\t*\t0\t0\t{seq}\t{'I'*len(seq)}{tags}\n")
    sam.write_text(header + "".join(records))
    command(SAM, "sort", "-o", path, sam)
    command(SAM, "index", "-c", path)
    return path


def test_explicit_groups_and_unintentional_collisions(tmp_path, monkeypatch):
    a, b, c = [tmp_path / x for x in ("a.bam", "b.bam", "c.bam")]
    # Fully specified explicit destination duplicates represent a group.
    result, renamed = _collect_named_bams([a, b, c], {a.name: "A", b.name: "A", c.name: "C"})
    assert result == {"A": (a, b), "C": c}
    assert renamed == [("A", "a.bam"), ("A", "b.bam"), ("C", "c.bam")]
    monkeypatch.setattr("ipyrad2.assembler.assemble.get_name_from_bam", lambda _: "A")
    with pytest.raises(IPyradError, match="explicitly assign every member"):
        _collect_named_bams([a, b], {a.name: "A"})
    monkeypatch.setattr("ipyrad2.assembler.assemble.get_names_from_bams",
                        lambda paths: {p: "A" for p in paths})
    with pytest.raises(IPyradError, match="Multiple input files"):
        _collect_named_bams([a, b], {})


def test_subsample_group_or_individual_and_ambiguous_alias(tmp_path, monkeypatch):
    a, b, c = [tmp_path / x for x in ("a.bam", "b.bam", "c.bam")]
    keep = tmp_path / "keep.tsv"
    renames = {a.name: "group", b.name: "group", c.name: "other"}
    keep.write_text("group\tpopulation\n")
    assert _normalize_bam_subsample_file(keep, [a,b,c], renames) == {a,b}
    keep.write_text("a.bam\n")
    assert _normalize_bam_subsample_file(keep, [a,b,c], renames) == {a}
    # A final group name must not absorb another BAM with a matching basename.
    keep.write_text("c.bam\n")
    with pytest.raises(IPyradError, match="ambiguous identifiers"):
        _normalize_bam_subsample_file(keep, [a,b,c], {a.name:"c.bam", b.name:"c.bam", c.name:"other"})


def test_grouping_preserves_reads_and_read_groups(tmp_path):
    seq = "ACGT" * 25
    a = make_bam(tmp_path/"a.bam", seq, rg="collision", sample="A", n=4, paired=True, low_mapq=True)
    b = make_bam(tmp_path/"b.bam", seq, rg="collision", sample="B", n=3, paired=True)
    c = make_bam(tmp_path/"c.bam", seq, rg="unique", sample="C", n=2, paired=True)
    d = make_bam(tmp_path/"d.bam", seq, rg=None, n=1, paired=True)
    members = (a,b,c,d)
    digests = [hashlib.sha256(p.read_bytes()).hexdigest() for p in members]
    merged = merge_grouped_bams(sname="merged", bam_paths=members, tmpdir=tmp_path, threads=1)
    assert Path(str(merged)+".csi").is_file()
    assert int(command(SAM,"view","-c",merged)) == 20
    header = command(SAM,"view","-H",merged)
    tags = [dict(x.split(":",1) for x in line.split("\t")[1:])
            for line in header.splitlines() if line.startswith("@RG\t")]
    assert len({d["ID"] for d in tags}) == len(tags)
    assert {d["SM"] for d in tags} == {"merged"}
    assert any(d["ID"]=="unique" and d["LB"]=="original-library" for d in tags)
    records = [l.split("\t") for l in command(SAM,"view",merged).splitlines()]
    assert len({x[0] for x in records}) == 10  # QNAMEs shared by mates, not libraries.
    assert len({next(t for t in x[11:] if t.startswith("RG:Z:")) for x in records}) == 4
    assert [hashlib.sha256(p.read_bytes()).hexdigest() for p in members] == digests
    assert not list(merged.parent.glob(".*"))  # member scratch cleaned
    result = prepare_filtered_analysis_bam(sname="merged", bam_file=merged, is_paired=True,
        tmpdir=tmp_path, min_map_q=30, max_tlen=1000, max_softclip=None, max_nm=None,
        min_aligned_len=None, threads=1)
    assert result.reads_before_filtering == 20
    assert result.reads_after_filtering == 18


def test_grouping_rejects_mixed_layout_or_reference(tmp_path):
    ref = tmp_path/"ref.fa"
    ref.write_text(">chr\n"+"A"*300+"\n")
    get_reference_sort_order(ref,tmp_path)
    a = make_bam(tmp_path/"a.bam", "A"*100)
    b = make_bam(tmp_path/"b.bam", "A"*100, paired=True)
    c = make_bam(tmp_path/"c.bam", "A"*100, ref_length=301)
    for other, msg in [(b,"mixes single-end and paired-end"),(c,"reference")]:
        with pytest.raises(IPyradError, match=msg):
            _prepare_grouped_input(sname="bad", bam_paths=(a,other), tmpdir=tmp_path,
                                   reference=ref, threads=1)
    assert not (tmp_path/"grouped_bams").exists()


@pytest.mark.parametrize("genotype_model", ["native", "sfs"])
def test_grouping_matches_manual_merge_at_assembly_outputs(tmp_path, genotype_model):
    # Each replicate has 3 reads: individually insufficient for 5x locus
    # discovery, but jointly sufficient. Compare to an independent manually
    # merged input, not a BAM produced by the implementation under test.
    rng = random.Random(42)
    sequence = "".join(rng.choice("ACGT") for _ in range(3000))
    ref = tmp_path/"ref.fa"
    ref.write_text(">chr\n"+sequence+"\n")
    loci = [(31 + 400*i, sequence[30+400*i:130+400*i]) for i in range(6)]
    alternate = [(start, piece[:49] + next(x for x in "ACGT" if x != piece[49]) + piece[50:])
                 for start,piece in loci]
    a = make_bam(tmp_path/"a.bam", "", rg="libA", sample="replicateA", loci=loci, ref_length=3000)
    b = make_bam(tmp_path/"b.bam", "", rg="libB", sample="replicateB", loci=loci, ref_length=3000)
    c = make_bam(tmp_path/"c.bam", "", rg="libC", sample="other", n=6, loci=alternate, ref_length=3000)
    # The fixture manually merges equivalent records with unique QNAMEs and
    # shared SM, independent of the grouping implementation.
    manual_sam = tmp_path/"manual.sam"
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr\tLN:3000",
             "@RG\tID:libA\tSM:specimen", "@RG\tID:libB\tSM:specimen"]
    for i, source in enumerate((a,b)):
        for line in command(SAM,"view",source).splitlines():
            lines.append(f"source{i}:"+line)
    manual_sam.write_text("\n".join(lines)+"\n")
    manual = tmp_path/"manual.bam"
    command(SAM,"sort","-o",manual,manual_sam)
    command(SAM,"index","-c",manual)
    rename = tmp_path/"rename.tsv"
    rename.write_text("a.bam specimen\nb.bam specimen\n")
    common = ["-r",ref,"-m","2","-a","2","-z","25","-g","10","-s","1","--min-locus-depth","5",
              "-c","2","-t","1","-qm","20","--evidence","loci",
              "--genotype-model",genotype_model,"--sfs-folds","2","--sfs-fit-loci","4"]
    def run(out, bams, extra):
        command(sys.executable,"-c","from ipyrad2.cli import main; main()","assemble","-d",*bams,"-o",out,*common,*extra)
    grouped, control = tmp_path/"grouped", tmp_path/"control"
    run(grouped,[a,b,c],["--rename",rename])
    run(control,[manual,c],[])
    # Compare output biological matrix and all variant genotypes/allele depths.
    gv = next(grouped.glob("*.vcf.gz")); cv = next(control.glob("*.vcf.gz"))
    assert command(BCF,"query","-l",gv)==command(BCF,"query","-l",cv)=="other\nspecimen\n"
    fmt = "%CHROM\t%POS\t%REF\t%ALT[\t%GT:%AD:%DP]\n"
    variants = command(BCF,"query","-f",fmt,gv)
    assert variants and variants==command(BCF,"query","-f",fmt,cv)
    assert not (grouped/"assembly_tmpdir").exists()
    report = json.loads(next(grouped.glob("*.stats.json")).read_text())
    assert report["run"]["inputs"]["RAD"]["specimen"] == [str(a),str(b)]
    control_report = json.loads(next(control.glob("*.stats.json")).read_text())
    for field in ("summary", "locus_filtering", "sample_masking", "alignment_summary",
                  "sample_summary", "read_filter_diagnostics", "locus_occupancy"):
        assert report[field] == control_report[field]
    assert report["summary"]["samples"] == 2
    assert report["summary"]["final_loci_written"] == 6
