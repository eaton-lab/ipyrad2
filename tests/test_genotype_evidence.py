from pathlib import Path
import subprocess

from ipyrad2.assembler.genotype_evidence import BIN_BCF, write_genotype_evidence


def test_genotype_evidence_retains_native_fields_without_editing_standard_vcf(tmp_path: Path):
    header = (
        '##fileformat=VCFv4.3\n'
        '##contig=<ID=chr1,length=200>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRAD\tWGS\n'
    )
    raw = tmp_path / "raw.vcf"
    raw.write_text(
        header
        + "chr1\t10\t.\tA\tC\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t0/1:8:4,4:30:30,0,30\t0/0:1:1,0:3:0,3,30\n"
        + "chr1\t20\t.\tA\tC,G\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4:20:60,20,30,30,0,20\t./.:0:0,0,0:.:.\n"
        + "chr1\t30\t.\tA\tG,C\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4:20:60,20,30,30,0,20\t./.:0:0,0,0:.:.\n"
    )
    final = tmp_path / "final.vcf"
    final.write_text(
        header
        + "chr1\t10\t.\tA\tC\t40\tPASS\t.\tGT\t0/1\t./.\n"
        + "chr1\t20\t.\tA\tC,G\t40\tPASS\t.\tGT\t1/2\t./.\n"
        + "chr1\t30\t.\tA\tC,G\t40\tPASS\t.\tGT\t1/2\t./.\n"
    )
    before = final.read_bytes()
    mask = tmp_path / "mask.bed"
    mask.write_text("chr1\t9\t10\n")
    out = tmp_path / "oak.evidence.bcf"
    write_genotype_evidence(raw_vcf=raw, final_vcf=final, out_bcf=out, hard_masks=({"RAD": mask},))
    assert final.read_bytes() == before
    assert out.with_suffix(".bcf.csi").is_file()
    observed = subprocess.check_output(
        [BIN_BCF, "query", "-r", "chr1", "-f", r"%POS[\t%GT:%DP:%AD:%GQ:%PL:%FT]\n", str(out)],
        text=True,
    )
    assert observed.splitlines() == [
        "10\t./.:.:.:.:.:HardMask\t./.:1:1,0:3:0,3,30:NoHardCall",
        "20\t1/2:6:0,2,4:20:60,20,30,30,0,20:PASS\t./.:0:.:.:.:NoRead",
    ]
