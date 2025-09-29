#!/usr/bin/env python

"""INDIVIDUAL AND/OR JOINT VARIANT CALLS

Given a set of bam files this step makes variant calls and writes the
output to a VCF, and identified the loci windows with sufficient
coverage across samples and writes to a BED file.

SUMMARY
-------
1. [bedtools] Get locus BEDs above depth = 1 for all samples.
2. [bedtools] Get merged locus BEDs across >=4 samples (including REF)
3. [bcftools] Get variants for each ind or pop in locus BEDs.
4. [bcftools] Get filtered variants using per-sample (DP/GQ) and per-site (QUAL)
5. [bcftools] Get Norm/Merge variants to separate indels from snps.

PARAMS TO EXPOSE TO USER
------------------------
min_locus_cov:
    Minimum across-sample coverage for delimiting loci.
min_locus_len:
    Minimum length of a delimited locus. Note this is applied both
    during locus delim, and after trimming edges.
min_call_depth:
    Minimum within-sample coverage below which samples are masked from
    making variant calls.
merge_distance:
    Merge locus beds that are within this distance from each other.

PER-SAMPLE STATS TO RECORD
--------------------------
- nloci_with_nonzero_mapping: int
- median_depth_per_locus_with_nonzero_mapping: float
- median_depth_per_locus_total: float

PROJECT STATS TO RECORD
-----------------------
locus_beds:
    number of delimited locus beds.
locus_beds_mean_sample_cov:
    mean sample coverage per locus
locus_beds_stdev_sample_cov:
    mean sample coverage per locus
nvariants_raw
    number of variants in raw vcf
nvariants_filtered
    number of variants in filtered vcf
nsnps_in_nvariants_filtered
    number of snp variants in filtered vcf
nindels_in_nvariants_filtered
    number of indel variants in filtered vcf
"""

from typing import List
import os
import sys
import shlex
from pathlib import Path
import subprocess as sp
from ..utils.parallel import run_pipeline, run_with_pool

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")

# ==========================================================================


def get_locus_and_snp_stats_in_loci_bed(outdir: Path, threads: int = 4):
    """Return dict with stats"""
    # file paths
    loci_bed = outdir / "beds" / "loci.bed"
    raw_vcf = outdir / "vcfs" / "loci.raw.vcf.gz"
    vcf = outdir / "vcfs" / "variants.resolved.vcf.gz"

    # get the number of loci beds
    cmd1 = ["wc", "-l", str(loci_bed)]
    p1 = sp.Popen(cmd1, stdout=sp.PIPE)

    # get the number of snps
    cmd2 = [
        BIN_BCF, "view",
        "-f", "PASS",
        "-H",
        "-v", "snps",
        "--threads", str(threads),
        str(vcf),
    ]
    cmd3 = ["wc", "-l"]
    p2 = sp.Popen(cmd2, stdout=sp.PIPE, stderr=sp.DEVNULL)
    p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=sp.PIPE, stderr=sp.PIPE)

    # get the number of indels
    cmd4 = [
        BIN_BCF, "view",
        "-f", "PASS",
        "-H",
        "-v", "indels",
        "--threads", str(threads),
        str(vcf),
    ]
    cmd5 = ["wc", "-l"]
    p4 = sp.Popen(cmd4, stdout=sp.PIPE, stderr=sp.DEVNULL)
    p5 = sp.Popen(cmd5, stdin=p4.stdout, stdout=sp.PIPE, stderr=sp.PIPE)

    # get the number of raw variants
    cmd6 = [
        BIN_BCF, "view",
        "-H",
        "-v", "snps,indels",
        "--threads", str(threads),
        str(raw_vcf),
    ]
    cmd7 = ["wc", "-l"]
    p6 = sp.Popen(cmd6, stdout=sp.PIPE, stderr=sp.DEVNULL)
    p7 = sp.Popen(cmd7, stdin=p6.stdout, stdout=sp.PIPE, stderr=sp.PIPE)

    # parse result nloci
    out, _ = p1.communicate()
    nloci = int(out.decode().strip().split()[0])
    if not nloci:
        raise RuntimeError("No locus beds passed filtering. Considering lowering the min sample coverage.")

    # parse nsnps
    out, err = p3.communicate()
    if p2.returncode:
        raise RuntimeError(f"Error in {cmd1}: {err.decode()}")
    nsnps = int(out.decode())

    # parse indels
    out, err = p5.communicate()
    if p2.returncode:
        raise RuntimeError(f"Error in {cmd1}: {err.decode()}")
    nindels = int(out.decode())

    # parse raw variants
    out, err = p7.communicate()
    if p2.returncode:
        raise RuntimeError(f"Error in {cmd1}: {err.decode()}")
    nvariants_raw = int(out.decode())

    return {"nloci": nloci, "nvariants": nsnps + nindels, "nsnps": nsnps, "nindels": nindels, "nvariants_raw": nvariants_raw}


def get_chunked_loci_beds(outdir: Path, nchunks: int) -> List[Path]:
    """Return a list of Paths from breaking loci.bed into chunks.
    """
    loci_bed = outdir / "beds" / "loci.bed"
    lines = loci_bed.read_text().split("\n")
    q, r = divmod(len(lines), nchunks)

    paths = []
    i = 0
    for k in range(nchunks):
        chunk_bed = outdir / "beds" / f"chunk-{i}.bed"
        size = q + (1 if k < r else 0)
        chunk = lines[i: i+size]
        with open(chunk_bed, 'w') as out:
            out.write("\n".join(chunk))
        paths.append(chunk_bed)
        i += size
    return paths


def get_group_called_variants_in_vcf_chunks(outdir: Path, reference: Path, bam_files: List[Path], locus_chunk: Path, threads: int):
    """Make variant calls for all samples using -G (groups).

    >>> $ bcftools mpileup \
    >>>       -f REF -q 20 -Q 20 -d 10000 -R loci.bed \
    >>>       -a FMT/DP,FMT/AD -Ou S1.bam S2.bam ... \
    >>>   | bcftools call -m -a GQ -G GROUPS.tsv -W -Oz -o out.vcf.gz
    """
    # file paths
    vcf_dir = outdir / "vcfs"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    out_vcf_gz = vcf_dir / locus_chunk.with_suffix(".vcf.gz").name

    # divide threads
    threads_mpileup = max(1, threads // 2)
    threads_call = max(1, threads - threads_mpileup)

    # get genotype likelihoods at all sites in Region with decent mapping.
    cmd1 = [
        BIN_BCF, "mpileup",
        "-f", str(reference),
        "-q", str(20),
        "-Q", str(20),
        "-d", str(10_000),
        "-a", "FMT/DP,FMT/AD",
        "-R", str(locus_chunk),
        "--threads", str(threads_mpileup),
        "-Ou",
    ] + [str(i) for i in bam_files]

    # call variants with GQ scores, write to tmp VCF and index it.
    cmd2 = [
        BIN_BCF, "call",
        "-m",
        "-a", "GQ",
        # "-W",  # index after concatenating
        "-G", "-",
        # "-G", str(groups_file),
        "-Oz",
        "-o", str(out_vcf_gz),
        "--threads", str(threads_call),
    ]

    # run cmds 1 and 2 and report errors
    e1 = open(log_dir / "mpileup.err", "wb")
    e2 = open(log_dir / "call.err", "wb")
    try:
        p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=e1)
        p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.DEVNULL, stderr=e2)
        if p1.stdout:
            p1.stdout.close()  # allow SIGPIPE to p1 if p2 dies
        rc2 = p2.wait()
        rc1 = p1.wait()

    finally:
        for fh in (e1, e2):
            try:
                fh.close()
            except Exception:
                pass

    if any(rc != 0 for rc in (rc1, rc2)):
        cmds = "\n".join(shlex.join(c) for c in (cmd1, cmd2))
        raise RuntimeError(
            f"bcftools pipeline failed: mpileup={rc1}, call={rc2}\n"
            f"Commands were:\n{cmds}\nLogs in: {log_dir}"
        )
    return out_vcf_gz


def get_concat_chunk_vcfs(outdir: Path, threads: int):
    """Concatenate filtered vcf chunks back into one large vcf.

    The vcfs should already be indexed. This will re-sort just in case.
    # bcftools concat V1 V2 V3 ... -Oz -o loci.vcf.gz --threads 8 -W
    """
    vcf_dir = outdir / "vcfs"
    out_vcf_gz = vcf_dir / "loci.raw.vcf.gz"
    chunk_vcfs = vcf_dir.glob("chunk-*.vcf.gz")
    sorted_vcfs = sorted(chunk_vcfs, key=lambda x: int(x.name.split(".")[0].split("-")[-1]))

    # build command
    cmd = [
        BIN_BCF, "concat",
        "--threads", str(threads),
        "-Oz", "-o", str(out_vcf_gz),
        "-W",
    ] + [str(i) for i in sorted_vcfs]

    # run command
    proc = sp.Popen(cmd, stdout=sp.DEVNULL, stderr=sp.STDOUT)
    _, err = proc.communicate()
    if proc.returncode:
        raise RuntimeError(f"Error in {cmd}: {err.decode()}")

    # clean up tmp chunk files
    for chunk in sorted_vcfs:
        if chunk.exists():
            chunk.unlink()
    return out_vcf_gz


def get_filtered_vcf(outdir: Path, min_read_depth: int, min_gq: int, min_qual: int, threads: int) -> Path:
    """Apply filtering to raw genotype calls by depth and quality

    $ bcftools +setGT VCF -- -t q -n . -i "FMT/DP<X | FMT/GQ<Y"
    $ bcftools +setGT VCF -- -t q -n . -i "QUAL<Z"
    $ bcftools +fill-tags VCF -- -t "AC,AN,AF,MAF,F_MISSING"

    TODO: record how many alleles are masked by [1] and how
    many sites are masked by [2].
    """
    in_vcf_gz = outdir / "vcfs" / "loci.raw.vcf.gz"
    out_vcf_gz = outdir / "vcfs" / "loci.filtered.vcf.gz"
    out_vcf_tmp = out_vcf_gz.with_suffix(out_vcf_gz.suffix + ".tmp")
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    dp_min: int = min_read_depth                   # DP<X (pser sample DP)
    gq_min: int = min_gq                           # GQ<X (per-sample GQ)
    # indel_snp_mask: int = 0                      # filter SNPs within N bases on an indel
    qual_min: int = min_qual                       # QUAL (across samples)
    threads = max(1, int(threads / 2))      # assign threads among piped jobs

    # filter per-sample genotypes by min depth and geno quality
    expr_gt_mask = f"FMT/DP<{dp_min} | FMT/GQ<{gq_min}"
    cmd1 = [
        BIN_BCF, "+setGT", str(in_vcf_gz),
        "--",
        "-t", "q",
        "-n", ".",
        "-i", expr_gt_mask,
    ]

    # filter site by QUAL across all samples
    # cmd2 = [BIN_BCF, "+setGT", "-",            "--", "-t", "q", "-n", ".", "-i", expr_site_mask]
    expr_site_mask = f"QUAL<{qual_min}"
    cmd2 = [
        BIN_BCF, "filter",
        "-S", ".",
        "-s", "lowQual",
        # "-g", "1",               # filter SNPs within N bases on an indel
        "-e", expr_site_mask,
        "--threads", str(threads),
        "-Ou", "-",
    ]

    # compute new tags
    cmd3 = [
        BIN_BCF, "+fill-tags", "-",
        "--",
        "-t", "AC,AN,AF,MAF,F_MISSING",
    ]

    # filter sites that are no longer variants
    # NO, we don't need to do this here. This is no the final VCF, this
    # is one used to create consensus alignments. We want to retain the
    # info that these sites are masked as Ns if all or some samples.
    # Monomorphic sites will be removed from the VCF in the last assembly step.

    # clean up tags to keep only minimal
    remove_tags = "FORMAT/PL,FORMAT/GQ,INFO/RPBZ,INFO/SCBZ,INFO/MQBZ,INFO/BQBZ,INFO/MQSBZ,INFO/DP4,INFO/VDB,INFO/MQ0F,INFO/SGB"
    cmd4 = [BIN_BCF, "annotate",
        "-x", remove_tags,
        "-Oz", "-o", str(out_vcf_tmp),
        "--threads", str(threads),
        "-",
    ]

    e1 = open(log_dir / f"{out_vcf_gz.stem}.setgt.err", "wb")
    p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=e1)

    e2 = open(log_dir / f"{out_vcf_gz.stem}.filter.err", "wb")
    p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, stderr=e2)
    if p1.stdout:
        p1.stdout.close()

    e3 = open(log_dir / f"{out_vcf_gz.stem}.filltags.err", "wb")
    p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=sp.PIPE, stderr=e3)
    if p2.stdout:
        p2.stdout.close()

    e4 = open(log_dir / f"{out_vcf_gz.stem}.annotate.err", "wb")
    p4 = sp.Popen(cmd4, stdin=p3.stdout, stdout=sp.PIPE, stderr=e4)
    if p3.stdout:
        p3.stdout.close()

    # ---- wait in downstream->upstream order ----
    rc4 = p4.wait()
    rc3 = p3.wait()
    rc2 = p2.wait()
    rc1 = p1.wait()

    # ---- close logs ----
    for fh in (e1, e2, e3, e4):
        try:
            fh.close()
        except Exception:
            pass

    # ---- error reporting ----
    bad = [(name, rc, cmd) for name, rc, cmd in (
        ("setGT/genotype", rc1, cmd1),
        ("filter",         rc2, cmd2),
        ("fill-tags",      rc3, cmd3),
        ("annotate",       rc4, cmd4),
    ) if rc != 0]

    if bad:
        import shlex
        detail = "\n".join(f"{n} rc={rc} :: $ {shlex.join(c)}" for n, rc, c in bad)
        raise RuntimeError(f"bcftools pipeline failed:\n{detail}\nLogs in: {log_dir}")

    # ---- finalize & index ----
    os.replace(out_vcf_tmp, out_vcf_gz)
    sp.run([BIN_BCF, "index", "-f", "-c", str(out_vcf_gz)], check=True)

    # clean up by removing raw SNPs file
    # if in_vcf_gz.exists():
    #     in_vcf_gz.unlink()
    # if in_vcf_gz.with_suffix(in_vcf_gz.suffix + ".csi").exists():
    #     in_vcf_gz.with_suffix(in_vcf_gz.suffix + ".csi").unlink()
    return out_vcf_gz


def get_vcf_with_indels_resolved(outdir: Path, reference: Path, threads: int) -> Path:
    """Resolve overlapping snps and indels. Keep indel type when overlapping.

    If no indels are present then it just renames the vcf.

    Steps:
      1) norm -m -both
      2) split to snps or indels
      3) make indel.regions.bed from REF/ALT lengths
      4) drop conflicting SNPs
      5) concat + sort
      6) collapse biallelic back to multiallelic + sort
    """
    ref_fa = Path(reference)
    in_vcf_gz = outdir / "vcfs" / "loci.filtered.vcf.gz"
    out_vcf_gz = outdir / "vcfs" / "variants.resolved.vcf.gz"
    vcf_dir = outdir / "vcfs"
    bed_dir = outdir / "beds"

    # ------------------------------------------------------------
    # 1) Normalize & decompose to primitives
    cmd1 = [
        BIN_BCF, "norm",
        "-f", str(ref_fa),
        "-m", "-both",
        "--threads", str(threads),
        "-W",
        "-Oz", "-o", str(vcf_dir / "norm.vcf.gz"),
        str(in_vcf_gz),
    ]
    run_pipeline(cmd1)

    # ------------------------------------------------------------
    # 2) Split into SNPs vs. INDELs
    cmd1 = [
        BIN_BCF, "view",
        "-v", "snps",
        "-Oz", "-o", str(vcf_dir / "snps.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF, "view",
        "-v", "indels",
        "-Oz", "-o", str(vcf_dir / "indels.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    run_pipeline([cmd1, cmd2], )

    # -----------------------------------------------------------
    # 3) Build an indel-affected BED (0-based, half-open)
    # bcftools query | awk (length-based) | sort | bedtools merge > bed
    awk_prog = (
        r'BEGIN{OFS="\t"}'
        r'{chrom=$1; pos0=$2; ref=$4; n=split($5,alts,",");'
        r' for(i=1;i<=n;i++){alt=alts[i];'
        r'  if(length(ref)>length(alt)){print chrom, pos0, pos0+length(ref);} '
        r'  else if(length(alt)>length(ref)){print chrom, pos0, pos0+1;} '
        r' }}'
    )
    cmd1 = [
        BIN_BCF, "query",
        "-f", r"%CHROM\t%POS0\t%POS\t%REF\t%ALT\n",
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = ["awk", awk_prog]
    cmd3 = ["sort", "-k1,1", "-k2,2n", "-T", str(vcf_dir)]
    cmd4 = [BIN_BED, "merge", "-i", "-"]

    # run pipeline
    run_pipeline()



    # ----------------------------------------------------------
    # 4) drop the conflicting SNPs
    cmd1 = [
        BIN_BCF, "view",
        "-T", f"^{str(bed_dir / 'indel.regions.bed')}",
        "-Oz", "-o", str(vcf_dir / "snps.clean.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "snps.vcf.gz"),
    ]
    sp.run(cmd1, check=True)

    # ----------------------------------------------------------
    # 6) keep REF rows but avoid duplicate POS lines
    # Build variant.pos.bed = 1-bp windows of SNP positions + indel anchors
    # We stream both queries into a single sort|uniq pipeline.
    cmd1 = ["sort", "-k1,1", "-k2,2n", "-T", str(vcf_dir)]
    cmd2 = ["uniq"]

    # get indel positions to rm variants from
    variant_pos_bed = bed_dir / "variant.pos.bed"
    with open(variant_pos_bed, "w") as o2:
        # open processors waiting on stdin
        sort_p = sp.Popen(cmd1, stdin=sp.PIPE, stdout=sp.PIPE, text=True)
        uniq_p = sp.Popen(cmd2, stdin=sort_p.stdout, stdout=o2, text=True)

        # iterate over each file to feed it in
        vfiles = [vcf_dir / "snps.clean.vcf.gz", vcf_dir / "indels.vcf.gz"]
        for vcf_file in vfiles:
            # Write both files into sort_p stdin
            cmd_ = [BIN_BCF, "query", "-f", r"%CHROM\t%POS0\t%POS\n", str(vcf_file)]
            q2 = sp.run(cmd_, check=True, capture_output=True, text=True)
            sort_p.stdin.write(q2.stdout)

        # wait for sort and uniq processes to finish
        sort_p.stdin.close()
        rc_sort = sort_p.wait()
        rc_uniq = uniq_p.wait()

        # check for errors
        if rc_sort != 0 or rc_uniq != 0:
            raise RuntimeError("Failed to create variant.pos.bed")

    # ----------------------------------------------------------
    # 7) Recombine (refs + clean SNPs + indels) and sort
    cmd1 = [
        BIN_BCF, "concat",
        "-a",
        "-Oz", "-o", str(vcf_dir / "combined.vcf.gz"),
        "--threads", str(threads),
        # str(vcf_dir / "refs.clean.vcf.gz"),
        str(vcf_dir / "snps.clean.vcf.gz"),
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF, "sort",
        "-Oz", "-o", str(vcf_dir / "combined.sorted.vcf.gz"),
        "-T", str(vcf_dir),
        "-W",
        str(vcf_dir / "combined.vcf.gz"),
    ]
    sp.run(cmd1, check=True, capture_output=True)
    sp.run(cmd2, check=True, capture_output=True)

    # 8) Collapse biallelic records at same POS back to multi-allelic; sort & index
    cmd1 = [
        BIN_BCF, "norm",
        "-m", "+both",
        "-Oz", "-o", str(vcf_dir / "combined.multi.vcf.gz"),
        "--threads", str(threads),
        str(vcf_dir / "combined.sorted.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF, "sort",
        "-Oz", "-o", str(out_vcf_gz),
        "-W",
        str(vcf_dir / "combined.multi.vcf.gz")
    ]
    sp.run(cmd1, check=True, capture_output=True)
    sp.run(cmd2, check=True, capture_output=True)

    # clean up
    # for path in vcf_dir.glob("*.vcf.gz"):
    #     if path.name != out_vcf_gz.name:
    #         if path.exists():
    #             print(f'removing {path}')
    #             path.unlink()
    #         ipath = path.with_suffix(path.suffix + ".csi")
    #         if ipath.exists():
    #             print(f'removing {ipath}')
    #             ipath.unlink()
    return out_vcf_gz


def old_get_vcf_with_indels_resolved(outdir: Path, reference: Path, threads: int) -> Path:
    """Resolve overlapping snps and indels. Keep indel type when overlapping.

    If no indels are present then it just renames the vcf.

    Steps:
      1) norm -m -both
      2) split to snps or indels
      3) make indel.regions.bed from REF/ALT lengths
      4) drop conflicting SNPs
      5) concat + sort
      6) collapse biallelic back to multiallelic + sort
    """
    ref_fa = Path(reference)
    in_vcf_gz = outdir / "vcfs" / "loci.filtered.vcf.gz"
    out_vcf_gz = outdir / "vcfs" / "variants.resolved.vcf.gz"
    vcf_dir = outdir / "vcfs"
    bed_dir = outdir / "beds"
    log_dir = outdir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1) Normalize & decompose to primitives
    cmd1 = [
        BIN_BCF, "norm",
        "-f", str(ref_fa),
        "-m", "-both",
        "--threads", str(threads),
        "-W",
        "-Oz", "-o", str(vcf_dir / "norm.vcf.gz"),
        str(in_vcf_gz),
    ]
    sp.run(cmd1, check=True, capture_output=True)

    # ------------------------------------------------------------
    # 2) Split into SNPs, INDELs, and REF rows
    cmd1 = [
        BIN_BCF, "view",
        "-v", "snps",
        "-Oz", "-o", str(vcf_dir / "snps.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    sp.run(cmd1, check=True)
    cmd2 = [
        BIN_BCF, "view",
        "-v", "indels",
        "-Oz", "-o", str(vcf_dir / "indels.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "norm.vcf.gz"),
    ]
    sp.run(cmd2, check=True)

    # -----------------------------------------------------------
    # 3) Build an indel-affected BED (0-based, half-open)
    # bcftools query | awk (length-based) | sort | bedtools merge > bed
    awk_prog = (
        r'BEGIN{OFS="\t"}'
        r'{chrom=$1; pos0=$2; ref=$4; n=split($5,alts,",");'
        r' for(i=1;i<=n;i++){alt=alts[i];'
        r'  if(length(ref)>length(alt)){print chrom, pos0, pos0+length(ref);} '
        r'  else if(length(alt)>length(ref)){print chrom, pos0, pos0+1;} '
        r' }}'
    )
    cmd1 = [
        BIN_BCF, "query",
        "-f", r"%CHROM\t%POS0\t%POS\t%REF\t%ALT\n",
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = ["awk", awk_prog]
    cmd3 = ["sort", "-k1,1", "-k2,2n", "-T", str(vcf_dir)]
    cmd4 = [BIN_BED, "merge", "-i", "-"]

    # run pipeline
    e1 = open(log_dir / "query.err", "wb")
    p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=e1)

    e2 = open(log_dir / "awk.err", "wb")
    p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, stderr=e2)
    if p1.stdout:
        p1.stdout.close()

    e3 = open(log_dir / "sort.err", "wb")
    p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=sp.PIPE, stderr=e3)
    if p2.stdout:
        p2.stdout.close()

    o4 = open(bed_dir / "indel.regions.bed", "wb")
    e4 = open(log_dir / "merge.err", "wb")
    p4 = sp.Popen(cmd4, stdin=p3.stdout, stdout=o4, stderr=e4)
    if p3.stdout:
        p3.stdout.close()

    # wait for jobs to finish
    r4 = p4.wait()
    r3 = p3.wait()
    r2 = p2.wait()
    r1 = p1.wait()
    for fh in (e1, e2, e3, e4, 4):
        try:
            fh.close()
        except Exception:
            pass
    # collect errors
    if any(rc != 0 for rc in (r1, r2, r3, r4)):
        raise RuntimeError(f"Error. See logs in {log_dir}")

    # ----------------------------------------------------------
    # 4) drop the conflicting SNPs
    cmd1 = [
        BIN_BCF, "view",
        "-T", f"^{str(bed_dir / 'indel.regions.bed')}",
        "-Oz", "-o", str(vcf_dir / "snps.clean.vcf.gz"),
        "--threads", str(threads),
        "-W",
        str(vcf_dir / "snps.vcf.gz"),
    ]
    sp.run(cmd1, check=True)

    # ----------------------------------------------------------
    # 6) keep REF rows but avoid duplicate POS lines
    # Build variant.pos.bed = 1-bp windows of SNP positions + indel anchors
    # We stream both queries into a single sort|uniq pipeline.
    cmd1 = ["sort", "-k1,1", "-k2,2n", "-T", str(vcf_dir)]
    cmd2 = ["uniq"]

    # get indel positions to rm variants from
    variant_pos_bed = bed_dir / "variant.pos.bed"
    with open(variant_pos_bed, "w") as o2:
        # open processors waiting on stdin
        sort_p = sp.Popen(cmd1, stdin=sp.PIPE, stdout=sp.PIPE, text=True)
        uniq_p = sp.Popen(cmd2, stdin=sort_p.stdout, stdout=o2, text=True)

        # iterate over each file to feed it in
        vfiles = [vcf_dir / "snps.clean.vcf.gz", vcf_dir / "indels.vcf.gz"]
        for vcf_file in vfiles:
            # Write both files into sort_p stdin
            cmd_ = [BIN_BCF, "query", "-f", r"%CHROM\t%POS0\t%POS\n", str(vcf_file)]
            q2 = sp.run(cmd_, check=True, capture_output=True, text=True)
            sort_p.stdin.write(q2.stdout)

        # wait for sort and uniq processes to finish
        sort_p.stdin.close()
        rc_sort = sort_p.wait()
        rc_uniq = uniq_p.wait()

        # check for errors
        if rc_sort != 0 or rc_uniq != 0:
            raise RuntimeError("Failed to create variant.pos.bed")

    # ----------------------------------------------------------
    # 7) Recombine (refs + clean SNPs + indels) and sort
    cmd1 = [
        BIN_BCF, "concat",
        "-a",
        "-Oz", "-o", str(vcf_dir / "combined.vcf.gz"),
        "--threads", str(threads),
        # str(vcf_dir / "refs.clean.vcf.gz"),
        str(vcf_dir / "snps.clean.vcf.gz"),
        str(vcf_dir / "indels.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF, "sort",
        "-Oz", "-o", str(vcf_dir / "combined.sorted.vcf.gz"),
        "-T", str(vcf_dir),
        "-W",
        str(vcf_dir / "combined.vcf.gz"),
    ]
    sp.run(cmd1, check=True, capture_output=True)
    sp.run(cmd2, check=True, capture_output=True)

    # 8) Collapse biallelic records at same POS back to multi-allelic; sort & index
    cmd1 = [
        BIN_BCF, "norm",
        "-m", "+both",
        "-Oz", "-o", str(vcf_dir / "combined.multi.vcf.gz"),
        "--threads", str(threads),
        str(vcf_dir / "combined.sorted.vcf.gz"),
    ]
    cmd2 = [
        BIN_BCF, "sort",
        "-Oz", "-o", str(out_vcf_gz),
        "-W",
        str(vcf_dir / "combined.multi.vcf.gz")
    ]
    sp.run(cmd1, check=True, capture_output=True)
    sp.run(cmd2, check=True, capture_output=True)

    # clean up
    # for path in vcf_dir.glob("*.vcf.gz"):
    #     if path.name != out_vcf_gz.name:
    #         if path.exists():
    #             print(f'removing {path}')
    #             path.unlink()
    #         ipath = path.with_suffix(path.suffix + ".csi")
    #         if ipath.exists():
    #             print(f'removing {ipath}')
    #             ipath.unlink()
    return out_vcf_gz


if __name__ == "__main__":
    # main()
    import ipyrad as ip

    DIR = "Ama"
    NAME = "COL"
    DIR = "Ped"
    NAME = "BIG2"

    # load assembly and select a sample
    data = ip.load_json(f"../../{DIR}/{NAME}.json")
    data.stepdir = data.json_file.parent / f"{NAME}_clusters_within"
    logs = data.stepdir / "logs"
    logs.mkdir(exist_ok=True)

    # get_vcf_with_indels_resolved(data, 10)

    print(get_locus_and_snp_stats_in_loci_bed(data, 10))

    # get_chunked_loci_beds()
    # raw_vcf = get_concat_chunk_vcfs(data, vchunks, 8)
    # get_filtered_vcf(data, 8)
    # get_vcf_with_indels_resolved(data, 8)

    # print(data.populations)

    # chunks = get_chunked_loci_bed(data, 10)
    # vchunks = []
    # for chunk in chunks:
    #     v = get_group_chunk_vcf(data, chunk, 2)
    #     vchunks.append(v)

    # vdir = data.stepdir / "vcfs"
    # vchunks = sorted(vdir.glob("chunk*.vcf.gz"))

    # raw_vcf = get_concat_chunk_vcfs(data, vchunks, 8)
    # filt_vcf = get_filtered_vcf(data, 8)
    # vcf = get_merged_vcf(data)
    # fvcf = get_filtered_vcf(data, 8)
