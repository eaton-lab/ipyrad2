#!/usr/bin/env python

"""

"""

from typing import List, Dict, Any
import sys
import shlex
import tempfile
from pathlib import Path
import subprocess as sp
import numpy as np
from ..utils.exceptions import IPyradError

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")


def get_reference_sort_order(reference: Path, outdir: Path):
    """Get scaff order from sam indexed REF file.
    """
    out_path = outdir / "REF_info.txt"
    fai_path = reference.with_suffix(reference.suffix + ".fai")
    cmd = ["cut", "-f", "1,2", str(fai_path)]
    with open(out_path, 'wb') as out:
        p = sp.Popen(cmd, stderr=sp.PIPE, stdout=out)
        _, err = p.communicate()
        if p.returncode:
            raise IPyradError(err.decode())
    return out_path


def get_fragment_beds(sname: str, bam_file: Path, outdir: Path) -> Path:
    """Produce a fragments BED (full inserts) from a coordinate-sorted BAM.

    Shell command
    -------------
    >>> $ samtools collate -u -@ 2 -O S1.sorted.bam \
    >>>   | bedtools bamtobed -bedpe -i - \
    >>>   | awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}' \
    >>>   | bedtools sort -i - > S1.fragments.bed
    """
    bed_dir = outdir / "beds"
    log_dir = outdir / "logs"
    out_path = bed_dir / f"{sname}.fragments.bed"
    bed_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Pipeline commands
    # 1) name-group for -bedpe (stdout)
    cmd1 = [BIN_SAM, "collate", "-u", "-@", "1", "-O", str(bam_file)]

    # 2) paired-end to BEDPE
    cmd2 = [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]

    # 3) collapse each pair to its fragment span [min(start), max(end)]
    awk_prog = r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}'
    cmd3 = ["awk", awk_prog]

    # 4) sort BED
    cmd4 = [BIN_BED, "sort", "-i", "-"]

    # Open logs (avoid stderr PIPE backpressure)
    e1 = open(log_dir / f"{sname}.collate.err", "wb")
    e2 = open(log_dir / f"{sname}.bamtobed.err", "wb")
    e3 = open(log_dir / f"{sname}.awk.err", "wb")
    e4 = open(log_dir / f"{sname}.bedsort.err", "wb")

    with open(out_path, "wb") as OUT:
        p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=e1)
        upstream = p1.stdout

        p2 = sp.Popen(cmd2, stdin=upstream, stdout=sp.PIPE, stderr=e2)
        if upstream is not None:
            upstream.close()  # let p1 see EPIPE if p2 exits early

        p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=sp.PIPE, stderr=e3)
        p2.stdout.close()

        p4 = sp.Popen(cmd4, stdin=p3.stdout, stdout=OUT, stderr=e4)
        p3.stdout.close()

        rc4 = p4.wait()
        rc3 = p3.wait()
        rc2 = p2.wait()
        rc1 = p1.wait()

    # Close logs
    for fh in (e1, e2, e3, e4):
        try:
            fh.close()
        except Exception:
            pass

    # check for errors
    if any(rc != 0 for rc in (rc1, rc2, rc3, rc4)):
        # For easy debugging, show the exact commands
        cmds = "\n".join(shlex.join(c) for c in (cmd1, cmd2, cmd3, cmd4))
        raise RuntimeError(
            f"Fragment BED pipeline failed: samtools-collate={rc1}, bamtobed={rc2}, awk={rc3}, bedtools-sort={rc4}\n"
            f"Commands were:\n{cmds}\n"
            f"See logs in: {log_dir}"
        )
    return out_path


def get_fragment_coverage_beds(sname: str, reference: Path, outdir: Path) -> Path:
    """write depth filtered bed for each sample.

    >>> $ bedtools genomecov -i BED -g REF.scaflens -bg > fragments.bedgraph
    """
    # create a tmp file with REF scaffold length
    bed_dir = outdir / "beds"
    fragment_bed = bed_dir / f"{sname}.fragments.bed"
    out_path = bed_dir / f"{sname}.fragments.bedgraph"
    fai_path = reference.with_suffix(reference.suffix + ".fai")
    assert fai_path.exists(), "must call `samtools faidx $REF`"

    # get bedgraph format for storing depths
    cmd = [
        BIN_BED, "genomecov",
        "-i", str(fragment_bed),
        "-g", str(fai_path),
        "-bg",
    ]

    # write to bed
    with open(out_path, 'wb') as out:
        p1 = sp.Popen(cmd, stderr=sp.PIPE, stdout=out)
        _, err = p1.communicate()
    # remove non-merged bed and return merged bed
    fragment_bed.unlink()
    return out_path


def get_fragment_merged_coverage_beds(sname: str, outdir: Path):
    """write bed with intervals of coverage above {min_depth_majrule}.

    >>> $ awk -v MIN=3 '$4>=MIN' sname.fragments.bedgraph \
    >>>   | bedtools merge -i - > sname.loci.min3.bed
    """
    # paths
    bed_dir = outdir / "beds"
    bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
    out_path = bed_dir / f"{sname}.fragments.merged.bed"

    # keep all RAD beds above depth=1 and merge
    cmd1 = ["awk", "-v", "MIN=1", r'$4>=MIN', bedgraph]
    cmd2 = [BIN_BED, "merge", "-i", "-"]

    # write to
    with open(out_path, 'wb') as out:
        p1 = sp.Popen(cmd1, stderr=sp.PIPE, stdout=sp.PIPE)
        p2 = sp.Popen(cmd2, stderr=sp.PIPE, stdout=out, stdin=p1.stdout)
        p1.stdout.close()
        _, err2 = p2.communicate()
        _, err1 = p1.communicate()
    # ...
    # Check in reverse order to surface the first failing stage
    if p2.returncode:
        raise IPyradError(f"bedtools merge failed ({p2.returncode}).\n{err2.decode(errors='ignore')}")
    if p1.returncode:
        raise IPyradError(f"bwa mem failed ({p1.returncode}).\n{err1.decode(errors='ignore')}")
    # bedgraph.unlink()
    return out_path


def get_across_sample_loci_bed(
    names: List[str],
    min_sample_coverage: int,
    min_merge_distance: int,
    min_locus_length: int,
    outdir: Path,
) -> Dict[str, Any]:
    """Merge beds across samples to get joint bed regions (loci)

    Require at least sample coverage of 3 (with the ref makes 4).
    - sort beds
    - count sample cov using multiinter
    - drop low cov regions
    - merge remaining nearbys
    """
    bed_files = [
        outdir / "beds" / f"{sname}.fragments.merged.bed"
        for sname in names
    ]

    # paths
    ref_info = outdir / "REF_info.txt"
    bed_path = outdir / "beds" / "loci.bed"
    log_dir = outdir / "logs"

    # write genome sorted copy of each bed file
    sorted_paths = []
    with tempfile.TemporaryDirectory(prefix="bedmerge_") as tmpd:
        for i, src in enumerate(bed_files):
            dst = Path(tmpd) / f"{i:04d}_{src.name}.sorted.bed"
            sort_cmd = [BIN_BED, "sort", "-g", ref_info, "-i", str(src)]

            # write sorted copy to tempdir
            with open(dst, "wb") as out, open(log_dir / f"sort_{i}.err", "wb") as err:
                rc = sp.run(sort_cmd, stdout=out, stderr=err).returncode
            if rc != 0:
                raise RuntimeError(f"bedtools sort failed on {src} (see {err.name})")
            sorted_paths.append(dst)

        # cmd1: bedtools multiinter
        cmd1 = [BIN_BED, "multiinter", "-i"] + [str(p) for p in sorted_paths] + ["-names"] + names

        # cmd2: threshold by K
        awk1 = f'BEGIN{{OFS="\\t"}} $4>={int(min_sample_coverage)} {{print $1,$2,$3,$4,$5}}'
        cmd2 = ["awk", awk1]

        # cmd3: merge sub-intervals, keeping min support and distinct sample list
        cmd3 = [
            BIN_BED, "merge",
            "-i", "-",
            "-d", str(int(min_merge_distance)),
            "-c", "4",
            "-o", "min",
        ]

        # cmd4: filter intervals shorter than min_len (default=20)
        awk2 = 'BEGIN{OFS=FS="\t"} ($3-$2) >= L'
        cmd4 = ["awk", "-v", f"L={min_locus_length}", awk2]

        # cmd4: add count of distinct samples (col5 is the csv list)
        # awk2 = 'BEGIN{OFS="\\t"}{n=($5==""?0:split($5,a,",")); print $1,$2,$3,$4,n,$5}'
        # cmd4 = ["awk", awk2]
        # logger.debug(" ".join(cmd4))

        # Open logs to avoid stderr PIPE backpressure
        e1 = open(log_dir / "multiinter.err", "wb")
        e2 = open(log_dir / "awk_k.err", "wb")
        e3 = open(log_dir / "merge.err", "wb")
        e4 = open(log_dir / "awk_loci.err", "wb")

        with open(bed_path, "wb") as OUT:
            p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=e1)
            p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, stderr=e2)
            if p1.stdout:
                p1.stdout.close()

            p3 = sp.Popen(cmd3, stdin=p2.stdout, stdout=sp.PIPE, stderr=e3)
            if p2.stdout:
                p2.stdout.close()

            p4 = sp.Popen(cmd4, stdin=p3.stdout, stdout=OUT, stderr=e4)
            if p3.stdout:
                p3.stdout.close()

            rc4 = p4.wait()
            rc3 = p3.wait()
            rc2 = p2.wait()
            rc1 = p1.wait()

        # close logs
        for fh in (e1, e2, e3, e4):
            try:
                fh.close()
            except Exception:
                pass

        if any(rc != 0 for rc in (rc1, rc2, rc3, rc4)):
            cmds = "\n".join(shlex.join(c) for c in (cmd1, cmd2, cmd3, cmd4))
            raise RuntimeError(
                "merge_beds_with_support failed with return codes: "
                f"multiinter={rc1}, awk1={rc2}, merge={rc3}, awk2={rc4}\n"
                f"Commands were:\n{cmds}\nLogs in: {log_dir}"
            )
    return bed_path


def get_sample_coverage_stats_in_loci_bed(bam_file: Path, outdir: Path) -> Dict[str, float]:
    """Return dict with stats of sampling mapping per locus bed.
    """
    loci_bed = outdir / "beds" / "loci.bed"

    # commands
    cmd1 = [
        BIN_BED, "coverage",
        "-a", str(loci_bed),
        "-b", str(bam_file),
        "-counts",
    ]
    cmd2 = ["cut", "-f", "5"]

    # run pipeline
    p1 = sp.Popen(cmd1, stdout=sp.PIPE, stderr=sp.DEVNULL)
    p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.DEVNULL)
    out, err = p2.communicate()
    if p2.returncode:
        raise RuntimeError(f"error in stats coverage {cmd1}: {err.decode()}")
    # coverage = [int]
    covs = np.array(list(map(int, out.decode().strip().split("\n"))))

    # get nloci with non-zero coverage
    results = {
        "nloci_with_nonzero_mapping": int(np.sum(covs > 0)),
        "median_depth_per_locus_with_nonzero_mapping": float(np.median(covs[covs > 0])),
        "median_depth_per_locus_total": float(np.median(covs)),
        "covs": covs > 0,
    }
    return results




if __name__ == "__main__":

    pass
