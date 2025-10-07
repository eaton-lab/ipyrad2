#!/usr/bin/env python

"""Functions to delimit loci beds by sample coverage

"""

from typing import List, Dict, Any
import sys
import tempfile
from pathlib import Path
import numpy as np
from loguru import logger
from ..utils.parallel2 import run_pipeline

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")


def get_name_from_bam(bam_file: Path) -> str:
    cmd = [BIN_SAM, "samples", bam_file]
    _, out, _ = run_pipeline([cmd])
    return out.decode().strip().split()[0]


def samtools_index_reference(reference: Path, threads: int) -> None:
    """Index reference with samtools."""
    cmd = [BIN_SAM, "faidx", reference]
    run_pipeline([cmd])
    return


def get_reference_sort_order(reference: Path, outdir: Path) -> Path:
    """Get scaff order from sam indexed REF file.
    """
    # destination file
    out_path = outdir / "REF_info.txt"

    # write fai file if it doesn't exist
    fai_path = reference.with_suffix(reference.suffix + ".fai")
    if not fai_path.exists():
        samtools_index_reference(reference, 4)

    # write REF_info file.
    cmd = ["cut", "-f", "1,2", str(fai_path)]
    run_pipeline([cmd], out_path)
    return out_path


def get_fragment_beds(sname: str, bam_file: Path, min_map_q: int, threads: int, outdir: Path) -> Path:
    """Produce a fragments BED (full inserts) from a coordinate-sorted BAM.

    Shell command
    -------------
    >>> $ samtools collate -u -@ 2 -O S1.sorted.bam \
    >>>   | bedtools bamtobed -bedpe -i - \
    >>>   | awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}' \
    >>>   | bedtools sort -i - > S1.fragments.bed
    """
    bed_dir = outdir / "beds"
    out_path = bed_dir / f"{sname}.fragments.bed"
    bed_dir.mkdir(parents=True, exist_ok=True)
    collate_path = outdir / f"{sname}.collate.bam"

    cmd1 = [
        BIN_SAM, "collate",
        "-@", str(threads),
        "-T", str(outdir / f"{sname}"),
        "-o", str(collate_path),
        str(bam_file),
    ]
    run_pipeline([cmd1])

    cmd2 = [BIN_SAM, "view", "-u", "-q", str(min_map_q), str(collate_path)]
    cmd3 = [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]
    cmd4 = ["awk", r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}']
    cmd5 = [BIN_BED, "sort", "-i", "-"]
    run_pipeline([cmd2, cmd3, cmd4, cmd5], out_path)
    collate_path.unlink()

    return out_path

    # # Pipeline commands
    # # sort into pairs
    # cmd1 = [BIN_SAM, "collate", "-u", "-@", str(threads), "-O", "-T", str(outdir / f"{sname}"), str(bam_file)]
    # logger.debug(" ".join(cmd1))
    # # apply read-level quality filter; NOTE: same filter is applied in variants.py
    # cmd2 = [BIN_SAM, "view", "-u", "-q", str(min_map_q), "-"]
    # logger.debug(" ".join(cmd2))
    # # measure beds including PE insert region. Note, if filter removed one of
    # # the pairs this will skip the other and report a warning to stderr, that
    # # is OK, we want both pairs excluded if one has low quality mapping.
    # cmd3 = [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]

    # # 3) collapse each pair to its fragment span [min(start), max(end)]
    # cmd4 = ["awk", r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}']

    # # 4) sort BED
    # cmd5 = [BIN_BED, "sort", "-i", "-"]

    # run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5], out_path)
    # return out_path


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
        "-g", str(fai_path),          # genome file to define chrom lens
        "-bg",                        # report depth in bedgraph format
    ]
    run_pipeline([cmd], out_path)
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
    run_pipeline([cmd1, cmd2], out_path)
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
    ref_info = outdir / "REF_info.txt"
    bed_dir = outdir / "beds"
    bed_files = [bed_dir / f"{sname}.fragments.merged.bed" for sname in names]
    bed_path = bed_dir / "loci.bed"

    # write genome sorted copy of each bed file
    sorted_paths = []
    with tempfile.TemporaryDirectory(prefix="bedmerge_") as tmpd:
        for i, src in enumerate(bed_files):
            dst = Path(tmpd) / f"{i:04d}_{src.name}.sorted.bed"
            sort_cmd = [BIN_BED, "sort", "-g", ref_info, "-i", str(src)]
            run_pipeline([sort_cmd], dst)
            sorted_paths.append(dst)

            # # write sorted copy to tempdir
            # with open(dst, "wb") as out, open(log_dir / f"sort_{i}.err", "wb") as err:
            #     rc = sp.run(sort_cmd, stdout=out, stderr=err).returncode
            # if rc != 0:
            #     raise RuntimeError(f"bedtools sort failed on {src} (see {err.name})")
            # sorted_paths.append(dst)

        # cmd1: bedtools multiinter
        cmd1 = [BIN_BED, "multiinter", "-i"] + [str(p) for p in sorted_paths] + ["-names"] + names

        # cmd2: threshold by K
        cmd2 = ["awk", f'BEGIN{{OFS="\\t"}} $4>={int(min_sample_coverage)} {{print $1,$2,$3,$4,$5}}']

        # cmd3: merge sub-intervals, keeping min support and distinct sample list
        cmd3 = [
            BIN_BED, "merge",
            "-i", "-",
            "-d", str(int(min_merge_distance)),
            "-c", "4",
            "-o", "min",
        ]

        # cmd4: filter intervals shorter than min_len (default=20)
        cmd4 = ["awk", "-v", f"L={min_locus_length}", 'BEGIN{OFS=FS="\t"} ($3-$2) >= L']

        # run pipeline
        run_pipeline([cmd1, cmd2, cmd3, cmd4], bed_path)
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
    _, out, _ = run_pipeline([cmd1, cmd2])

    stats = {
        "nloci_with_nonzero_mapping": 0,
        "median_depth_per_locus_with_nonzero_mapping": 0,
        "median_depth_per_locus_total": 0,
    }

    # parse stdout of cut
    coverages = out.decode().strip().split("\n")

    if not coverages:
        return stats

    # get nloci with non-zero coverage
    covs = np.array(list(map(int, coverages)))
    stats["nloci_with_nonzero_mapping"] = int(np.sum(covs > 0))
    stats["median_depth_per_locus_with_nonzero_mapping"] = float(np.median(covs[covs > 0]))
    stats["mean_depth_per_locus_with_nonzero_mapping"] = float(np.mean(covs[covs > 0]))
    stats["std_depth_per_locus_with_nonzero_mapping"] = float(np.std(covs[covs > 0]))
    stats["median_depth_per_locus_total"] = float(np.median(covs))
    stats["mean_depth_per_locus_total"] = float(np.mean(covs))
    stats["std_depth_per_locus_total"] = float(np.std(covs))
    return stats




if __name__ == "__main__":

    # tmpdir = Path("/home/deren/Documents/ipyrad-tests/Ama-out/tmpdir/")
    # bam = Path("/home/deren/Documents/ipyrad-tests/Ama-map/SLH_AL_3065.marked.sorted.bam")
    # print(get_sample_coverage_stats_in_loci_bed(bam, tmpdir))
    a = "TTGAAGACTGCTCTGTGCACAACCATCTAATAGTCGATTGTCCGACGTCGAGTGTGCAGTTTCTCGAGAAACAGCTCGTATCACGGGCCGGTTTCTTAGCATGCAATATGTGGGCATAATTCTCCTACCTTCTTCCGTTAACTGGTAACGTGACACAACAGGTGGCGAGTGTTTACCATCCAT"
    print(len(a))
