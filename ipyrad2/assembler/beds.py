#!/usr/bin/env python

"""Functions to delimit loci beds by sample coverage

"""

from typing import List, Dict, Any
import sys
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from loguru import logger
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline

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


def get_reference_sort_order(reference: Path, tmpdir: Path) -> Path:
    """Get scaff order from sam indexed REF file.
    """
    # destination file
    out_path = tmpdir / "REF_info.txt"

    # write fai file if it doesn't exist
    fai_path = reference.with_suffix(reference.suffix + ".fai")
    if not fai_path.exists():
        samtools_index_reference(reference, 4)

    # write REF_info file.
    cmd = ["cut", "-f", "1,2", str(fai_path)]
    run_pipeline([cmd], out_path)
    return out_path


def get_coverage_bed_graphs(sname: str, bam_file: Path, reference: Path, tmpdir: Path, min_map_q: int, min_sample_depth: int, min_merge_distance: int, threads: int):
    r"""Produce a fragments BED (full inserts) from a coordinate-sorted BAM.

    Shell command
    -------------
    >>> $ samtools collate -u -@ 2 -O S1.sorted.bam \
    >>>   | bedtools bamtobed -bedpe -i - \
    >>>   | awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}' \
    >>>   | bedtools sort -i - \
    >>>   | bedtools genomecov -i - -g REF.info -bg > S1.bedgraph
    """
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True, exist_ok=True)
    coll_dir = tmpdir / f"{sname}.collate"
    coll_dir.mkdir(exist_ok=True)
    out_bed_count = bed_dir / f"{sname}.fragments.bedgraph"    # has counts
    out_bed_merge = bed_dir / f"{sname}.fragments.merged.bed"  # merged, no counts
    fai_path = reference.with_suffix(reference.suffix + ".fai")

    # Test for pe vs se bam file by checking PAIRED bit of the first few reads
    cmd1 = [BIN_SAM, "view", "-f", "0x1", bam_file]
    cmd2 = ["head", "-n", "1000"]
    cmd3 = ["wc", "-l"]
    _, ct, _ = run_pipeline([cmd1, cmd2, cmd3])

    is_paired = False
    # This is very simple. Count the number of lines in the samtools view call.
    if int(ct.strip()) > 0:
        is_paired = True

    # CHECKED that this properly pipes on large files. It does.
    # Note: collate has heavy I/O writing tmp files here on large WGS files.
    # this puts F/R pairs next to each other.
    cmd1 = [
        BIN_SAM, "collate",
        "-@", str(min(threads, 4)),             # doesn't benefit from >4
        "-T", str(coll_dir / f"{sname}"),
        "-r", "1000000",
        "-u",
        "-O",
        str(bam_file),
    ]
    # compute bedpe table for each pair record (skipping pairs if one was filtered out), e.g.,
    # Chr1    60908424        60908533        Chr1    60908434        60908543        LH00150:341:22HGMLLT3:3:1160:35876:20034        49      -       +
    # Chr3    5915120         5915253         Chr3    5915131         5915265         LH00150:341:22HGMLLT3:3:2236:23034:9892         57      +       -
    # Chr2    109898149       109898235       Chr2    109898287       109898367       LH00150:341:22HGMLLT3:3:1287:23015:18000        54      +       -
    cmd2 = [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]
    if is_paired:
        cmd2.insert(2, "-bedpe")
    # extract only the records from this table where the mean mapq passes this steps filter
    # Chr1    60908424        60908533        Chr1    60908434        60908543        LH00150:341:22HGMLLT3:3:1160:35876:20034        49      -       +
    # Chr3    5915120         5915253         Chr3    5915131         5915265         LH00150:341:22HGMLLT3:3:2236:23034:9892         57      +       -
    # Chr2    109898149       109898235       Chr2    109898287       109898367       LH00150:341:22HGMLLT3:3:1287:23015:18000        54      +       -
    cmd3 = ["awk", "-v", f'q={min_map_q}']
    if is_paired:
        cmd3 += [r'BEGIN{OFS="\t"} ($8+0) >= q']
    else:
        cmd3 += [r'BEGIN{OFS="\t"} ($5+0) >= q']
    # check start/end values and extract only (chrom, start, end)
    # Chr4    108577598       108577715
    # Chr4    106107228       106107344
    # Chr1    45721044        45721144
    if is_paired:
        cmd4 = ["awk", r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}']
    else:
        cmd4 = ["awk", r'BEGIN{OFS="\t"} {print $1,$2,$3}']
    # sort beds by genome coordinates
    # Chr1    825268  825547
    # Chr1    825268  825547
    # Chr1    825268  825547
    # Chr1    833321  833418
    # Chr1    833321  833418
    # Chr1    833321  833418
    cmd5 = [BIN_BED, "sort", "-i", "-", "-faidx", str(fai_path)]
    # get coverage in bedgraph format with counts within regions
    # Chr1    833321  833418  14
    # Chr1    833419  833520  2
    # Chr1    837052  837165  4
    # Chr1    837165  837230  18
    # Chr1    837230  837240  14
    cmd6 = [BIN_BED, "genomecov", "-i", "-", "-g", str(fai_path), "-bg"]
    # filter out sites below MIN_DEPTH coverage
    # Chr1    833321  833418  14
    # Chr1    837052  837165  4
    # Chr1    837165  837230  18
    # Chr1    837230  837240  14
    cmd7 = ["awk", "-v", f"MIN={min_sample_depth}", r'$4>=MIN', "-"]
    # pipe one stream forward and save another to file.
    cmd8 = ["tee", str(out_bed_count)]
    # merge beds within MIN_MERGE_DISTANCE of each other
    cmd9 = [BIN_BED, "merge", "-d", str(min_merge_distance), "-i", "-"]
    # Chr1    833321  833418
    # Chr1    837052  837240
    cmds = []
    logger.trace(repr(cmds))
    ret = run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9], out_bed_merge)
    logger.trace(ret)

    shutil.rmtree(coll_dir)
    logger.debug(f"wrote bed graph for {sname}")
    return out_bed_merge


# ======= deprecated SE code now incorporated above.
#     # stream the bed for each read/read pair. toggle `-bedpe` only for PE data
#     cmd2 = [BIN_BED, "bamtobed"]
#     if is_paired:
#         cmd2.append("-bedpe")
#     cmd2.extend(["-i", "-"])

#     # filter pairs on mapq: applies same mapq min as in variants.py
#     # mapq in bed for se is column 5; for pe is column 8
#     qcol = 5
#     if is_paired:
#         qcol = 8
#     cmd3 = ["awk", "-v", f'q={min_map_q}', r'BEGIN{OFS="\t"} ' + f'(${qcol}+0) >= q']

#     if is_paired:
#         # check and pull out only chrom, start, end
#         cmd4 = ["awk", r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}']
#     else:
#         # For SE data we just need the first 3 columns of the input bed
#         cmd4 = ["awk", r'BEGIN{OFS="\t"} {print $1,$2,$3}']

#     # sort beds
#     cmd5 = [BIN_BED, "sort", "-i", "-"]
#     # get coverage for each site from overlapping beds
#     cmd6 = [BIN_BED, "genomecov", "-i", "-", "-g", str(fai_path), "-bg"]
#     # pipeline
#     for cmd in [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6]:
#         logger.debug(" ".join(cmd))
#     run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6], out_path)
#     shutil.rmtree(coll_dir)
#     logger.debug(f"wrote bed graph for {sname}")
#     return out_path



# def get_fragment_beds(sname: str, bam_file: Path, threads: int, tmpdir: Path) -> Path:
#     r"""Produce a fragments BED (full inserts) from a coordinate-sorted BAM.

#     Shell command
#     -------------
#     >>> $ samtools collate -u -@ 2 -O S1.sorted.bam \
#     >>>   | bedtools bamtobed -bedpe -i - \
#     >>>   | awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}' \
#     >>>   | bedtools sort -i - > S1.fragments.bed
#     """
#     bed_dir = tmpdir / "beds"
#     out_path = bed_dir / f"{sname}.fragments.bed"
#     bed_dir.mkdir(parents=True, exist_ok=True)

#     # CHECKED that this properly pipes on large files. It does.
#     # Note: collate has heavy I/O writing tmp files here.
#     coll_dir = tmpdir / f"{sname}.collate"
#     coll_dir.mkdir(exist_ok=True)
#     cmd1 = [
#         BIN_SAM, "collate",
#         "-@", str(min(threads, 4)),             # doesn't benefit from >4
#         "-T", str(coll_dir / f"{sname}"),
#         "-r", "1000000",
#         "-u",
#         "-O",
#         str(bam_file),
#     ]
#     cmd2 = [BIN_BED, "bamtobed", "-bedpe", "-i", "-"]
#     cmd3 = ["awk", r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}']
#     cmd4 = [BIN_BED, "sort", "-i", "-"]
#     run_pipeline([cmd1, cmd2, cmd3, cmd4], out_path)
#     shutil.rmtree(coll_dir)
#     logger.debug(f"wrote fragment beds for {sname}")
#     return out_path


# def get_fragment_coverage_beds(sname: str, reference: Path, tmpdir: Path) -> Path:
#     """write depth filtered bed for each sample.

#     >>> $ bedtools genomecov -i BED -g REF.scaflens -bg > fragments.bedgraph
#     """
#     # create a tmp file with REF scaffold length
#     bed_dir = tmpdir / "beds"
#     fragment_bed = bed_dir / f"{sname}.fragments.bed"
#     out_path = bed_dir / f"{sname}.fragments.bedgraph"
#     fai_path = reference.with_suffix(reference.suffix + ".fai")
#     assert fai_path.exists(), "must call `samtools faidx $REF`"

#     # get bedgraph format for storing depths
#     cmd = [
#         BIN_BED, "genomecov",
#         "-i", str(fragment_bed),
#         "-g", str(fai_path),          # genome file to define chrom lens
#         "-bg",                        # report depth in bedgraph format
#     ]
#     run_pipeline([cmd], out_path)
#     fragment_bed.unlink()
#     return out_path


# def get_fragment_merged_coverage_beds(sname: str, tmpdir: Path):
#     """write bed with intervals of coverage above {min_depth_majrule}.

#     >>> $ awk -v MIN=3 '$4>=MIN' sname.fragments.bedgraph \
#     >>>   | bedtools merge -i - > sname.loci.min3.bed
#     """
#     # paths
#     bed_dir = tmpdir / "beds"
#     bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
#     out_path = bed_dir / f"{sname}.fragments.merged.bed"

#     # keep all RAD beds above depth=1 and merge
#     # TODO: The docstring says this honors min_depth_majrule, but here it doesn't
#     cmd1 = ["awk", "-v", "MIN=1", r'$4>=MIN', bedgraph]
#     cmd2 = [BIN_BED, "merge", "-i", "-"]
#     # -d merge beds within this many sites of each other.
#     run_pipeline([cmd1, cmd2], out_path)
#     return out_path
# >>>>>>> 24fdcf86e91a9aad562701d4ab2a9706ef0f1595


def get_across_sample_loci_bed(
    snames: List[str],
    min_sample_coverage: int,
    min_merge_distance: int,
    min_locus_length: int,
    tmpdir: Path,
) -> Dict[str, Any]:
    """Merge beds across samples to get joint bed regions (loci)

    Require at least sample coverage of 3 (with the ref makes 4).
    - sort beds
    - count sample cov using multiinter
    - drop low cov regions
    - merge remaining nearbys
    """
    ref_info = tmpdir / "REF_info.txt"
    bed_dir = tmpdir / "beds"
    bed_paths = [bed_dir / f"{sname}.fragments.merged.bed" for sname in snames]
    out_bed = bed_dir / "loci.bed"

    # for each bed get [chrom, start, end, nsamples, samplenames, A-present, B-present, ...]
    # Chr1    1789726 1789745 1       C       0       0       1       0
    # Chr1    1790639 1790658 1       B       0       1       0       0
    # Chr1    1792068 1792357 4       A,B,C,D 1       1       1       1
    # Chr1    1792357 1792384 3       A,B,C   1       1       1       0
    # Chr1    1792384 1792386 2       B,C     0       1       1       0
    # Chr1    1799627 1799701 1       B       0       1       0       0
    # Chr1    1810262 1810282 1       D       0       0       0       1
    cmd1 = [
        BIN_BED, "multiinter",
        "-g", str(ref_info)
    ]
    cmd1 += ["-i"] + [str(p) for p in bed_paths]
    cmd1 += ["-names"] + snames

    # require at least MIN_SAMPLES_COVERAGE in each bed and print only first 5 cols
    # Chr1    1792068 1792357 4       A,B,C,D
    # Chr1    1792357 1792384 3       A,B,C
    # Chr1    2344873 2344902 3       A,B,D
    # Chr1    2665674 2665760 3       B,C,D
    # Chr1    2824851 2824932 4       A,B,C,D
    # Chr1    3045768 3045944 3       A,B,D
    cmd2 = ["awk", f'BEGIN{{OFS="\\t"}} $4>={int(min_sample_coverage)} {{print $1,$2,$3,$4,$5}}']

    # merge sub-intervals by MIN_MERGE_DISTANCE
    # Chr1    1792068 1792384 4
    # Chr1    2344873 2344902 3
    # Chr1    2665674 2665760 3
    # Chr1    2824851 2824932 4
    # Chr1    3045768 3045944 3
    cmd3 = [
        BIN_BED, "merge",
        "-i", "-",
        "-d", str(int(min_merge_distance)),
        "-c", "4",
        "-o", "min",
    ]

    # cmd4: filter intervals shorter than MIN_LOCUS_LENGTH
    # Chr1    1792068 1792384 4
    # Chr1    2665674 2665760 3
    # Chr1    2824851 2824932 4
    # Chr1    3045768 3045944 3
    cmd4 = ["awk", "-v", f"L={min_locus_length}", 'BEGIN{OFS=FS="\t"} ($3-$2) >= L']

    # run pipeline
    run_pipeline([cmd1, cmd2, cmd3, cmd4], out_bed)
    return out_bed



# def get_fragment_merged_coverage_beds(sname: str, tmpdir: Path):
#     """write bed with intervals of coverage above {min_depth_majrule}.

#     >>> $ awk -v MIN=3 '$4>=MIN' sname.fragments.bedgraph \
#     >>>   | bedtools merge -i - > sname.loci.min3.bed
#     """
#     # paths
#     bed_dir = tmpdir / "beds"
#     bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
#     out_path = bed_dir / f"{sname}.fragments.merged.bed"

#     # keep all RAD beds above depth=1 and merge
#     # Chr1    833321  833418  14
#     # Chr1    837052  837165  4
#     # Chr1    837165  837230  18
#     # Chr1    837230  837240  14
#     cmd7 = ["awk", "-v", f"MIN={min_sample_depth}", r'$4>=MIN', "-"]
#     # merge book-ended windows with same depth into more concise interval format, e.g.,
#     # note that we will merge on larger distances later across-samples.
#     cmd8 = [BIN_BED, "merge", "-i", "-"]

#     cmd1 = ["awk", "-v", "MIN=1", r'$4>=MIN', bedgraph]
#     cmd2 = [BIN_BED, "merge", "-i", "-"]
#     # -d merge beds within this many sites of each other.
#     run_pipeline([cmd1, cmd2], out_path)
#     return out_path


def old_get_across_sample_loci_bed(
    names: List[str],
    min_sample_coverage: int,
    min_merge_distance: int,
    min_locus_length: int,
    tmpdir: Path,
) -> Dict[str, Any]:
    """Merge beds across samples to get joint bed regions (loci)

    Require at least sample coverage of 3 (with the ref makes 4).
    - sort beds
    - count sample cov using multiinter
    - drop low cov regions
    - merge remaining nearbys
    """
    ref_info = tmpdir / "REF_info.txt"
    bed_dir = tmpdir / "beds"
    bed_files = [bed_dir / f"{sname}.fragments.merged.bed" for sname in names]
    bed_path = bed_dir / "loci.bed"

    # write genome sorted copy of each bed file
    sorted_paths = []
    with tempfile.TemporaryDirectory(prefix="bedmerge_", delete=False) as tmpd:
        for i, src in enumerate(bed_files):
            dst = Path(tmpd) / f"{i:04d}_{src.name}.sorted.bed"
            sort_cmd = [
                BIN_BED, "sort",
                "-g", ref_info,
                "-i", str(src),
            ]
            run_pipeline([sort_cmd], dst)
            sorted_paths.append(dst)

        # cmd1: bedtools multiinter
        cmd1 = [
            BIN_BED, "multiinter",
            "-g", str(ref_info),
        "-i"] + [str(p) for p in sorted_paths] + ["-names"] + names

        # cmd2: threshold by K
        cmd2 = ["awk", f'BEGIN{{OFS="\\t"}} $4>={int(min_sample_coverage)} {{print $1,$2,$3,$4,$5}}']

        # cmd3: Sort again to ensure proper order for the merge
        cmd3 = [
            BIN_BED, "sort",
            "-g", str(ref_info),
            "-i", "-",
        ]

        # cmd4: merge sub-intervals, keeping min support and distinct sample list
        cmd4 = [
            BIN_BED, "merge",
            "-i", "-",
            "-d", str(int(min_merge_distance)),
            "-c", "4",
            "-o", "min",
        ]

        # cmd5: filter intervals shorter than min_len (default=20)
        cmd5 = ["awk", "-v", f"L={min_locus_length}", 'BEGIN{OFS=FS="\t"} ($3-$2) >= L']

        # run pipeline
        for cmd in [cmd1, cmd2, cmd3, cmd4, cmd5]:
            logger.debug(" ".join(cmd))
        run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5], bed_path)
    return bed_path



def get_sample_coverage_in_loci_bed(sname: str, loci_bed: Path, tmpdir: Path):
    """
    Get sample coverage from its bedgraph but filtered to only the
    loci that are in the loci_bed file, and compute stats.
    """
    cov_bed = tmpdir / "beds" / f"{sname}.fragments.bedgraph"

    # 1) Intersect with -wao (write all overlaps, including zero-overlap rows)
    s1 = [
        BIN_BED, "intersect",
        "-a", str(loci_bed),
        "-b", str(cov_bed),
        "-wao",
    ]
    _, out, _ = run_pipeline([s1])

    # parse intersect table as dataframe, label columns and convert types
    data = pd.DataFrame([i.split("\t") for i in out.decode().strip().split("\n")])
    data[3] = data[0] + ":" + data[1] + "-" + data[2]
    data.columns = ["chrom", "start", "end", "name", "c", "s", "e", "overlap", "coverage"]
    data.loc[data["overlap"] == '.', "overlap"] = 0
    data["overlap"] = data["overlap"].astype(int)
    data["coverage"] = data["coverage"].astype(int)
    data["weighted"] = data["overlap"] * data["coverage"]

    # aggregate by locus to get length-weighted coverage (of sites where there is coverage)
    covs = data.groupby("name").apply(lambda x: (sum(x["weighted"]) / sum(x["coverage"])) if sum(x["coverage"]) else 0, include_groups=False)

    # return stats dict and coverages
    nonzero = covs[covs > 0]
    stats = {
        "nloci": sum(covs > 0),
        "nloci_zero_cov": covs.size - sum(covs > 0),
        "mean_depth_per_nonzero_cov_locus": float(np.mean(nonzero)),
        "median_depth_per_nonzero_cov_locus": float(np.median(nonzero)),
        "stdev_depth_per_nonzero_cov_locus": float(np.std(nonzero)),
        "mean_depth_per_locus_total": float(np.mean(covs)),
        "median_depth_per_locus_total": float(np.median(covs)),
        "stdev_depth_per_locus_total": float(np.std(covs)),
    }
    return stats




def get_sample_coverage_stats_in_loci_bed(bam_file: Path, loci_bed: Path, min_map_q: int, ref_info: Path) -> Dict[str, float]:
    """Return dict with stats of sampling mapping per locus bed.
    """
    # this shouldn't happen, but sanity check.
    if not bam_file.exists():
        raise IPyradError(f"bam file {bam_file} does not exist.")

    # apply mapq filter again in case it was increased since the map step.
    cmd1 = [
        BIN_SAM, "view",
        "-q", str(min_map_q),
        "-u",
        str(bam_file),
    ]

    # compute coverage table ...
    cmd2 = [
        BIN_BED, "coverage",
        "-a", str(loci_bed),
        "-b", "-",
        "-g", str(ref_info),
        "-sorted",
        "-counts",  # TODO: do we want counts or means?
    ]
    cmd3 = ["cut", "-f", "5"]
    _, out, _ = run_pipeline([cmd1, cmd2, cmd3])

    stats = {
        "nloci": 0,
        "mean_depth_per_locus_with_nonzero_mapping": 0,
        "median_depth_per_locus_with_nonzero_mapping": 0,
        "std_depth_per_locus_with_nonzero_mapping": 0,
        "mean_depth_per_locus_total": 0,
        "median_depth_per_locus_total": 0,
        "std_depth_per_locus_total": 0,
    }

    # parse stdout of cut
    coverages = out.decode().strip().split("\n")
    del out

    # no loci has sufficient sample coverage
    # TODO: raise warning and handle this instaed.
    if (coverages[0] == ""):
        raise IPyradError(f"{bam_file.name} has no regions in {loci_bed}")

    # get nloci with non-zero coverage
    covs = np.array(list(map(int, coverages)))
    if not sum(covs):
        return stats, np.zeros(len(coverages))
    del coverages
    stats["nloci"] = int(np.sum(covs > 0))
    stats["median_depth_per_locus_with_nonzero_mapping"] = float(np.median(covs[covs > 0]))
    stats["mean_depth_per_locus_with_nonzero_mapping"] = float(np.mean(covs[covs > 0]))
    stats["std_depth_per_locus_with_nonzero_mapping"] = float(np.std(covs[covs > 0]))
    stats["median_depth_per_locus_total"] = float(np.median(covs))
    stats["mean_depth_per_locus_total"] = float(np.mean(covs))
    stats["std_depth_per_locus_total"] = float(np.std(covs))
    _cmeans = covs[covs > 0].mean()
    _cstds = np.clip(covs[covs > 0].std(), a_min=1.0, a_max=None)
    read_depth_zscores = abs(covs - _cmeans) / _cstds
    return stats, read_depth_zscores




if __name__ == "__main__":

    # tmpdir = Path("/home/deren/Documents/ipyrad-tests/Ama-out/tmpdir/")
    # bam = Path("/home/deren/Documents/ipyrad-tests/Ama-map/SLH_AL_3065.marked.sorted.bam")
    # print(get_sample_coverage_stats_in_loci_bed(bam, tmpdir))
    a = "TTGAAGACTGCTCTGTGCACAACCATCTAATAGTCGATTGTCCGACGTCGAGTGTGCAGTTTCTCGAGAAACAGCTCGTATCACGGGCCGGTTTCTTAGCATGCAATATGTGGGCATAATTCTCCTACCTTCTTCCGTTAACTGGTAACGTGACACAACAGGTGGCGAGTGTTTACCATCCAT"
    print(len(a))
