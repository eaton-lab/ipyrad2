

import sys
from pathlib import Path
from loguru import logger
import pandas as pd
import numpy as np
from ipyrad2.utils.parallel import run_pipeline


BIN_BED = str(Path(sys.prefix) / "bin" / "bedtools")


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


def make_lowdepth_mask(sname: str, min_sample_depth: int, tmpdir: Path):
    """Build a per-bp mask of positions inside `loci_bed` where bedGraph depth < min_depth.

    Output mask contains only the A (loci) columns and is split into minimal sub-intervals
    where coverage is below threshold (including 0-coverage gaps).
    """
    bed_dir = tmpdir / "beds"
    loci_bed = bed_dir / "loci.bed"
    ref_info = tmpdir / "REF_info.txt"
    sample_bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
    good_bed = bed_dir / f"{sname}.goodcov.bed"
    out_bed = bed_dir / f"{sname}.mask.bed"

    # 1) Threshold bedGraph: keep depth >= min_depth, drop depth column for set ops
    cmd1 = [
        "awk",
        f'BEGIN{{OFS="\\t"}} $4>={min_sample_depth} {{print $1,$2,$3}}',
        str(sample_bedgraph),
    ]
    cmd2 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    run_pipeline([cmd1, cmd2, cmd3], good_bed)

    # subtract good loci positions from all loci positions
    cmd1 = [
        BIN_BED, "subtract",
        "-a", str(loci_bed),
        "-b", str(good_bed),
        "-sorted",
    ]
    run_pipeline([cmd1], out_bed)
    return out_bed




if __name__ == "__main__":

    sname = "chumbia-JJ86-plate_J2"
    sname = "alaschanica-DE237-plate_J2"
    sname = "bella-JJ85-plate_J2"
    tmpdir = Path("/home/deren/Documents/ipyrad-tests/OUT4/tmpdir/")
    loci_bed = Path("/home/deren/Documents/ipyrad-tests/OUT4/tmpdir/beds/loci.bed")

    # s = get_sample_coverage_in_loci_bed(sname, loci_bed, tmpdir)
    # print(s)

    print(make_lowdepth_mask(sname, 5, tmpdir))
