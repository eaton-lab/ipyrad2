#!/usr/bin/env python

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from ipyrad2.utils.parallel import run_pipeline
from ipyrad2.utils.seqs import revcomp
import json

BIN_SAMTOOLS = Path(sys.prefix) / "bin" / "samtools"


DIR = Path("/home/deren/Documents/ipyrad-tests/tmp_map/")
bam_file = DIR / "SLH_AL_0001-restricted.filtered.bam"
F1 = DIR / "SLH_AL_0001-restricted.tmp.stats1.json"
FD = DIR / "SLH_AL_0001-restricted.tmp.stats_dups.json"
F2 = DIR / "SLH_AL_0001-restricted.tmp.stats2.json"
F3 = DIR / "SLH_AL_0001-restricted.tmp.stats3.json"


if __name__ == "__main__":

    with F1.open('r') as indata:
        d1 = json.loads(indata.read())
    with FD.open('r') as indata:
        a = d1['records_filter_accepted']
        b = [int(i.split()[-1]) for i in indata.readlines() if i.startswith("DUPLICATE TOTAL")][0]
        dd = {
            'records_processed': a,
            'records_filter_accepted': a - b,
            'records_filter_rejected': b,
        }
    with F2.open('r') as indata:
        d2 = json.loads(indata.read())
    with F3.open('r') as indata:
        d3 = json.loads(indata.read())

    # get mean, std of mapq
    cmd1 = [BIN_SAMTOOLS, "view", str(bam_file)]
    cmd2 = ["cut", "-f", "5"]
    _, out, _ = run_pipeline([cmd1, cmd2])
    out = np.array(list(map(int, out.decode().strip().split())))
    if out.size:
        mean_mapq = np.mean(out)
        median_mapq = np.median(out)
        stdev_mapq = np.std(out)
    else:
        mean_mapq = float('nan')
        median_mapq = float('nan')
        stdev_mapq = float('nan')

    data = {
        "nreads_processed": d1["records_processed"],
        "nreads_filtered_by_not_primary": d1["records_filter_rejected"],
        "nreads_filtered_by_min_mapq": d2["records_filter_rejected"],
        "nreads_filtered_by_bad_pairing": d3["records_filter_rejected"],
        "nreads_passed_filters": d3["records_filter_accepted"],
        "mapq_mean_after_filters": float(mean_mapq),
        "mapq_median_after_filters": float(median_mapq),
        "mapq_stdev_after_filters": float(stdev_mapq),
    }


    j = {'a': data, 'b': data}
    print(pd.DataFrame(j).T)


