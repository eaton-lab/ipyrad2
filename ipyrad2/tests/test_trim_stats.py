

from typing import List
from pathlib import Path
import pandas as pd
import json


def get_stats_dicts(snames: List[Path], outdir: Path):
    """..."""
    jdata = {}
    snames = sorted(snames)
    for sname in snames:
        stats_file = outdir / f"{sname}.stats.json"
        if stats_file.exists():
            with open(stats_file, 'r', encoding="utf-8") as indata:
                jdata[sname] = json.loads(indata.read())
    df = pd.DataFrame(index=snames, columns=[
        "total_reads_before", "total_bases_before", "q20_rate_before", "q30_rate_before",
        "read1_mean_length_before", "read2_mean_length_before",
        "total_reads_after", "total_bases_after", "q20_rate_after", "q30_rate_after",
        "read1_mean_length_after", "read2_mean_length_after",
        "reads_filtered_by_low_quality",
        "reads_filtered_by_too_many_N",
        "reads_filtered_by_low_complexity",
        "reads_filtered_by_too_short",
        "adapter_trimmed_reads",
        "adapter_trimmed_bases",
    ])

    for sname in snames:
        # print(sname, )
        j = jdata[sname]
        df.loc[sname, "total_reads_before"] = j["summary"]["before_filtering"]["total_reads"]
        df.loc[sname, "total_bases_before"] = j["summary"]["before_filtering"]["total_bases"]
        df.loc[sname, "q20_rate_before"] = j["summary"]["before_filtering"]["q20_rate"]
        df.loc[sname, "q30_rate_before"] = j["summary"]["before_filtering"]["q30_rate"]
        df.loc[sname, "read1_mean_length_before"] = j["summary"]["before_filtering"]["read1_mean_length"]
        df.loc[sname, "read2_mean_length_before"] = j["summary"]["before_filtering"]["read2_mean_length"]
        df.loc[sname, "total_reads_after"] = j["summary"]["after_filtering"]["total_reads"]
        df.loc[sname, "total_bases_after"] = j["summary"]["after_filtering"]["total_bases"]
        df.loc[sname, "q20_rate_after"] = j["summary"]["after_filtering"]["q20_rate"]
        df.loc[sname, "q30_rate_after"] = j["summary"]["after_filtering"]["q30_rate"]
        df.loc[sname, "read1_mean_length_after"] = j["summary"]["after_filtering"]["read1_mean_length"]
        df.loc[sname, "read2_mean_length_after"] = j["summary"]["after_filtering"]["read2_mean_length"]
        df.loc[sname, "reads_filtered_by_low_quality"] = j["filtering_result"]["low_quality_reads"]
        df.loc[sname, "reads_filtered_by_too_many_N"] = j["filtering_result"]["too_many_N_reads"]
        df.loc[sname, "reads_filtered_by_low_complexity"] = j["filtering_result"]["low_complexity_reads"]
        df.loc[sname, "reads_filtered_by_too_short"] = j["filtering_result"]["too_short_reads"]
        df.loc[sname, "adapter_trimmed_reads"] = j["adapter_cutting"]["adapter_trimmed_reads"]
        df.loc[sname, "adapter_trimmed_bases"] = j["adapter_cutting"]["adapter_trimmed_bases"]
    return df


DIR = Path("/home/deren/Documents/ipyrad-tests/Ped_trim_Oct8")
snames = ["cinerascens-DE719-plate_J2", "chumbia-JJ86-plate_J2", "bella-JJ85-plate_J2"]
jdata = get_stats_dicts(snames, DIR)
print(jdata.to_string("/tmp/test.tsv", float_format=lambda x: f"{x:.6f}"))