

from typing import Dict, Tuple
from pathlib import Path
from collections import defaultdict
import pandas as pd
from loguru import logger
from ipyrad2.utils.parallel import run_pipeline


def concat_tech_reps_into_tmpdir(popfile: Path, tmpdir: Path, fastq_dict: Dict[str, Tuple[Path, Path]]) -> Dict[str, Path]:
    """Return fastq_dict pointing to updated concat paths in tmpdir"""
    # parse the population file
    df = pd.read_csv(popfile, header=None, sep=r"\s+")

    # fill pdict and warn if names don't match any samples.
    snames = set(fastq_dict)
    pop2tups = defaultdict(list)
    pop2snames = defaultdict(list)
    for idx in df.index:
        sname, pname, *_ = df.loc[idx]
        if sname in snames:
            pop2tups[pname].append(fastq_dict.pop(sname))
            pop2snames[pname].append(sname)
        else:
            logger.warning(f"sample name '{sname}' from popfile was not found in data. Skipping.")

    # report to logger
    logger.debug(f"merging/renaming samples according to popfile {popfile}")
    maxlen = max(len(i) for i in pop2snames)
    for pname, tups in pop2tups.items():
        snames = pop2snames[pname]
        logger.debug(f"{pname}{' ' * (maxlen - len(pname))} <- {' + '.join(snames)}")
        out1 = tmpdir / f"{pname}.tmp.R1.fastq.gz"
        cmd = ["cat"] + [str(i[0]) for i in tups]
        run_pipeline([cmd], out1)

        if tups[0][1] is not None:
            out2 = tmpdir / f"{pname}.tmp.R2.fastq.gz"
            cmd = ["cat"] + [str(i[1]) for i in tups]
            run_pipeline([cmd], out2)
            fastq_dict[pname] = (out1, out2)
        else:
            fastq_dict[pname] = (out1, None)
    return fastq_dict


if __name__ == "__main__":

    DIR = Path("/home/deren/Documents/ipyrad-tests/Ped2_TRIM")
    fastq_dict = {
        "L1": (DIR / "longiflora-39348-plate_8.R1.trimmed.fastq.gz", DIR / "longiflora-39348-plate_8.R2.trimmed.fastq.gz"),
        "L2": (DIR / "longiflora-JJ60-plate_8.R1.trimmed.fastq.gz", DIR / "longiflora-JJ60-plate_8.R2.trimmed.fastq.gz"),
        "L3": (DIR / "longiflora-var-tubiformis-36071-plate_8.R1.trimmed.fastq.gz", DIR / "longiflora-var-tubiformis-36071-plate_8.R1.trimmed.fastq.gz"),
    }

    # skip/warn on header; merge A1 and A2, rename B
    popfile = Path("/tmp/popfile.tsv")
    with open(popfile, 'w') as out:
        out.write("sample\tpopulation\n")
        out.write("L1\tL\n")
        out.write("L2\tL\n")
        out.write("L3\tL33\n")

    fdict = concat_tech_reps_into_tmpdir(popfile, Path("/tmp"), fastq_dict)
    for i, j in fdict.items():
        print(i, j)