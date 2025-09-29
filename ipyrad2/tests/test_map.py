#!/usr/bin/env python


import subprocess as sp
from pathlib import Path
# from ipyrad2.mapper.mapper import map_filter_sort_dedup

R1 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R1.trimmed.fastq.gz")
R2 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R2.trimmed.fastq.gz")
REF = Path("/home/deren/Documents/ipyrad-tests/examples/Atub-genome/AmaTu_v01_no00_renamed.fa")
OUT = Path("/tmp/")


import subprocess as sp
from pathlib import Path

def run_pipeline(cmds, outfile: Path):
    """
    cmds: list[list[str]] — each command must read from stdin and write to stdout ("-o -" or "-" where needed).
    outfile: final (seekable) file written by the last command.
    """
    procs = []
    prev = None

    # open final output for the LAST step
    with outfile.open("wb") as fout:
        for i, argv in enumerate(cmds):
            p = sp.Popen(
                argv,
                stdin=None if prev is None else prev.stdout,
                stdout=fout if i == len(cmds)-1 else sp.PIPE,
                stderr=sp.PIPE,
            )
            # the parent must close its copy of the previous stdout to send EOF downstream
            if prev is not None and prev.stdout is not None:
                prev.stdout.close()
            procs.append(p)
            prev = p

        # Now wait on every process and capture stderrs
        stderrs = []
        for p in procs:
            # communicate() waits and drains pipes for this proc
            # (for non-last procs, stdout was connected to next proc, so there's nothing to read here)
            _, err = p.communicate()
            stderrs.append(err.decode(errors="replace"))

        # Check return codes after waiting
        failures = [(i, cmds[i], p.returncode) for i, p in enumerate(procs) if p.returncode != 0]
        if failures:
            msg = []
            for i, cmd, rc in failures:
                msg.append(
                    f"[step {i}] cmd: {' '.join(cmd)}\nexit {rc}\nstderr:\n{stderrs[i]}"
                )
            raise RuntimeError("Pipeline failed:\n" + "\n\n".join(msg))


# Example: BWA → primaries-only → name sort → fixmate → coord sort → markdup → filter
T = "8"
TMP = "/tmp/"  # put temp on SSD/NVMe if you can

cmds = [
    ["bwa-mem2","mem","-t",T, "-K", "50000000", str(REF), str(R1), str(R2)],
    ["samtools","view","-b", "-u", "-F","0x900","-@", "1"],                                       # drop secondary+supp
    ["samtools","sort","-n","-@", "1","-m","50M","-T",f"{TMP}/name","-o","-"],              # name-sort → stdout
    ["samtools","fixmate","-m","-","-","-@","1"],                                # stdin "-", stdout to "-" by default
    ["samtools","sort","-@",T,"-m","50M","-T",f"{TMP}/coord","-o","-"],                  # coord-sort
    # add UMI if needed: ["samtools","markdup",,"-@",T,"-","-"],
    ["samtools","markdup", "-","-", "--barcode-rgx","UMI_([ACGTN]+)","-@",T],             # mark dups
    ["samtools","view","-b","-f","0x2","-q","20","-@",T, "-o", "/tmp/final.bam"],                               # filter after fixmate/markdup
]

run_pipeline(cmds, Path("final.bam"))
