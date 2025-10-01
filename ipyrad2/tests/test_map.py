#!/usr/bin/env python


from pathlib import Path
from ipyrad2.utils.parallel import run_pipeline

R1 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R1.trimmed.fastq.gz")
R2 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R2.trimmed.fastq.gz")
REF = Path("/home/deren/Documents/ipyrad-tests/examples/Atub-genome/AmaTu_v01_no00_renamed.fa")
OUT = Path("/tmp/")

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
