#!/usr/bin/env python

"""
"""

from pathlib import Path
from ipyrad3.assembler.loci import build_locus_fasta_database

if __name__ == "__main__":

    REF = Path("/home/deren/Documents/tool/ipyrad2/examples/LiuLiu-genome/Pcr.genome.1.0.fasta")
    OUT = Path("/tmp/LOCI")
    snames = [
        "alaschanica-DE237",
        "axillaris-DE37",
        "bella-JJ85",
        "axillaris-JJ125",
        "bracteosa-RC2",
    ]
    kwargs = dict(
        prefix="TEST",
        snames=snames,
        reference=REF,
        outdir=OUT,
        exclude_reference=False,
        masks=[],
    )
    build_locus_fasta_database(**kwargs)
