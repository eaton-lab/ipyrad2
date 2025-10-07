#!/usr/bin/env python


from pathlib import Path
import numpy as np
from ipyrad2.utils.parallel import run_pipeline
from ipyrad2.utils.seqs import revcomp

R1 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R1.trimmed.fastq.gz")
R2 = Path("/home/deren/Documents/ipyrad-tests/Ama-trim-umi/SLH_AL_3065.R2.trimmed.fastq.gz")
REF = Path("/home/deren/Documents/ipyrad-tests/examples/Atub-genome/AmaTu_v01_no00_renamed.fa")
OUT = Path("/tmp/")


def write_small_genome_with_pair_merged_pair_and_clipped_pair(seed: int=123):
    """...
    ---------R1----R2----------RM--------..R1----R2..------
    """
    R1 = "AAACCCTTTGGGAAA" * 10
    R2 = "AACCTTGGAACCTTG" * 10

    RM = "ACTG" * 50

    rng = np.random.default_rng(seed)
    rng.choice("ATCG", size=500)




def map():
    pass





          #                     ...CTTGGAACTCAGTTAACTGTTCAAGTTGGGCAAGATCAAGTCGTCCCCTTAGCCCCCGCTATCTCAGGCG
R1_10 = "CCCCTATGTGTCCGGCACCCCAACGCCTTGGAACTCAGTTAACTGTTCAAGTTGGGCAAGATCAAGTCGTCCCCTTAGCCCCCGCT"
R2_10 = "CCAGCGGGAAGATGGTGCACTATCCGCAGACAATAAGTTCGGCGAGTGATACGTTCTCCATCCAAATGAAACTAACTACCTTCATTCTGGGGAGCCA"


R1_999 = "ATCGGCCTTGCAGGACTGATTGTGTATTAGTTTGCTCGTGTCATACCCTGAGGGAGCCGCCAGACGGGAGCGTCTAGCATACGAAA"
R2_999 = "GCTTTATAACATTTAGTAGGCTAAGAGATGATCCAGCTGGCAGATCTATCAGAGGAGTCGTAAAAAAGAATTCTTCGAAATAGGGAGGGGTCTCATC"



print(revcomp(R2_999))



if __name__ == "__main__":
    pass
