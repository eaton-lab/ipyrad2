

from pathlib import Path
import h5py
import numpy as np


def main(database):

    with h5py.File(database, "r", rdcc_nbytes=512*1024*1024, rdcc_nslots=2_000_003) as h5:
        # All SNPs (streaming)
        # G = h5["genos"][:]                     # shape (nsamples, nsnps, 3)

        # Window of SNPs (e.g., first 100k)
        # w = slice(0, 100_000)
        # Gw = h5["genos"][:, w, :]
        # M  = h5["snpsmap"][w, :]
        # R  = h5["reference"][w]

        # write an alignment
        print(h5.attrs["names"])
        print(h5['genos'].shape)
        arr = h5['genos'][:, :50, 2]
        print(arr.shape)
        for i in range(arr.shape[0]):
            print(arr[i, :].tobytes().decode())
        # print(R)

        # identify bi-allelic SNPS


    # with h5py.File(H5, 'r') as io5:
    #     print(io5.attrs.keys())
    #     print(io5.keys())

    #     genos = io5["genos"]

    #     # get snp aligment
    #     arr = genos[:, :, 2].T

    #     # get heterozygous sites


    #     # get bi-allelic sites
    #     genos


        # hetero_mask = genos[:, 0, ...] != genos[:, 1, ...]
        # print(hetero_mask)
        # genos[:, hetero_mask]


if __name__ == "__main__":

    DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    DATA = DIR / "assembly.hdf5"
    main(DATA)


    # DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    # DIR = Path("/tmp/")
    # DATA = DIR / "TEST.hdf5"


