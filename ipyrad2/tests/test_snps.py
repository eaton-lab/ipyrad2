

from pathlib import Path
import h5py
import numpy as np



    # All SNPs (streaming)
    # G = h5["genos"][:]                     # shape (nsamples, nsnps, 3)

    # Window of SNPs (e.g., first 100k)
    # w = slice(0, 100_000)
    # Gw = h5["genos"][:, w, :]
    # M  = h5["snpsmap"][w, :]
    # R  = h5["reference"][w]


def stream_process_hdf5(database: Path, ):
    """Return filtered database data.

    Streams through the hdf5 using its stored chunk sizes and keeps
    only ...
    """
    with h5py.File(database, "r", rdcc_nbytes=512*1024*1024, rdcc_nslots=2_000_003) as io5:
        snpsmap = io5["snpsmap"]
        genos = io5["genos"]
        ref = io5["reference"]

        print(snpsmap.chunks, snpsmap.shape)
        print(genos.chunks, genos.shape)
        print(ref.chunks, ref.shape)

        # the full size of the stepped dimension
        size = genos.shape[1]
        # the step size to sample along that axis
        step = genos.chunks[1]

        # iterate over chunks of size step
        for i0 in range(0, size, step):
            # slice along the step dimension
            s0 = slice(i0, min(i0 + step, size))

            # create empty arr of correct size and fill it
            genos_buff = np.empty(genos[:, s0, :].shape, dtype=genos.dtype)
            genos.read_direct(genos_buff, (slice(None), s0, slice(None)))

            snpsmap_buff = np.empty(snpsmap[s0, :].shape, dtype=snpsmap.dtype)
            snpsmap.read_direct(snpsmap_buff, (s0, slice(None)))

            ref_buff = np.empty(ref[s0].shape, ref.dtype)
            ref.read_direct(ref_buff, (s0,))

            # PROCESS CHUNK
            missing = np.sum(genos[:, :, 2] == 255)
            # print(missing)
            # ntotal = genos.size

            # get filter masks and diploid genotypes
            masks = np.zeros((genos.shape[1], 6), dtype=np.bool_)

            # mask0 is True if an indel is present
            masks[:, 0] = np.any(genos[:, :, 2] == 45, axis=0)

            # mask1 is True if >2 alleles present
            print(genos[:10, 0, :2] == 2)
            print(np.any(genos[:10, :2, :2] == 2, axis=2))
            # print(np.isin(genos[:10, :10, :2], [2, 3], axis=0))
            # print(np.any(genos[:, :, :2] > 1, axis=0))
            # masks[:, 1] = np.any(genos[:, :, :2] > 1, axis=0)

            # mask2 is True if >2 alleles present
            # masks[:, 2] = np.any(genos[:, :, :2] == 2, axis=2)
            # print(masks)




def main(database):
    stream_process_hdf5(database)


if __name__ == "__main__":

    DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    DATA = DIR / "assembly.hdf5"
    main(DATA)


    # DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    # DIR = Path("/tmp/")
    # DATA = DIR / "TEST.hdf5"


