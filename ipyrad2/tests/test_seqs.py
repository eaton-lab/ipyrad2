


from pathlib import Path
import h5py
import numpy as np


def main(H5):
    with h5py.File(H5, 'r') as io5:
        print(io5.attrs.keys())
        print(io5.keys())
        phy = io5["phy"]
        phymap = io5["phymap"]

        # print(phymap[:10])
        print(phy.shape)
        # help(io5.create_dataset)


if __name__ == "__main__":

    # DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    # DATA = DIR / "assembly.seqs2.hdf5"

    DIR = Path("/home/deren/Documents/ipyrad-tests/Ama-out/")
    DATA = DIR / "assembly.hdf5"
    main(DATA)

