#!/usr/bin/env python

import numpy as np

AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)


def max_heteros_count(seqs: np.ndarray) -> int:
    """Return max number of samples with a shared polymorphism.
    """
    counts = np.zeros(seqs.shape[1], dtype=np.uint16)
    for fidx in range(seqs.shape[1]):
        subcount = 0
        for ambig in AMBIGARR:
            subcount += np.sum(seqs[:, fidx] == ambig)
        counts[fidx] = subcount
    return counts.max()


def snp_count(seqs: np.ndarray) -> np.ndarray:
    """Return the SNP array (see get_snps_array docstring).

    Parameters
    ----------
    seqs: ndarray
        A locus sequence array shape (ntaxa, nsites) in np.uint8.
    rowstart: int
        Taxon row to start on. Default if 0 (iter over all taxa),
        but when excluding the reference as counting towards
        identifying variants then the first row is skipped (the
        reference sample is always first row).
    """
    # record for every site as 0, 1, or 2, where 0 indicates the site
    # is invariant, 1=autapomorphy, and 2=synapomorphy.
    snpsarr = np.zeros(seqs.shape[1], dtype=np.uint8)

    # iterate over all loci
    for site in range(seqs.shape[1]):

        # count Cs As Ts and Gs at each site (up to 65535 sample depth)
        catg = np.zeros(4, dtype=np.uint16)

        # select the site column (potentially skipping first sample if ref.)
        ncol = seqs[:, site]

        # iterate over bases in the site column recording
        for idx in range(ncol.shape[0]):
            if ncol[idx] == 67:    # C
                catg[0] += 1
            elif ncol[idx] == 65:  # A
                catg[1] += 1
            elif ncol[idx] == 84:  # T
                catg[2] += 1
            elif ncol[idx] == 71:  # G
                catg[3] += 1
            elif ncol[idx] == 82:  # R
                catg[1] += 1       # A
                catg[3] += 1       # G
            elif ncol[idx] == 75:  # K
                catg[2] += 1       # T
                catg[3] += 1       # G
            elif ncol[idx] == 83:  # S
                catg[0] += 1       # C
                catg[3] += 1       # G
            elif ncol[idx] == 89:  # Y
                catg[0] += 1       # C
                catg[2] += 1       # T
            elif ncol[idx] == 87:  # W
                catg[1] += 1       # A
                catg[2] += 1       # T
            elif ncol[idx] == 77:  # M
                catg[0] += 1       # C
                catg[1] += 1       # A

        # sort counts so we can find second most common site.
        catg.sort()

        # if invariant      [0, 0, 0, 9] -> 0
        # if autapomorphy   [0, 0, 1, 8] -> 1
        # if synapomorphy   [0, 0, 2, 7] -> 2
        if catg[2] == 0:
            pass
        elif catg[2] == 1:
            snpsarr[site] = 1
        else:
            snpsarr[site] = 2
    return snpsarr

