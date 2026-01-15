#!/usr/bin/env python

"utility functions for the analysis tools"

import datetime
import glob
import time
import sys
import os

import numpy as np
import pandas as pd

from loguru import logger

from ..utils.exceptions import IPyradError


def parse_version_string(version):
    return tuple(map(int, version.split('.')))


def subsample_snps(snpsmap, seed):
    "Subsample snps, one per locus, using snpsmap"
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    sidxs = np.unique(snpsmap[:, 0])
    subs = np.zeros(sidxs.size, dtype=np.int64)
    idx = 0
    for sidx in sidxs:
        sites = snpsmap[snpsmap[:, 0] == sidx, 1]
        site = rng.choice(sites)
        subs[idx] = site
        idx += 1
    return subs


def subsample_loci(snpsmap, seed):
    """
    Return SNPs from re-sampled loci (shape = (nsample, ...can change)
    """
    np.random.seed(seed)

    # the number of unique loci with SNPs in this subset
    lidxs = np.unique(snpsmap[:, 0])

    # resample w/ replacement N loci
    lsample = np.random.choice(lidxs, len(lidxs))

    # the size of array to fill
    size = 0
    for lidx in lsample:
        size += snpsmap[snpsmap[:, 0] == lidx].shape[0]

    # fill with data
    subs = np.zeros(size, dtype=np.int64)
    idx = 0
    for lidx in lsample:
        block = snpsmap[snpsmap[:, 0] == lidx, 1]
        subs[idx: idx + block.size] = block
        idx += block.size
    return len(lidxs), subs


def get_spans(maparr, spans):
    """ 
    Get span distance for each locus in original seqarray. This
    is used to create re-sampled arrays in each bootstrap to sample
    unlinked SNPs from. Used on snpsphy or str or ...
    """
    start = 0
    end = 0
    for idx in prange(1, spans.shape[0] + 1):
        lines = maparr[maparr[:, 0] == idx]
        if lines.size:
            end = lines[:, 3].max()
            spans[idx - 1] = [start, end]
        else: 
            spans[idx - 1] = [end, end]
        start = spans[idx - 1, 1]

    # drop rows with no span (invariant loci)
    spans = spans[spans[:, 0] != spans[:, 1]]
    return spans


def count_snps(seqarr):
    nsnps = 0
    for site in range(seqarr.shape[1]):
        # make new array
        catg = np.zeros(4, dtype=np.int16)

        ncol = seqarr[:, site]
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
        # get second most common site
        catg.sort()

        # if invariant e.g., [0, 0, 0, 9], then nothing (" ")
        if catg[2] > 1:
            nsnps += 1
    return nsnps



class Params(object):
    "A dict-like object for storing params values with a custom repr"
    def __init__(self):
        self._i = 0

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __iter__(self):
        return self

    def __next__(self):
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if self._i > len(keys) - 1:
            self._i = 0
            raise StopIteration
        else:
            self._i += 1
            return keys[self._i - 1]
    next = __next__  # Python 2        

    def update(self, dict):
        self.__dict__.update(dict)


    def __repr__(self):
        "return simple representation of dict with ~ shortened for paths"
        _repr = ""
        keys = [i for i in sorted(self.__dict__.keys()) if i != "_i"]
        if keys:
            _printstr = "{:<" + str(2 + max([len(i) for i in keys])) + "} {:<20}\n"
            for key in keys:
                _val = str(self[key]).replace(os.path.expanduser("~"), "~")
                _repr += _printstr.format(key, _val)
        return _repr
