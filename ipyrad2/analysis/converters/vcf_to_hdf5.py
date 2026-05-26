#!/usr/bin/env python
"""
convert VCF to database format for SNP analyses"

snpsmap (ip2 specification)
-------
description: The map of SNP positions on RAD loci and genome scaffolds.
Loci, scaffolds, and positions are all stored 0-indexed.
dtype: np.uint8
shape: (nsnps, 5)
attrs["columns"]: ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
attrs["indexing"]: [0, 0, 0, 0, 0]

reference
---------
description: The ordered REF allele at every SNP position.
dtype: np.uint8
shape: (nsnps,)

sample_dp
---------
description: The ordered per-sample FORMAT/DP value at every SNP position.
dtype: np.uint32
shape: (nsamples, nsnps)

site_qual
---------
description: The ordered VCF QUAL value at every SNP position.
dtype: np.float32
shape: (nsnps,)
"""

import os
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from loguru import logger

from ...utils.progress import ProgressBar
from ...utils.exceptions import IPyradError
from ..extracters.snps_extracter import _MISSING_GENO

class VCFToHDF5(object):
    """
    Creates a temporary snps.hdf5 file conversion of the VCF file.
    For ipyrad assembled RAD seq data this will use RAD loci as the
    grouping of SNPs into linkage blocks for subsampling. If VCF is from
    some other assembly (e.g., WGS) then you can use the ld_block_size arg to
    group SNPs into linkage blocks for subsampling analyses.
    """
    def __init__(
        self,
        data,
        name="test",
        workdir="./analysis-vcf2hdf5",
        ld_block_size=10000,
        quiet=False,
        ):

        # attrs
        self.data = data
        self.name = (name if name else "test")
        self.workdir = (workdir if workdir else tempfile.gettempdir())
        self.names = []
        self.nsamples = 0
        self.nsnps = 0
        self.hlines = 0
        self.ld_block_size = ld_block_size
        self.database = ""
        self.quiet = quiet

        # check for data file
        self.database = os.path.join(self.workdir, self.name + ".hdf5")
        if not os.path.exists(self.data):
            raise IPyradError(f"VCF file not found: {self.data}")
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # vcf format info
        self.source = ""
        self.reference = ""


    def run(self,
            force=False):
        """
        Parse and convert data to HDF5 file format.
        """
        # print message
        self._print("Indexing VCF to HDF5 database file")

        if os.path.exists(self.database) and not force:
            raise IPyradError(
                f"HDF5 file already exists: {self.database}. Use --force to overwrite."
            )
        else:
            pass

        # get sample names, count header lines and nsnps
        self.get_meta()

        # init the database to fill
        self.init_database()

        # fill snps matrix
        self.build_chunked_matrix()

        # Modify the database to conform to ipy2 hdf5 standard
        self._finalize_database()

        # report on new database
        with h5py.File(self.database, 'r') as io5:
            self.nscaffolds = io5["snpsmap"][-1, 0]
            # self.nlinkagegroups = io5["snpsmap"][-1, 3]

        self._print(
            "HDF5: {} SNPs; {} linkage group"
            .format(
                self.nsnps,
                self.nscaffolds,
            )
        )
        self._print(
            "SNP database written to {}"
            .format(self.database)
        )


    def _print(self, msg):
        if not self.quiet:
            logger.info(msg)


    def get_meta(self):
        """
        Skip and count ## lines, and get names from first # line in VCF.
        """
        # store a list of chrom names
        chroms = set()

        if self.data.endswith(".gz"):
            import gzip
            infile = gzip.open(self.data)
        else:
            infile = open(self.data)

        # get data header line
        for dat in infile:

            # split on space
            try:
                data = dat.decode().split()
            except AttributeError:
                data = dat.split()

            # parse meta data lines
            if data[0][0] == "#":
                if data[0][1] != "#":
                    # parse names from header line in their order
                    self.names = data[9:]
                    self.nsamples = len(self.names)
                else:
                    # store header line count and look for format str
                    self.hlines += 1
                    if "source=" in data[0].lower():
                        self.source = data[0].split("source=")[-1].strip()
                    if "reference=" in data[0].lower():
                        self.reference = data[0].split("reference=")[-1].strip()

            # meta snps data
            else:
                self.nsnps += 1
                chroms.add(data[0])

        # close file handle
        infile.close()

        # convert chroms into a factorized list
        self._print("VCF: {} SNPs; {} scaffolds".format(self.nsnps, len(chroms)))        


    def init_database(self):
        """
        # load vcf file as a pandas dataframe in chunks.
        """
        # init the database file
        with h5py.File(self.database, 'w') as io5:
            io5.attrs["version"] = 2.0
            io5.attrs["names"] = np.array(self.names, dtype=h5py.string_dtype(encoding="utf-8"))
            io5.attrs["nsnps"] = int(self.nsnps)

            # core data sets (should SNPs be S1?)
            io5.create_dataset("genos", (self.nsnps, self.nsamples, 3), np.uint8)
            io5.create_dataset("snps", (self.nsamples, self.nsnps), np.uint8)
            io5.create_dataset("snpsmap", (self.nsnps, 5), np.uint32)
            io5.create_dataset("reference", (self.nsnps,), np.uint8)
            io5.create_dataset("sample_dp", (self.nsnps, self.nsamples), np.uint32)
            io5.create_dataset("site_qual", (self.nsnps,), np.float32)
            io5["snps"].attrs["names"] = [i.encode() for i in self.names]
            io5["genos"].attrs["names"] = [i.encode() for i in self.names]
            io5["snpsmap"].attrs["columns"] = [
                b"loc", b"loc_idx", b"loc_pos", b"scaff", b"pos",
            ]


    def _finalize_database(self):
        """ This merges the genos and snps datasets into one samps x snps x 3
        dataset, to conform to the new ipyrad2 hdf5 format."""

        # Reshape genos to align with new ip2 hdf5 format
        with h5py.File(self.database, 'a') as io5:
            io5["tmp"] = io5["genos"]  
            io5["tmp_dp"] = io5["sample_dp"]
            del io5["genos"]
            del io5["sample_dp"]
            del io5["snps"]
        with h5py.File(self.database, 'a') as io5:
            io5.create_dataset("genos", (self.nsamples, self.nsnps, 3), np.uint8)
            io5.create_dataset("sample_dp", (self.nsamples, self.nsnps), np.uint32)
            io5["genos"][:] = np.ascontiguousarray(io5["tmp"][:].transpose(1, 0, 2))
            io5["sample_dp"][:] = np.ascontiguousarray(io5["tmp_dp"][:].transpose(1, 0))
            io5.attrs["nsnps"] = int(self.nsnps)
            io5.attrs["names"] = np.array(self.names, dtype=h5py.string_dtype(encoding="utf-8"))
            del io5["tmp"]
            del io5["tmp_dp"]


    def build_chunked_matrix(self):
        """
        Fill HDF5 database with VCF data in chunks at a time.
        """

        # chunk retriever
        self.df = pd.read_csv(
            self.data, 
            sep="\t", 
            skiprows=self.hlines, 
            chunksize=int(1e5),
            index_col=False,  # prevent from interpreting int CHROM as index
        )

        # open h5
        with h5py.File(self.database, 'a') as io5:
            prog = ProgressBar(self.nsnps, 0, "converting VCF to HDF5")
            prog.finished = 0
            prog.update()
            try:
                # iterate over chunks of the file
                xx = 0
                lastchrom = "NULL"
                e0 = -1  # 0-indexed new-locus index, will advance in get_snps/lastchrom
                e1 = 0  # 0-indexed snps-per-loc index
                e2 = 0  # 0-indexed snps-per-loc position
                e3 = 0  # 0-indexed original-locus index, TODO, advancer
                e4 = 0  # 0-indexed global snps counter
                for chunkdf in self.df:

                    # get sub arrays
                    genos, snps, reference, sample_dp, site_qual = chunk_to_arrs(
                        chunkdf,
                        self.nsamples,
                    )

                    # get sub snpsmap
                    snpsmap, lastchrom = self.get_snpsmap(
                        chunkdf, lastchrom=lastchrom, e0=e0, e1=e1, e2=e2, e4=e4)

                    # store sub arrays
                    e0 = snpsmap[-1, 0].astype(int)
                    e1 = snpsmap[-1, 1].astype(int) + 1
                    e2 = snpsmap[-1, 2].astype(int) + 1
                    e4 = snpsmap[-1, 4].astype(int) + 1

                    # write to HDF5
                    io5['snps'][:, xx:xx + chunkdf.shape[0]] = snps.T
                    io5['genos'][xx:xx + chunkdf.shape[0], :, :2] = genos
                    io5['snpsmap'][xx:xx + chunkdf.shape[0], :] = snpsmap
                    io5['reference'][xx:xx + chunkdf.shape[0]] = reference
                    io5['sample_dp'][xx:xx + chunkdf.shape[0], :] = sample_dp
                    io5['site_qual'][xx:xx + chunkdf.shape[0]] = site_qual
                    xx += chunkdf.shape[0]

                    # print progress
                    prog.finished = xx
                    prog.update()

                # return with last chunk
                self.df = chunkdf
            finally:
                prog.close()

            # close h5 handle
            self._print("")

    def get_snpsmap(self, chunkdf, lastchrom, e0, e1, e2, e4):

        # convert snps back to "S1" view to enter data...
        nsnps = chunkdf.shape[0]
        snpsmap = np.zeros((nsnps, 5), dtype=np.uint32)

        # check whether locus is same as end of last chunk
        currchrom = chunkdf.iloc[0, 0]
        if currchrom != lastchrom:
            e0 += 1
            e1 = 0
            e2 = 0

        # hotfix b/c of a TypeError 4/1/2025
        e1 = np.uint32(e1)

        # snpsmap: if ipyrad denovo it's easy, and they should just use hdf5.
        if ("ipyrad" in self.source) and ("pseudo-ref" in self.reference):

            # print warning that we're not ussng ld_block_size
            if self.ld_block_size:
                self._print(
                    "\nThis appears to be a denovo assembly, "
                    "ld_block_size arg is being ignored.")

            # get locus index
            snpsmap[:, 0] = (
                chunkdf["#CHROM"].factorize()[0].astype(np.uint32) + e0)

            # get snp index counter possibly continuing from last chunk
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in chunkdf.groupby("#CHROM")])
            snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1

            # get snp pos counter possibly continuing from last chunk
            snpsmap[:, 2] = chunkdf.POS + e2
            snpsmap[:, 3] = snpsmap[:, 0]

            # get total snp counter, always continuing from any previous chunk
            snpsmap[:, 4] = range(e4, snpsmap.shape[0] + e4)


        # snpsmap: if ipyrad ref VCF the per-RAD loc info is available too
        elif "ipyrad" in self.source:

            # skip to generic vcf method if ld_block_size is set:
            if not self.ld_block_size:
                snpsmap[:, 0] = (
                    chunkdf["#CHROM"].factorize()[0].astype(np.uint32) + e0)

                # add ldx counter from last chunk
                snpsmap[:, 1] = np.concatenate(
                    [range(i[1].shape[0]) for i in chunkdf.groupby("#CHROM")])
                snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1
                snpsmap[:, 2] = chunkdf.POS + e2
                snpsmap[:, 3] = snpsmap[:, 0]
                snpsmap[:, 4] = range(e4, snpsmap.shape[0] + e4)

        else:
            # snpsmap: for other program's VCF's we need ldsize arg to chunk.
            if not self.ld_block_size:
                raise IPyradError(
                    "You must enter an ld_block_size estimate for this VCF.")

        # cut it up by block size (unless it's denovo, then skip.)
        if (self.ld_block_size) and ("pseudo-ref" not in self.reference):

            # create a BLOCK column to keep track of original chroms
            chunkdf["BLOCK"] = 0

            # block and df index counters
            original_e0 = e0
            dfidx = e4

            # iterate over existing scaffolds (e.g., could be one big chrom)
            for _, scaff in chunkdf.groupby("#CHROM"):

                # current start and end POS of this scaffold before breaking
                gpos = e2
                end = scaff.POS.max()

                # iterate to break scaffold into linkage blocks
                while 1:

                    # grab a block 
                    mask = (scaff.POS >= gpos) & (scaff.POS < gpos + self.ld_block_size)
                    block = scaff[mask]

                    # check for data and sample a SNP
                    if block.size:

                        # enter new block into dataframe
                        chunkdf.loc[dfidx:dfidx + block.shape[0], "BLOCK"] = e0
                        dfidx += block.shape[0]
                        e0 += 1

                    # advance counter
                    gpos += self.ld_block_size

                    # break on end of scaff
                    if gpos > end:
                        break

            # store it (CHROMS/BLOCKS are stored 0-indexed !!!!!!)
            snpsmap[:, 0] = chunkdf.BLOCK
            snpsmap[:, 1] = np.concatenate(
                [range(i[1].shape[0]) for i in chunkdf.groupby("BLOCK")])
            # add ldx counter from last chunk
            snpsmap[snpsmap[:, 0] == snpsmap[:, 0].min(), 1] += e1
            snpsmap[:, 2] = chunkdf.POS
            snpsmap[:, 3] = chunkdf["#CHROM"].factorize()[0] + original_e0
            snpsmap[:, 4] = range(e4, chunkdf.shape[0] + e4)
        return snpsmap, currchrom


def run_vcf_to_hdf5(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    ld_block_size: int,
    force: bool = False,
) -> Path:
    outdir = Path(outdir).expanduser().absolute()
    outdir.mkdir(parents=True, exist_ok=True)
    tool = VCFToHDF5(
        data=str(Path(data).expanduser().absolute()),
        name=name,
        workdir=str(outdir),
        ld_block_size=ld_block_size,
        quiet=False,
    )
    tool.run(force=force)
    logger.info("wrote SNP HDF5 database to {}", tool.database)
    return Path(tool.database)


def _safe_parse_qual(value) -> np.float32:
    """Return one float QUAL value or NaN when unavailable."""
    if pd.isna(value):
        return np.float32(np.nan)
    text = str(value).strip()
    if text in {"", "."}:
        return np.float32(np.nan)
    try:
        return np.float32(float(text))
    except ValueError:
        return np.float32(np.nan)


def _parse_gt_alleles(gt_field: str) -> tuple[np.uint8, np.uint8]:
    """Return diploid allele indexes from one GT token or missing sentinels."""
    if not isinstance(gt_field, str):
        return _MISSING_GENO, _MISSING_GENO
    token = gt_field.strip()
    if token in {"", ".", "./.", ".|."}:
        return _MISSING_GENO, _MISSING_GENO
    token = token.replace("|", "/")
    parts = token.split("/")
    if len(parts) != 2:
        return _MISSING_GENO, _MISSING_GENO
    alleles = []
    for part in parts:
        if part in {"", "."}:
            return _MISSING_GENO, _MISSING_GENO
        try:
            parsed = int(part)
        except ValueError:
            return _MISSING_GENO, _MISSING_GENO
        alleles.append(parsed if parsed >= 0 else _MISSING_GENO)
    return np.uint8(alleles[0]), np.uint8(alleles[1])


def _parse_dp_value(dp_field: str) -> np.uint32:
    """Return one non-negative DP value, defaulting missing or malformed to 0."""
    if not isinstance(dp_field, str):
        return np.uint32(0)
    token = dp_field.strip()
    if token in {"", "."}:
        return np.uint32(0)
    try:
        parsed = int(token)
    except ValueError:
        return np.uint32(0)
    return np.uint32(max(0, parsed))


def _extract_format_arrays(chunkdf, nsamples):
    """Return GT-derived genotype indexes and per-sample DP arrays for one chunk."""
    nsnps = int(chunkdf.shape[0])
    g0 = np.full((nsnps, nsamples), _MISSING_GENO, dtype=np.uint8)
    g1 = np.full((nsnps, nsamples), _MISSING_GENO, dtype=np.uint8)
    sample_dp = np.zeros((nsnps, nsamples), dtype=np.uint32)

    for ridx, row in enumerate(chunkdf.itertuples(index=False, name=None)):
        format_field = row[8] if len(row) > 8 else ""
        format_keys = str(format_field).split(":") if format_field else []
        try:
            gt_idx = format_keys.index("GT")
        except ValueError:
            gt_idx = -1
        try:
            dp_idx = format_keys.index("DP")
        except ValueError:
            dp_idx = -1

        for sidx, sample_field in enumerate(row[9 : 9 + nsamples]):
            parts = str(sample_field).split(":") if sample_field else []
            gt_field = parts[gt_idx] if 0 <= gt_idx < len(parts) else ""
            dp_field = parts[dp_idx] if 0 <= dp_idx < len(parts) else ""
            g0[ridx, sidx], g1[ridx, sidx] = _parse_gt_alleles(gt_field)
            sample_dp[ridx, sidx] = _parse_dp_value(dp_field)

    return g0, g1, sample_dp


def chunk_to_arrs(chunkdf, nsamples):
    """
    In development...
    Read in chunk of VCF and convert to numpy arrays
    """
    # nsnps in this chunk
    nsnps = chunkdf.shape[0]

    # base calls as int8 (0/1/2/3/255)
    ref = np.frombuffer(''.join(chunkdf.iloc[:, 3]).encode('ascii'), dtype=np.uint8)
    alts = chunkdf.iloc[:, 4].astype(bytes)
    sas = np.char.replace(alts, b",", b"")
    alts1 = np.zeros(alts.size, dtype=np.uint8)
    alts2 = np.zeros(alts.size, dtype=np.uint8)
    alts3 = np.zeros(alts.size, dtype=np.uint8)
    lens = np.array([len(i) for i in sas])
    alts1[lens == 1] = [i[0] for i in sas[lens == 1]]
    alts2[lens == 2] = [i[1] for i in sas[lens == 2]]
    alts3[lens == 3] = [i[2] for i in sas[lens == 3]]

    # genotypes and per-sample depth from row-specific FORMAT layouts.
    g0, g1, sample_dp = _extract_format_arrays(chunkdf, nsamples)
    genos = np.zeros((nsnps, nsamples, 2), dtype=np.uint8)
    genos[:, :, 0] = g0
    genos[:, :, 1] = g1

    site_qual = np.array(
        [_safe_parse_qual(value) for value in chunkdf.iloc[:, 5]],
        dtype=np.float32,
    )

    # numba func to fill
    snps = jfill_snps(nsnps, nsamples, ref, g0, g1, alts1, alts2, alts3)
    return genos, snps, ref, sample_dp, site_qual


def jfill_snps(nsnps, nsamples, ref, g0, g1, alts1, alts2, alts3):

    # fill snps
    snps = np.zeros((nsnps, nsamples), dtype=np.uint8)

    # fill snps in rows by indexing genos from ref,alt with g0,g1
    for ridx in range(snps.shape[0]):

        # get it
        tmpr = ref[ridx]
        tmp0 = g0[ridx]
        tmp1 = g1[ridx]

        # missing set to 78
        tmpsnps = snps[ridx]
        tmpsnps[tmp0 == _MISSING_GENO] = 78
        snps[ridx] = tmpsnps

        # 0/0 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 + tmp1) == 0] = tmpr
        snps[ridx] = tmpsnps

        # 1/1 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 1) & (tmp1 == 1)] = alts1[ridx]
        snps[ridx] = tmpsnps   

        # 2/2 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 2) & (tmp1 == 2)] = alts2[ridx]
        snps[ridx] = tmpsnps   

        # 3/3 put to ref allele
        tmpsnps = snps[ridx]
        tmpsnps[(tmp0 == 3) & (tmp1 == 3)] = alts3[ridx]
        snps[ridx] = tmpsnps 

    # fill ambiguity sites 
    ambs = np.where(g0 != g1)
    for idx in range(ambs[0].size):

        # row, col indices of the ambiguous site in the snps mat
        row = ambs[0][idx]
        col = ambs[1][idx]

        # get genos (0123) from the geno matrices
        a0 = g0[row, col]
        a1 = g1[row, col]        
        alls = sorted([a0, a1])

        # get the alleles (CATG) from the ref/alt matrices
        if alls[0] == 0:
            b0 = ref[row]
            if alls[1] == 1:
                b1 = alts1[row]
            elif alls[1] == 2:
                b1 = alts2[row]
            else:
                b1 = alts3[row]
        elif alls[0] == 1:
            b0 = alts1[row]
            if alls[1] == 2:
                b1 = alts2[row]
            else:
                b1 = alts3[row]
        elif alls[0] == 2:
            b0 = alts2[row]
            b1 = alts3[row]

        # convert allele tuples into an ambiguity byte
        fill = np.argmax((GETCONS[:, 2] == b0) & (GETCONS[:, 1] == b1))
        snps[row, col] = GETCONS[fill, 0]

    # return the three arrays
    return snps

# # used in write_outfiles.write_geno
TRANSFULL = {
    ('G', 'A'): "R",
    ('G', 'T'): "K",
    ('G', 'C'): "S",
    ('T', 'C'): "Y",
    ('T', 'A'): "W",
    ('C', 'A'): "M",
    ('A', 'C'): "M",
    ('A', 'T'): "W",
    ('C', 'T'): "Y",
    ('C', 'G'): "S",
    ('T', 'G'): "K",
    ('A', 'G'): "R",
}


# used in baba.py / write_outfiles..py
## with N and - masked to 255
GETCONS = np.array([
    [82, 71, 65],
    [75, 71, 84],
    [83, 71, 67],
    [89, 84, 67],
    [87, 84, 65],
    [77, 67, 65],
    [78, 255, 255],
    [45, 255, 255],
    ], dtype=np.uint8)
