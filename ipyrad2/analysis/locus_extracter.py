#!/usr/bin/env python

"""Extract/subset sequences from HDF5 database and write loci to individual files.

Note that genome coordinates are 1-based, closed (inclusive): both
start and end are included. This is a general standard in concordance
with samtools. If a user requests scaff:100-250 the extracted region
is 151 bp long. The wex tool handles this internally and will return
the 151 bp region by slicing from the phymap using 0-based indices.

Note that each row of the h5 phymap represent a delimited locus.

Command
-------
$ ipyrad lex -d ... -w ... -o ... -m ... -f phy

Output file example
-------------------
10 100
A    NNNNNATCCGAGA...
B    NNNNNNNNNNNNN...
C    CCAGGATCCGAAA...
D    CCAGGATCCGAAA...
...

Stats file example
------------------
CMD: ipyrad wex -d ... -o ... ...
windows: Chr1:X-Y Chr1:A-B ...
populations: A B C
min_sample_coverage: A=1 B=2 C=3
max_sample_missing: 1.0
nsamples_before_filtering: 29
nsites_in_windows_before_filtering: 1000
nvariant_sites_in_windows_before_filtering: 100
nsamples_after_filtering: 29
nsites_in_windows_after_filtering: 300
nvariant_sites_in_windows_after_filtering: 20
outfile: alignment.phy
"""

import h5py
import itertools
import numpy as np
import pandas as pd
import random
import sys
import tempfile
from loguru import logger
from pathlib import Path
from typing import List, Dict, Tuple
from ipyrad2.utils.exceptions import IPyradError
from .window_extracter import WindowExtracter
from ..utils.parallel import run_pipeline


BIN = Path(sys.prefix) / "bin"
BIN_BED = str(BIN / "bedtools")


class LocusExtracter:
    """Tool to extract sequences from one or more loci and write to a
    concatenated sequence file in phylip or nexus format.
    """
    def __init__(
        self,
         **kwargs: Dict[str, int | float | str | Path | None]
    ):

        # store params
        self.data = kwargs["data"]
        self.name = kwargs["name"]
        self.outdir = Path(kwargs["outdir"]).expanduser().absolute()
        self.out_format = kwargs["out_format"]
        # Wex doesn't care about nloci or length so pop them
        self.nloci = kwargs.pop("nloci")
        self.length = kwargs.pop("length")
        # If no windows selected, then set to a regex that will match everything
        # in the scaffold_table
        if kwargs["windows"] == None:
            kwargs["windows"] = [r".*"]
            logger.debug("lex: No windows specified. Sampling from full seq array.")

        if kwargs["exclude"] == None:
            kwargs["exclude"] = ["assembly_reference_sequence"]

        # Pass the rest of the args into a wex and retrieve the phymap
        logger.info(kwargs)
        self.wex = WindowExtracter(**kwargs)
        self.wex._get_phymap_windows()
        self.wex._get_phymap()
        self.phymap_windows = self.wex.phymap_windows
        self.phymap = self.wex.phymap
        self._validate_length()

        # Create the outdir
        self.outdir.mkdir(exist_ok=True)
        self.loci = None


    def _validate_length(self):
        """Ensure the requested length doesn't exceed the length of the
        longest locus, or else this is a non-sensical request."""
        max_len = max(self.phymap.apply(lambda x: x["pos1"] - x["pos0"], axis=1))
        if self.length > max_len:
            logger.info(f"Requested locus length {self.length} exceeds max locus size of the data {max_len}.")
            logger.info(f"  Forcing `length={max_len}`.")
            self.length = max_len


    def _run(self):
        self._get_loci()
        self._write_loci()


    def _get_loci(self) -> None:
        """Get random loci of specified length from the phymap.
        """
        logger.info("Entering _get_loci()")

        # Format phymap as a 'fai'-style file for bedtools random
        with tempfile.NamedTemporaryFile('w', delete=False) as fp:
            for chrom, windows in self.phymap.groupby("scaff"):
                for widx, (_, win) in enumerate(windows.iterrows()):
                    fp.write(f"{chrom}-{widx}\t{win['pos1'] - win['pos0']}\n")
            fp.close()

            # Get max chrom length to constrain samtools random. If this is RAD
            # data and you ask for a `-l` value that is longer than the longest
            # rad locus then samtools will spin forever.
            g = pd.read_csv(fp.name, header=None, names=["Chrom", "Length"], sep="\t")
            maxlen = max(g["Length"])

        length = self.length
        if maxlen < length:
            logger.info(f"Requested length ({self.length}) > max locus length of the data ({maxlen}). Constraining loci to {maxlen}bp")
            length = maxlen

        cmd1 = [
            BIN_BED, "random",
            "-n", str(self.nloci),
            "-l", str(length),
            "-g", fp.name
        ]
        # run the command in subprocess
        logger.debug(f"CMD: {' '.join(cmd1)}")
        rc, out, err = run_pipeline([cmd1])

        # The list comprehension drops the last (blank) line, and selects only
        # the first 3 columns of the return (the remaining columns aren't useful).
        self.loci = pd.DataFrame([x.split("\t")[:3] for x in out.decode().split("\n")[:-1]],
            columns=["chrom", "startpos", "endpos"])


    def _write_loci(self) -> None:

        if self.loci is None:
            msg = "No loci selected, run _get_loci() first"
            logger.info(msg)
            raise IPyradError(msg)

        for _, locus in self.loci.iterrows():
            # Get chrom and window id for indexing into phymap
            cidx, widx = locus["chrom"].split("-")
            wstart, wend = self.phymap_windows[int(cidx)][int(widx)]
            scaf = self.wex.scaffold_table.iloc[int(cidx)]["scaffold_name"]
            wend = int(wstart) + int(locus["endpos"])
            # Plus 1 because wex windows are 1-based inclusive
            wstart = int(wstart) + int(locus["startpos"]) + 1
            window = [f"{scaf}:{wstart}-{wend}"]
            self.wex.windows = window
            self.wex.name = window[0]
            if self.out_format == "phy":
                self.wex._write_to_phy()
            elif self.out_format == "nex":
                self.wex._write_to_nex()
            else:
                logger.error(f"Unrecognized output format: {self.out_format}")


def run_locus_extracter(**kwargs):
    """command line wrapper for locus-extracter.

    Parameters:
    -----------
    data: Path | str
        A 'seqs.hdf5' database file from ipyrad2.
    name: str
        Prefix name used for outfiles. If None it is automatically set.
    outdir: Path | str
        Dir for output files. Created if it doesn't exist.
    windows: str | List[str]:
        Subsample scaffold(s) by index number. If unsure, leave this
        empty when loading a file and then check the .scaffold_table
        to view the indices of scaffolds. Scaffolds are ordered by
        their order in the reference genome file.
    min_sample_coverage: int | float:
        Min number of individuals that must have data at a site
        for it to be included in the alignment (def=4).
    max_sample_missing: float
        Max proportion of sites that can be missing (N) in a sample.
        (def=1.0)
    exclude: List[str]
        A list of sample names to exclude from the data set. Samples
        can also be excluded by using an imap dictionary and not
        including them.
    imap: Dict
        A dictionary mapping group names (keys) to lists of sample
        names (values) to be included in the analysis. This can be
        used for 3 things: (1) to select samples to extract data for;
        (2) to filter based on sample coverage in groups (minmap);
        or (3) to use consensus_reduce=True to reduce the dataset to a
        consensus sequence for each group.
    minmap: Dict
        A dictionary mapping group names (keys) to integers or floats
        to act as a filter requiring that at least N (or N%) of samples
        in this group have data for a locus to be retained in the
        dataset. When using consensus_reduce=True the minmap applies to
        the reduced data set, i.e., it applies to the groups (keys) so
        that all values must be <= 1.
    stdout: bool
        ...
    force: bool
        ...
    """
    request_table = kwargs.pop("print_scaffold_table")
    if request_table:
        # pop args wex doesn't care about
        _ = kwargs.pop("nloci")
        _ = kwargs.pop("length")
        tool = WindowExtracter(**kwargs)
        tool.scaffold_table.to_csv(sys.stdout, sep="\t")
        sys.exit(0)

    lex = LocusExtracter(**kwargs) 
    lex._run()
    sys.exit(0)


if __name__ == "__main__":
    pass
    #with h5py.File(h5, 'r') as io5:
    #    print(io5["phymap"][:])
    #    print(io5["phy"].shape)

        # help(io5.create_dataset)


    # tool = WindowExtracter(
    #     data=h5,
    #     name='test',
    #     outdir=Path("/tmp/WEX"),
    #     windows=r"MT",
    #     min_sample_coverage=4,
    #     max_sample_missing=1.0,
    #     exclude=[],
    #     imap=None,
    #     minmap=None,
    #     stdout=True,
    #     force=True,
    # )
    # tool._write_to_phy()

    # print(tool.scaffold_table)
    # arr, stats = tool.run(return_data=True)
    # print(stats.T)
    # print(arr)

