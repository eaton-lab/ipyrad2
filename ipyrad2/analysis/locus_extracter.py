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

from typing import List, Dict, Tuple
import sys
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import h5py
from .window_extracter import WindowExtracter
from loguru import logger
from ipyrad2.utils.exceptions import IPyradError


NEXHEADER = """#nexus
begin data;
  dimensions ntax={ntax} nchar={nchar};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""


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
        # Wex doesn't care about nloci or length so pop them
        self.nloci = kwargs.pop("nloci")
        self.length = kwargs.pop("length")
        # If no windows selected, then set to a regex that will match everything
        # in the scaffold_table
        if kwargs["windows"] == None:
            kwargs["windows"] = [r".*"]
            logger.debug("lex: No windows specified. Sampling from full seq array.")

        # Pass the rest of the args into a wex and retrieve the phymap_windows
        self.wex = WindowExtracter(**kwargs)
        self.wex._get_phymap_windows()
        self.phymap_windows = self.wex.phymap_windows


    def _run(self):
        self._get_loci()


    def _get_loci(self):
        """Get random loci of specified length from the wex"""
        logger.debug("Entering _get_loci()")


    def _write_to_phy(self) -> None:
        """Writes the .seqarr matrix as a string to .outfile."""
        # get the filtered alignment
        fnames, fseqarr = self._run()

        # get padded names
        longname = max(len(i) for i in fnames)
        pnames = [i.ljust(longname + 5) for i in fnames]

        # build phy
        phy = []
        for idx, _ in enumerate(fnames):
            seq = fseqarr[idx].tobytes().decode("utf-8")
            phy.append(f"{pnames[idx]} {seq}")

        # write to temp file
        ntaxa = len(fnames)
        nsites = fseqarr.shape[1]

        # write to stdout
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            sys.stdout.write(f"{ntaxa} {nsites}\n{'\n'.join(phy)}\n")
            outfile = "STDOUT"
        else:
            self.outdir.mkdir(exist_ok=True)
            outfile = self.outdir / f"{self.name}.phy"
            with open(outfile, 'w') as out:
                out.write(f"{ntaxa} {nsites}\n")
                out.write("\n".join(phy))
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")
        # write stats
        self._write_stats(fnames, fseqarr, outfile)

    def _write_to_nex(self, seqarr, names):
        """Writes concatenated alignment to nex format..."""
        # get the filtered alignment
        fnames, fseqarr = self._run()

        # get padded names
        longname = max(len(i) for i in fnames)
        pnames = [i.ljust(longname + 5) for i in fnames]

        # write to temp file
        ntaxa = len(fnames)
        nsites = fseqarr.shape[1]

        # write the header
        lines = []
        lines.append(NEXHEADER.format(ntaxa, nsites))

        # grab a big block of data
        for block in range(0, fseqarr.shape[1], 100):
            # store interleaved seqs 100 chars with longname+2 before
            stop = min(block + 100, seqarr.shape[1])
            for idx, name in enumerate(pnames):
                seq = fseqarr[idx, block:stop].tobytes().decode()
                lines.append(f"  {name}{seq}\n")
            lines.append("\n")
        lines.append("  ;\nend;")

        # write to stdout
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            sys.stdout.write("".join(lines))
            outfile = "STDOUT"
        else:
            self.outdir.mkdir(exist_ok=True)
            outfile = self.outdir / f"{self.name}.phy"
            with open(outfile, 'w') as out:
                out.write("".join(lines))
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")
        # write stats
        self._write_stats(fnames, fseqarr, outfile)

    def _write_stats(self, fnames, fseqarr, outfile):
        """Write stats for the extracted windows."""
        stats_file = self.outdir / f"{self.name}.stats.tsv"
        stats_dict = {
            "nsamples_before_filtering": len(self.snames),
            "nsites_in_windows_before_filtering": self.seqarr.shape[1],
            "nvariants_in_windows_before_filtering": count_snps(self.seqarr),
            "nsamples_after_filtering": len(fnames),
            "nsites_in_windows_after_filtering": fseqarr.shape[1],
            "nvariants_in_windows_afater_filtering": count_snps(fseqarr),
            "infile": self.data,
            "outfile": outfile,
            "windows": self.windows, #" ".join(self.windows),
            "imap": self.imap,
            "min_sample_coverage_filter": self.minmap,
            "max_sample_missing_filter": self.max_sample_missing,
        }
        with open(stats_file, "w") as out:
            for key, val in stats_dict.items():
                out.write(f"{key}\t{val}\n")
        logger.info(f"wrote stats/log to: {stats_file}")


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

