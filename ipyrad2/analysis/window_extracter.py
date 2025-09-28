#!/usr/bin/env python

"""Extract/subset sequences from HDF5 database and write to a supermatrix.

Command
-------
$ ipyrad wex -d ... -w ... -o ... -m ... -f phy

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
from loguru import logger

NEXHEADER = """#nexus
begin data;
  dimensions ntax={ntax} nchar={nchar};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""


class WindowExtracter:
    """Tool to extract sequences from one or more loci and write to a
    concatenated sequence file in phylip or nexus format.
    """
    def __init__(
        self,
        data: str,
        name: str,
        outdir: Path | str,
        windows: str | List[str],
        min_sample_coverage: int | float,
        max_sample_missing: float,
        exclude: List[str] | None,
        imap: Dict[str, List[str]] | None,
        minmap: Dict[str, int | float] | None,
        stdout: bool,
        force: bool,
    ):
        # rmincov must be float
        assert windows, "must select one or more windows."

        # store params
        self.data = data
        self.name = name
        self.outdir = Path(outdir).expanduser().absolute()
        self.windows = [windows] if isinstance(windows, str) else list(windows)
        self.exclude = set(exclude if exclude else [])
        self.imap = imap
        self.minmap = minmap
        self.min_sample_coverage = min_sample_coverage
        self.max_sample_missing = min(1.0, max(0, max_sample_missing))
        self.stdout = stdout

        # data parsed from h5
        self.scaffold_table: pd.DataFrame = None
        self.snames: List[str] = None
        self.sidxs: List[str] = None
        self.pnames: Dict[str, str] = None
        self.phymap: pd.DataFrame = None
        self.phymap_windows: Dict[int, List[Tuple[int, int]]] = None
        self._imap: Dict[str, List[str]] = {}
        self._minmap: Dict[str, int] = {}

        # ...
        self.stats: pd.DataFrame = None

        # fills: snames, sidxs, scaffold_table
        self._get_scaffold_table()
        self._get_snames_and_sidxs_subset()
        self._get_imap_minmap()

        # run commands
    def _run(self):
        self._get_phymap_windows()
        self._get_phymap()
        self._get_seqarr()
        return self._filter_seqarr()

    def _get_imap_minmap(self):
        """Set _imap and _minmap for seqarr filtering."""
        # if no imap was entered then group all samples into one group
        # and use the global mincov as the min coverage of that group.
        if not self.imap:
            self._imap = {'all': self.snames}
            self._minmap = {'all': int(self.min_sample_coverage)}
            logger.debug(f"sample coverage minmap = {self._minmap}")

        # if imap was provided, then (1) check the names; (2) apply a
        # min value to each group from minmap; or (3) raise errors.
        else:
            assert self.minmap is not None, "must provide a minmap when using imap."
            assert set(self.minmap) == set(self.imap), "imap and minmap keys must match"
            self._imap = self.imap.copy()
            self._minmap = {}
            for key in self._imap:
                self._minmap[key] = self.minmap[key]
            logger.debug(f"sample coverage minmap = {self._minmap}")

    def _get_scaffold_table(self) -> None:
        """Store table with scaffold names and lengths in the order they are stored in H5."""
        with h5py.File(self.data, 'r') as io5:
            scaff_names = io5.attrs["scaffold_names"]
            scaff_lengths = io5.attrs["scaffold_lengths"]
            self.scaffold_table = pd.DataFrame(
                columns=["scaffold_name", "scaffold_length"],
                data={"scaffold_name": scaff_names, "scaffold_length": scaff_lengths},
            )

    def _get_snames_and_sidxs_subset(self) -> None:
        with h5py.File(self.data, 'r') as io5:
            # get sample names and get them as padded names
            snames = io5.attrs["names"]

            # auto-update exclude from imap difference
            if self.imap:  # is not None:
                imapset = set(itertools.chain(*self.imap.values()))
                self.exclude.update(set(snames).difference(imapset))
                logger.debug(
                    "dropping samples that are either not in the imap dict, "
                    f"or are in the exclude list: {self.exclude}")

            # filter to only the included samples, store their new indices (sidxs)
            self.sidxs = [i for (i, j) in enumerate(snames) if j not in self.exclude]
            self.snames = [j for (i, j) in enumerate(snames) if i in self.sidxs]

    def _get_phymap_windows(self) -> None:
        """Check each window for a matching scaffold name, and position within its bounds."""
        windows: Dict[str, List[Tuple(int, int)]] = {}

        # set names in index for easy fetching
        t = self.scaffold_table.set_index("scaffold_name")

        # iterate over user-entered windows
        for window in self.windows:

            # sub-scaffold window
            if ":" in window:
                scaff, region = window.split(":")
                assert region.count("-") == 1, f"malformatted window '{window}'. Must be {{scaff}} or {{scaff}}:{{start}}-{{end}}"
                start, end = [int(i) for i in region.split("-")]
                if scaff not in windows:
                    windows[scaff] = [(start, end)]
                else:
                    # check for overlap with other windows
                    for (s, e) in windows[scaff].items():
                        if (start < e) & (end > start):
                            raise RuntimeError(f"windows cannot be overlapping. {window} & {windows}")
                    windows[scaff] = [(start, end)]

            # full scaffold window
            else:
                mask = t.index.str.fullmatch(pat=window, na=False)
                scaffs = t.index[mask].values
                for scaff in scaffs:
                    if scaff not in windows:
                        length = int(t.loc[scaff, "scaffold_length"])
                        windows[scaff] = [(0, length)]
                    else:
                        raise RuntimeError(f"windows cannot be overlapping. {window} & {windows}")

        # log to INFO and DEBUG
        logger.debug(f"windows: {windows}")
        logger.info(f"selected {sum(len(i) for i in windows.values())} windows from {len(windows)} scaffolds")

        # store as dict mapping {scaff_index: window, ...}
        scaff_names = t.index.tolist()
        self.phymap_windows = {scaff_names.index(i): j for (i, j) in windows.items()}

    def _get_phymap(self) -> None:
        """Load the phymap for selecting windows from the seqs array."""
        with h5py.File(self.data, 'r') as io5:
            colnames = io5.attrs["columns"]
            mask = np.isin(io5["phymap"][:, 0], list(self.phymap_windows))
            phymap = pd.DataFrame(io5["phymap"][mask], columns=colnames)
        self.phymap = phymap

    def _get_seqarr(self) -> None:
        """Fill .seqarr with data from .phymap_windows for samples in .sidx."""
        phy_windows = []
        for scaff_idx, windows in self.phymap_windows.items():
            for (start, end) in windows:

                # get subset rows of phymap containing the window
                smap = self.phymap[self.phymap["scaff"] == scaff_idx]
                mask1 = smap.pos1 >= start
                mask2 = smap.pos0 <= end
                mask = mask1 & mask2
                block = smap.loc[mask, :]

                # skip if block is empty
                if not block.size:
                    continue

                # get how far past pos0 the window start is
                wmin_offset = max(0, start - int(block.iloc[0, 3]))
                # get phy start as phy0 + offset
                wmin = int(block.iloc[0, 1]) + wmin_offset

                # get how far back from pos1 the window end is
                wmax_offset = max(0, int(block.iloc[-1, 4]) - end)
                # get phy end as phy1 - offset
                wmax = int(block.iloc[-1, 2]) - wmax_offset

                # store phy start, end indices
                phy_windows.append((wmin, wmax))

        # if no windows then raise error
        if not phy_windows:
            raise ValueError("Selected windows contain zero data in the assembly. Try larger/different windows.")

        # extract sequences
        with h5py.File(self.data, 'r') as io5:
            lengths = [i[1] - i[0] for i in phy_windows]
            nsites = sum(lengths)
            seqarr = np.zeros((len(self.sidxs), nsites), dtype=np.uint8)
            x = 0
            for wlen, (start, end) in zip(lengths, phy_windows):
                seqarr[:, x:x + wlen] = io5["phy"][self.sidxs, start:end]
                x += wlen
        logger.debug(f"Extracted {nsites} sites from {len(phy_windows)} windows.")
        self.seqarr = seqarr

    def _filter_seqarr(self) -> Tuple[np.ndarray, List[str]]:
        """..."""
        # create a copy and convert - to N
        seqs = self.seqarr.copy()
        seqs[seqs == 45] = 78

        # create and apply mask for sites (columns) that fail minmap filter
        mask = np.zeros(seqs.shape[1], dtype=np.bool_)
        for pop in self._imap:
            pop_snames = self._imap[pop]
            pop_mincov = self._minmap[pop]
            pop_sidxs = [self.snames.index(i) for i in pop_snames]
            pop_arr = seqs[pop_sidxs, :]
            mask += np.sum(pop_arr != 78, axis=0) <= pop_mincov
        seqs = seqs[:, np.invert(mask)]
        if not seqs.size:
            raise ValueError("Selected windows contain zero data after filtering for coverage.")

        # get list of ordered names remaining in the array
        fnames = []
        row_missing = np.sum(seqs == 78, axis=1) / seqs.shape[1]
        mask = row_missing <= self.max_sample_missing
        for sidx in np.where(mask)[0]:
            fnames.append(self.snames[sidx])
        # todo: log.debug the dropped samples
        if not fnames:
            raise ValueError("No samples passed max_sample_missing filter.")
        return fnames, seqs

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
            "imap": self._imap,
            "min_sample_coverage_filter": self._minmap,
            "max_sample_missing_filter": self.max_sample_missing,
        }
        with open(stats_file, "w") as out:
            for key, val in stats_dict.items():
                out.write(f"{key}\t{val}\n")
        logger.info(f"wrote stats/log to: {stats_file}")


def count_snps(arr):
    """Count variants to report in the stats for an alignment."""
    m = np.ma.masked_equal(arr, 78)
    multi_cols = (np.ma.ptp(m, axis=0) > 0).filled(False)
    return int(np.sum(multi_cols))


def run_window_extracter(**kwargs):
    """command line wrapper for window-extracter.

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
        tool = WindowExtracter(**kwargs)
        tool.scaffold_table.to_csv(sys.stdout, sep="\t")
        sys.exit(0)

    tool = WindowExtracter(**kwargs)
    tool._write_to_phy()
    sys.exit(0)


if __name__ == "__main__":

    h5 = Path("/tmp/OUT_klmnop/assembly.seqs.hdf5")
    assert h5.exists(), "h5 doesn't exist"

    with h5py.File(h5, 'r') as io5:
        pass

    tool = WindowExtracter(
        data=h5,
        name='test',
        outdir=Path("/tmp/WEX"),
        windows=r"Chr[2]",
        min_sample_coverage=4,
        max_sample_missing=1.0,
        exclude=[],
        imap=None,
        minmap=None,
        stdout=True,
        force=True,
    )
    tool._write_to_phy()
    # print(tool.scaffold_table)
    # arr, stats = tool.run(return_data=True)
    # print(stats.T)
    # print(arr)

