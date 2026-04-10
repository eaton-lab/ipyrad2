#!/usr/bin/env python

"""Extract/subset sequences from HDF5 database and write to a supermatrix.

Command
-------
$ ipyrad2 analysis wex -d ... -w ... -o ... -O phy

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
CMD: ipyrad2 analysis wex -d ... -o ... ...
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

from typing import Dict, List, Tuple
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from loguru import logger

from ...utils.exceptions import IPyradError
from .sequence_common import build_sequence_imap_minmap
from .sequence_common import load_sequence_scaffold_table
from .sequence_common import normalize_sequence_population_inputs
from .sequence_common import resolve_sequence_sample_subset


NEXHEADER = """#nexus
begin data;
  dimensions ntax={} nchar={};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""
REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"


class WindowExtracter:
    """Tool to extract sequences from one or more loci and write to a
    concatenated sequence file in phylip or nexus format.
    """
    def __init__(
        self,
        data: str,
        name: str,
        outdir: Path | str,
        out_format: str,
        windows: str | List[str],
        min_sample_coverage: int | float,
        max_sample_missing: float,
        exclude: List[str] | None = None,
        include_reference: bool = False,
        imap: Dict[str, List[str]] | None = None,
        minmap: Dict[str, int | float] | None = None,
        stdout: bool = False,
        force: bool = False,
    ):
        # store params
        imap, minmap = normalize_sequence_population_inputs(imap, minmap)
        self.data = data
        self.name = name
        self.outdir = Path(outdir).expanduser().absolute()
        self.out_format = out_format
        self.windows = [] if windows is None else [windows] if isinstance(windows, str) else list(windows)
        self.exclude = set(exclude if exclude else [])
        self.include_reference = include_reference
        self.min_sample_coverage = min_sample_coverage
        self.max_sample_missing = min(1.0, max(0, max_sample_missing))
        self.stdout = stdout
        self.force = force

        # data parsed from h5
        self.scaffold_table: pd.DataFrame = None
        self.snames: List[str] = None
        self.sidxs: List[str] = None
        self.pnames: Dict[str, str] = None
        self.phymap: pd.DataFrame = None
        self.phymap_windows: Dict[int, List[Tuple[int, int]]] = None
        self.imap: Dict[str, List[str]] = {}
        self.minmap: Dict[str, int] = {}
        self.selected_windows: List[str] = []

        # fills: snames, sidxs, scaffold_table
        self.scaffold_table = load_sequence_scaffold_table(self.data)
        self.snames, self.sidxs, self.exclude, imap = resolve_sequence_sample_subset(
            self.data,
            exclude=self.exclude,
            include_reference=self.include_reference,
            imap=imap,
        )
        self.imap, self.minmap = build_sequence_imap_minmap(
            self.snames,
            min_sample_coverage=self.min_sample_coverage,
            imap=imap,
            minmap=minmap,
        )

        # run commands
    def _run(self):
        # First two are fast
        self._get_phymap_windows()
        self._get_phymap()
        # This call is slow as it is accessing the full hdf5 data
        self._get_seqarr()
        return self._filter_seqarr()

    def _get_scaffold_table(self) -> None:
        """Store table with scaffold names and lengths in the order they are stored in H5."""
        self.scaffold_table = load_sequence_scaffold_table(self.data)

    def _get_snames_and_sidxs_subset(self, imap) -> None:
        self.snames, self.sidxs, self.exclude, self.imap = resolve_sequence_sample_subset(
            self.data,
            exclude=self.exclude,
            include_reference=self.include_reference,
            imap=imap,
        )

    def _parse_imap_minmap_inputs(self, imap, minmap):
        """Normalize imap/minmap inputs from dicts or files."""
        return normalize_sequence_population_inputs(imap, minmap)

    def _get_imap_minmap(self, imap, minmap):
        """Set _imap and _minmap for seqarr filtering."""
        self.imap, self.minmap = build_sequence_imap_minmap(
            self.snames,
            min_sample_coverage=self.min_sample_coverage,
            imap=imap,
            minmap=minmap,
        )

    @staticmethod
    def _windows_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
        """Return True if 1-based inclusive windows overlap."""
        return not (end1 < start2 or end2 < start1)

    def _add_window(
        self,
        windows: Dict[str, List[Tuple[int, int]]],
        selected_windows: List[str],
        scaff: str,
        start: int,
        end: int,
        source: str,
    ) -> None:
        """Validate and store a 1-based inclusive window."""
        if start < 1 or end < start:
            raise IPyradError(
                f"Malformed window '{source}'. Windows must use valid positive coordinates."
            )

        existing = windows.setdefault(scaff, [])
        for existing_start, existing_end in existing:
            if self._windows_overlap(start, end, existing_start, existing_end):
                raise IPyradError(
                    f"windows cannot overlap. {source} overlaps "
                    f"{scaff}:{existing_start}-{existing_end}"
                )
        existing.append((start, end))
        selected_windows.append(f"{scaff}:{start}-{end}")

    def _get_phymap_windows(self) -> None:
        """Check each window for a matching scaffold name, and position within its bounds."""
        windows: Dict[str, List[Tuple[int, int]]] = {}
        selected_windows: List[str] = []

        # set names in index for easy fetching
        t = self.scaffold_table.set_index("scaffold_name")

        if not self.windows:
            logger.info(
                "No windows specified; selecting the full length of all scaffolds. "
                "Use -w to subset scaffold windows and -P to view scaffold names."
            )
            self.windows = [r".*"]

        # Load windows from bed file if they are passed in this way
        if len(self.windows) == 1:
            bedfile = Path(self.windows[0])
            if bedfile.exists():
                logger.info(f"Loading windows from bed file: '{bedfile}'")
                for scaff, start, end in self._get_windows_from_bed(bedfile):
                    if scaff not in t.index:
                        raise IPyradError(
                            f"'{scaff}' from {bedfile} does not match to any scaffold names. Check with '-P'."
                        )
                    self._add_window(
                        windows,
                        selected_windows,
                        scaff,
                        start,
                        end,
                        f"{scaff}:{start}-{end}",
                    )
                self.selected_windows = selected_windows
                logger.debug(f"windows: {windows}")
                nwindows = sum(len(i) for i in windows.values())
                ws = 's' if nwindows > 1 else ''
                ss = 's' if len(windows) > 1 else ''
                logger.info(f"selected {nwindows} window{ws} from {len(windows)} scaffold{ss}")
                scaff_names = t.index.tolist()
                scaff_to_idx = {name: idx for idx, name in enumerate(scaff_names)}
                self.phymap_windows = {scaff_to_idx[i]: j for i, j in windows.items()}
                return

            logger.debug("Loading windows from command line arguments")

        # iterate over user-entered windows
        for window in self.windows:

            # sub-scaffold window
            if ":" in window:
                scaff, region = window.split(":", 1)
                mask = t.index.str.fullmatch(pat=scaff, na=False)
                scaffs = t.index[mask].tolist()
                if not scaffs:
                    raise IPyradError(
                        f"No scaffold names match '{window.split(':')[0]}'. Use -P to view scaffold names."
                    )
                if len(scaffs) > 1:
                    raise IPyradError("Cannot use regex with ':'. List windows separately: -w Chr1:1-1000 Chr2:1-1000")
                if region.count("-") != 1:
                    raise IPyradError(f"malformatted window '{window}'. Must be {{scaff}} or {{scaff}}:{{start}}-{{end}}")
                start, end = [int(i) for i in region.split("-")]
                self._add_window(windows, selected_windows, scaffs[0], start, end, window)

            # full scaffold window
            else:
                mask = t.index.str.fullmatch(pat=window, na=False)
                scaffs = t.index[mask].values
                if not scaffs.size:
                    raise IPyradError(f"'{window}' does not match to any scaffold names. Check with '-P'.")
                for scaff in scaffs:
                    if scaff not in windows:
                        length = int(t.loc[scaff, "scaffold_length"])
                        self._add_window(windows, selected_windows, scaff, 1, length, window)
                    else:
                        raise IPyradError(f"windows cannot overlap. {window} & {windows}")

        # log to INFO and DEBUG
        self.selected_windows = selected_windows
        logger.debug(f"windows: {windows}")
        nwindows = sum(len(i) for i in windows.values())
        ws = 's' if nwindows > 1 else ''
        ss = 's' if len(windows) > 1 else ''
        logger.info(f"selected {nwindows} window{ws} from {len(windows)} scaffold{ss}")

        # store as dict mapping {scaff_index: window, ...}
        scaff_names = t.index.tolist()
        scaff_to_idx = {name: idx for idx, name in enumerate(scaff_names)}
        self.phymap_windows = {
            scaff_to_idx[i]: j
            for i, j in windows.items()
        }


    def _get_windows_from_bed(self, bedfile: Path) -> List[Tuple[str, int, int]]:
        """Read windows from a BED file.

        BED uses 0-based, half-open coordinates. These are converted to the
        extracter's 1-based inclusive region semantics.
        """
        windows: List[Tuple[str, int, int]] = []
        with open(bedfile, "r", encoding="utf-8") as infile:
            for lineno, line in enumerate(infile, start=1):
                # Ignore comments and blank lines
                if line.startswith("#") or line.strip() == "":
                    continue

                chrom, start, end, *rest = line.rstrip("\t\n").split()
                start0 = int(start)
                end0 = int(end)
                if start0 < 0 or end0 <= start0:
                    raise IPyradError(
                        f"Malformed BED window at line {lineno}: {bedfile}"
                    )
                windows.append((chrom, start0 + 1, end0))
        return windows


    def _get_phymap(self) -> None:
        """Load the phymap for selecting windows from the seqs array."""
        with h5py.File(self.data, 'r') as io5:
            phymap = io5["phymap"]
            colnames = phymap.attrs["columns"]
            mask = np.isin(phymap[:, 0], list(self.phymap_windows))
            phymap = pd.DataFrame(phymap[mask], columns=colnames)
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
            raise IPyradError("Selected windows contain zero data in the assembly. Try larger/different windows.")
        logger.debug(phy_windows)
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
        for pop in self.imap:
            pop_snames = self.imap[pop]
            pop_mincov = self.minmap[pop]
            pop_sidxs = [self.snames.index(i) for i in pop_snames]
            pop_arr = seqs[pop_sidxs, :]
            mask += np.sum(pop_arr != 78, axis=0) < pop_mincov
        seqs = seqs[:, np.invert(mask)]
        if not seqs.size:
            raise IPyradError("Selected windows contain zero data after filtering for coverage.")

        # get list of ordered names remaining in the array
        fnames = []
        row_missing = np.sum(seqs == 78, axis=1) / seqs.shape[1]
        mask = row_missing <= self.max_sample_missing
        for sidx in np.where(mask)[0]:
            fnames.append(self.snames[sidx])
        # todo: log.debug the dropped samples
        if not fnames:
            raise IPyradError("No samples passed max_sample_missing filter.")
        return fnames, seqs[mask, :]

    def _get_output_path(self, suffix: str) -> Path:
        return self.outdir / f"{self.name}.{suffix}"

    def _get_stats_path(self) -> Path:
        return self.outdir / f"{self.name}.stats.tsv"

    def _prepare_output_paths(
        self,
        suffix: str | None,
        *,
        write_stats: bool,
        return_locus: bool = False,
    ) -> Path | str:
        paths: List[Path] = []

        if suffix and not self.stdout and not return_locus:
            paths.append(self._get_output_path(suffix))
        if write_stats:
            paths.append(self._get_stats_path())

        if paths:
            self.outdir.mkdir(exist_ok=True)
            if not self.force:
                existing = next((path for path in paths if path.exists()), None)
                if existing is not None:
                    raise IPyradError(
                        f"Output file already exists: {existing}. Use --force to overwrite."
                    )

        if self.stdout:
            return "STDOUT"
        if return_locus:
            return "RETURN"
        if suffix is None:
            raise IPyradError("Internal error: missing output suffix.")
        return self._get_output_path(suffix)

    def _build_stats_dict(self, fnames, fseqarr, outfile):
        """Build stats for the extracted windows without writing them."""
        return {
            "nsamples_before_filtering": len(self.snames),
            "nsites_in_windows_before_filtering": self.seqarr.shape[1],
            "nvariants_in_windows_before_filtering": count_snps(self.seqarr),
            "nsamples_after_filtering": len(fnames),
            "nsites_in_windows_after_filtering": fseqarr.shape[1],
            "nvariants_in_windows_after_filtering": count_snps(fseqarr),
            "infile": self.data,
            "outfile": outfile,
            "windows": self.selected_windows,
            "imap": self.imap,
            "min_sample_coverage_filter": self.minmap,
            "max_sample_missing_filter": self.max_sample_missing,
        }

    def _write_stats_dict(self, stats_dict, stats_file: Path | None = None) -> None:
        """Write a precomputed stats dictionary to disk."""
        stats_path = self._get_stats_path() if stats_file is None else Path(stats_file)
        self.outdir.mkdir(exist_ok=True)
        with open(stats_path, "w", encoding="utf-8") as out:
            for key, val in stats_dict.items():
                out.write(f"{key}\t{val}\n")
        logger.info(f"wrote stats/log to: {stats_path}")

    def _write_to_phy(self, 
                      write_stats: bool = True,
                      prefix: str = None,
                      bpp_format: bool = False,
                      return_locus: bool = False,
                      return_alignment: bool = False,
                      return_stats: bool = False):
        """Writes the .seqarr matrix as a string to .outfile."""
        # get the filtered alignment
        fnames, fseqarr = self._run()
        outfile = self._prepare_output_paths(
            "phy",
            write_stats=write_stats,
            return_locus=return_locus,
        )

        # get padded names
        longname = max(len(i) for i in fnames)
        pnames = [i.ljust(longname + 5) for i in fnames]

        # build phy
        phy = []
        prefix = prefix if prefix else ""
        for idx, _ in enumerate(fnames):
            seq = fseqarr[idx].tobytes().decode("utf-8")
            phy.append(f"{prefix}{pnames[idx]} {seq}")

        # write to temp file
        ntaxa = len(fnames)
        nsites = fseqarr.shape[1]

        bpp_sep = "\n" if bpp_format else ""
        phy_text = "\n".join(phy)
        alignment = f"{ntaxa} {nsites}\n{bpp_sep}{phy_text}\n"
        stats_dict = self._build_stats_dict(fnames, fseqarr, outfile)

        if return_alignment:
            if return_stats:
                return alignment, stats_dict
            return alignment

        # write to stdout
        if return_locus:
            pass
        elif self.stdout:
            logger.debug("wrote alignment to stdout")
            sys.stdout.write(alignment)
        else:
            with open(outfile, 'w', encoding="utf-8") as out:
                out.write(alignment.rstrip("\n"))
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")

        if write_stats:
            self._write_stats_dict(stats_dict)

        if return_locus and return_stats:
            return alignment, stats_dict
        if return_locus:
            return alignment
        if return_stats:
            return stats_dict

    def _write_to_nex(
        self,
        write_stats: bool = True,
        return_alignment: bool = False,
        return_stats: bool = False,
    ):
        """Writes concatenated alignment to nex format..."""
        # get the filtered alignment
        fnames, fseqarr = self._run()
        outfile = self._prepare_output_paths("nex", write_stats=write_stats)

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
            stop = min(block + 100, fseqarr.shape[1])
            for idx, name in enumerate(pnames):
                seq = fseqarr[idx, block:stop].tobytes().decode()
                lines.append(f"  {name}{seq}\n")
            lines.append("\n")
        lines.append("  ;\nend;")
        alignment = "".join(lines)
        stats_dict = self._build_stats_dict(fnames, fseqarr, outfile)

        if return_alignment:
            if return_stats:
                return alignment, stats_dict
            return alignment

        # write to stdout
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            sys.stdout.write(alignment)
        else:
            with open(outfile, 'w', encoding="utf-8") as out:
                out.write(alignment)
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")
        if write_stats:
            self._write_stats_dict(stats_dict)
        if return_stats:
            return stats_dict

    def _write_to_fa(self, write_stats: bool = True, return_stats: bool = False):
        """Write the extracted alignment as FASTA."""
        fnames, fseqarr = self._run()
        outfile = self._prepare_output_paths("fa", write_stats=write_stats)

        records = []
        for idx, name in enumerate(fnames):
            seq = fseqarr[idx].tobytes().decode("utf-8")
            records.append(f">{name}\n{seq}")

        contents = "\n".join(records) + "\n"
        stats_dict = self._build_stats_dict(fnames, fseqarr, outfile)
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            sys.stdout.write(contents)
        else:
            with open(outfile, "w", encoding="utf-8") as out:
                out.write(contents)
            logger.info(
                f"wrote alignment ({len(fnames)}, {fseqarr.shape[1]}) to: {outfile}"
            )

        if write_stats:
            self._write_stats_dict(stats_dict)
        if return_stats:
            return stats_dict

    def _write_stats(self, fnames, fseqarr, outfile):
        """Write stats for the extracted windows."""
        self._write_stats_dict(self._build_stats_dict(fnames, fseqarr, outfile))


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
    out_format: str
        Format to write the alignments phy (default), nex, or fa.
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

    tool = WindowExtracter(**kwargs)

    if request_table:
        tool.scaffold_table.to_csv(sys.stdout, sep="\t")
        sys.exit(0)

    if tool.out_format == "phy":
        tool._write_to_phy()
    elif tool.out_format == "nex":
        tool._write_to_nex()
    elif tool.out_format == "fa":
        tool._write_to_fa()
    else:
        raise IPyradError(f"Unrecognized output format: {tool.out_format}")

    sys.exit(0)


if __name__ == "__main__":

    h5 = Path("/tmp/OUT_klmnop/assembly.seqs.hdf5")
    h5 = Path("/home/deren/Documents/ipyrad-tests/OUT/assembly.seqs.hdf5")
    h5 = Path("/home/deren/Documents/ipyrad-tests/Ped2_OUT/assembly.hdf5")
    assert h5.exists(), "h5 doesn't exist"

    with h5py.File(h5, 'r') as io5:
        print(len(io5.attrs['names']))
        print(io5["phymap"][:])
        print(io5["phy"].shape)

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
