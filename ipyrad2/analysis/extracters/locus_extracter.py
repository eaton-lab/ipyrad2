#!/usr/bin/env python

"""Extract/subset sequences from HDF5 database and write loci to individual files.

Command
-------
$ ipyrad2 lex -d ... -w ... -o ... -O phy

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
CMD: ipyrad2 lex -d ... -o ... ...
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

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from typing import Dict

from loguru import logger

from ...utils.exceptions import IPyradError
from .window_extracter import NEXHEADER, WindowExtracter


class LocusExtracter:
    """Extract whole loci to individual files or one concatenated alignment."""
    def __init__(
        self,
         **kwargs: Dict[str, int | float | str | Path | None]
    ):

        # store params
        self.data = kwargs["data"]
        self.name = kwargs["name"]
        self.outdir = Path(kwargs["outdir"]).expanduser().absolute()
        self.out_format = kwargs["out_format"]
        self.force = kwargs["force"]
        self.concatenate = bool(kwargs.pop("concatenate", False))
        self.random_seed = kwargs.pop("random_seed", None)
        if self.random_seed is not None and self.random_seed < 0:
            raise IPyradError("Random seed must be a non-negative integer.")
        # Wex doesn't care about nloci or length so pop them
        self.nloci = kwargs.pop("nloci")
        self.min_length = kwargs.pop("min_length", kwargs.pop("length", None))
        if self.min_length is None:
            raise IPyradError("Internal error: lex requires a minimum locus length.")
        if self.min_length < 1:
            raise IPyradError("Minimum locus length must be at least 1.")

        # Pass the rest of the args into a wex and retrieve the phymap
        self.wex = WindowExtracter(**kwargs)
        self.wex._get_phymap_windows()
        self.wex._get_phymap()
        self.phymap = self.wex.phymap.reset_index(drop=True)

        # Create the outdir
        self.outdir.mkdir(exist_ok=True)
        self.loci = None
        self.eligible_loci_before_filtering = 0
        self.rejected_after_filtering = 0
        self.concatenated_sites = 0

        # Sample name delimiter (optional)
        self._DELIM = ""

    def _get_stats_path(self) -> Path:
        return self.outdir / f"{self.name}.stats.txt"

    def _prepare_stats_path(self) -> Path:
        stats_path = self._get_stats_path()
        self.outdir.mkdir(exist_ok=True)
        if stats_path.exists() and not self.force:
            raise IPyradError(
                f"Output file already exists: {stats_path}. Use --force to overwrite."
            )
        return stats_path

    @staticmethod
    def _parse_window_label(window: str) -> tuple[str, int, int]:
        scaff, region = window.split(":", 1)
        start, end = region.split("-", 1)
        return scaff, int(start), int(end)

    def _build_locus_stats_record(
        self,
        locus_index: int,
        locus_name: str,
        outfile,
        stats_dict: dict,
    ) -> dict:
        scaff, start, end = self._parse_window_label(locus_name)
        return {
            "locus_index": locus_index,
            "locus_name": locus_name,
            "scaffold": scaff,
            "start": start,
            "end": end,
            "outfile": outfile,
            "nsamples_before_filtering": stats_dict["nsamples_before_filtering"],
            "nsites_in_windows_before_filtering": stats_dict["nsites_in_windows_before_filtering"],
            "nvariants_in_windows_before_filtering": stats_dict["nvariants_in_windows_before_filtering"],
            "nsamples_after_filtering": stats_dict["nsamples_after_filtering"],
            "nsites_in_windows_after_filtering": stats_dict["nsites_in_windows_after_filtering"],
            "nvariants_in_windows_after_filtering": stats_dict["nvariants_in_windows_after_filtering"],
        }

    def _write_stats_summary(self, stats_records: list[dict]) -> None:
        stats_path = self._get_stats_path()
        if self.concatenate or self.out_format == "bpp":
            summary_outfile = stats_records[0]["outfile"]
        else:
            summary_outfile = "multiple"
        header = {
            "tool": "lex",
            "name": self.name,
            "infile": self.data,
            "out_format": self.out_format,
            "concatenate": self.concatenate,
            "random_seed": self.random_seed,
            "outfile": summary_outfile,
            "concatenated_sites": self.concatenated_sites if self.concatenate else 0,
            "nloci_requested": self.nloci,
            "nloci_written": len(stats_records),
            "min_length_requested": self.min_length,
            "eligible_loci_before_filtering": self.eligible_loci_before_filtering,
            "loci_rejected_after_filtering": self.rejected_after_filtering,
            "windows": self.wex.selected_windows,
            "imap": self.wex.imap,
            "min_sample_coverage_filter": self.wex.minmap,
            "max_sample_missing_filter": self.wex.max_sample_missing,
        }
        columns = [
            "locus_index",
            "locus_name",
            "scaffold",
            "start",
            "end",
            "outfile",
            "nsamples_before_filtering",
            "nsites_in_windows_before_filtering",
            "nvariants_in_windows_before_filtering",
            "nsamples_after_filtering",
            "nsites_in_windows_after_filtering",
            "nvariants_in_windows_after_filtering",
        ]

        with open(stats_path, "w", encoding="utf-8") as out:
            out.write("Summary\n")
            out.write("-------\n")
            for key, val in header.items():
                out.write(f"{key}: {val}\n")

            out.write("\n")
            out.write("Accepted loci\n")
            out.write("-------------\n")

            table_rows = [{col: str(record[col]) for col in columns} for record in stats_records]
            widths = {
                col: max(
                    len(col),
                    max((len(row[col]) for row in table_rows), default=0),
                )
                for col in columns
            }

            out.write("  ".join(col.ljust(widths[col]) for col in columns) + "\n")
            for row in table_rows:
                out.write("  ".join(row[col].ljust(widths[col]) for col in columns) + "\n")
        logger.info(f"wrote stats/log to: {stats_path}")

    @staticmethod
    def _windows_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
        return not (end1 < start2 or end2 < start1)

    def _get_candidate_loci(self) -> pd.DataFrame:
        """Return shuffled whole-locus candidates that overlap selected windows."""
        records = []
        for rowidx, row in self.phymap.iterrows():
            scaff_idx = int(row["scaff"])
            row_start = int(row["pos0"])
            row_end = int(row["pos1"])
            raw_length = row_end - row_start + 1
            if raw_length < self.min_length:
                continue

            selected_windows = self.wex.phymap_windows.get(scaff_idx, [])
            if not any(
                self._windows_overlap(row_start, row_end, win_start, win_end)
                for win_start, win_end in selected_windows
            ):
                continue

            records.append(
                {
                    "chrom": rowidx,
                    "startpos": 0,
                    "endpos": raw_length,
                    "raw_length": raw_length,
                }
            )

        if not records:
            raise IPyradError(
                "No loci met the minimum length requirement before filtering. "
                "Try reducing --min-length or selecting different windows."
            )

        loci = pd.DataFrame(records)
        order = np.random.default_rng(self.random_seed).permutation(len(loci))
        return loci.iloc[order].reset_index(drop=True)

    def _get_concatenated_output_path(self) -> Path:
        suffix = {"phy": "phy", "nex": "nex", "bpp": "bpp"}[self.out_format]
        return self.outdir / f"{self.name}.{suffix}"

    def _get_individual_output_stem(self, locus_name: str) -> str:
        scaff, start, end = self._parse_window_label(locus_name)
        return f"{self.name}.{scaff}_{start}-{end}"

    @staticmethod
    def _as_sequence_bytes(sequence: bytes | bytearray | np.ndarray) -> bytes:
        if isinstance(sequence, np.ndarray):
            return sequence.tobytes()
        return bytes(sequence)

    @classmethod
    def _format_phylip_alignment(
        cls,
        names: list[str],
        sequences: list[bytes | bytearray | np.ndarray],
        *,
        prefix: str = "",
        bpp_format: bool = False,
    ) -> str:
        sequence_bytes = [cls._as_sequence_bytes(seq) for seq in sequences]
        nsites = len(sequence_bytes[0])
        longname = max(len(name) for name in names)
        padded_names = [name.ljust(longname + 5) for name in names]
        rows = [
            f"{prefix}{padded_names[idx]} {sequence_bytes[idx].decode('utf-8')}"
            for idx in range(len(names))
        ]
        separator = "\n" if bpp_format else ""
        return f"{len(names)} {nsites}\n{separator}{chr(10).join(rows)}\n"

    @classmethod
    def _format_nexus_alignment(
        cls,
        names: list[str],
        sequences: list[bytes | bytearray | np.ndarray],
    ) -> str:
        sequence_bytes = [cls._as_sequence_bytes(seq) for seq in sequences]
        nsites = len(sequence_bytes[0])
        longname = max(len(name) for name in names)
        padded_names = [name.ljust(longname + 5) for name in names]
        lines = [NEXHEADER.format(len(names), nsites)]
        for block in range(0, nsites, 100):
            stop = min(block + 100, nsites)
            for idx, name in enumerate(padded_names):
                seq = sequence_bytes[idx][block:stop].decode("utf-8")
                lines.append(f"  {name}{seq}\n")
            lines.append("\n")
        lines.append("  ;\nend;")
        return "".join(lines)

    def _write_concatenated_alignment(
        self,
        names: list[str],
        sequences: list[bytearray],
    ) -> None:
        if self.out_format == "nex":
            contents = self._format_nexus_alignment(names, sequences)
        else:
            contents = self._format_phylip_alignment(
                names,
                sequences,
                prefix=self._DELIM if self.out_format == "bpp" else "",
                bpp_format=self.out_format == "bpp",
            )

        if self.wex.stdout:
            sys.stdout.write(contents)
        else:
            with open(self.outfile, "w", encoding="utf-8") as out:
                out.write(contents)

    def _run(self, postfix: str = None):
        self._get_loci()
        self._write_loci(postfix)


    def _get_loci(self) -> None:
        """Get whole-locus candidates that meet the pre-filter length threshold."""
        logger.info("Selecting whole-locus candidates for lex.")
        self.loci = self._get_candidate_loci()
        self.eligible_loci_before_filtering = len(self.loci)


    def _write_loci(self, postfix: str = None) -> None:
        if self.loci is None:
            msg = "No loci selected, run _get_loci() first"
            logger.info(msg)
            raise IPyradError(msg)

        if not self.eligible_loci_before_filtering:
            self.eligible_loci_before_filtering = len(self.loci)

        locus_data: list[str] = []
        stats_records: list[dict] = []
        concat_buffers = {name: bytearray() for name in self.wex.snames}
        concat_present: set[str] = set()
        self.rejected_after_filtering = 0
        self.concatenated_sites = 0
        self._prepare_stats_path()

        shared_outfile: Path | str | None = None
        if self.concatenate:
            self.outfile = self._get_concatenated_output_path()
            shared_outfile = "STDOUT" if self.wex.stdout else self.outfile
            if not self.wex.stdout and self.outfile.exists() and not self.force:
                raise IPyradError(
                    f"Output file already exists: {self.outfile}. Use --force to overwrite."
                )
            if self.out_format == "bpp":
                logger.warning(
                    "BPP is intended for multi-locus analyses; --concatenate writes "
                    "one locus and is likely unsuitable for BPP analysis."
                )
        elif self.out_format == "bpp":
            fpost = f"-{postfix}" if postfix else ""
            self.outfile = self.outdir / f"{self.name}{fpost}.phy"
            if self.outfile.exists() and not self.force:
                raise IPyradError(
                    f"Output file already exists: {self.outfile}. Use --force to overwrite."
                )

        for _, locus in self.loci.iterrows():
            if len(stats_records) >= self.nloci:
                break
            rowidx = int(locus["chrom"])
            row = self.phymap.iloc[rowidx]
            scaff_idx = int(row["scaff"])
            scaff_name = str(self.wex.scaffold_table.iloc[scaff_idx]["scaffold_name"])
            startpos = int(locus["startpos"])
            endpos = int(locus["endpos"])
            window_start = int(row["pos0"]) + startpos
            window_end = int(row["pos0"]) + endpos - 1
            window = [f"{scaff_name}:{window_start}-{window_end}"]
            self.wex.windows = window
            locus_label = window[0]
            self.wex.name = self._get_individual_output_stem(locus_label)
            locus_index = len(stats_records) + 1

            if self.concatenate:
                fnames, fseqarr, stats_dict = self.wex._get_filtered_alignment_data(
                    shared_outfile
                )
                block_length = stats_dict["nsites_in_windows_after_filtering"]
                if block_length < self.min_length:
                    self.rejected_after_filtering += 1
                    continue

                seq_by_name = {
                    name: fseqarr[idx].tobytes() for idx, name in enumerate(fnames)
                }
                missing_block = b"N" * block_length
                for name, buffer in concat_buffers.items():
                    buffer.extend(seq_by_name.get(name, missing_block))
                concat_present.update(fnames)
                self.concatenated_sites += block_length
                stats_records.append(
                    self._build_locus_stats_record(
                        locus_index,
                        locus_label,
                        shared_outfile,
                        stats_dict,
                    )
                )
                continue

            if self.out_format == "phy":
                alignment, stats_dict = self.wex._write_to_phy(
                    write_stats=False,
                    return_alignment=True,
                    return_stats=True,
                )
                if stats_dict["nsites_in_windows_after_filtering"] < self.min_length:
                    self.rejected_after_filtering += 1
                    continue
                if self.wex.stdout:
                    logger.debug("wrote alignment to stdout")
                    sys.stdout.write(alignment)
                else:
                    outfile = Path(stats_dict["outfile"])
                    with open(outfile, "w", encoding="utf-8") as out:
                        out.write(alignment.rstrip("\n"))
                    logger.info(
                        "wrote alignment ({}, {}) to: {}",
                        stats_dict["nsamples_after_filtering"],
                        stats_dict["nsites_in_windows_after_filtering"],
                        outfile,
                    )
                stats_records.append(
                    self._build_locus_stats_record(
                        locus_index,
                        locus_label,
                        stats_dict["outfile"],
                        stats_dict,
                    )
                )
            elif self.out_format == "nex":
                alignment, stats_dict = self.wex._write_to_nex(
                    write_stats=False,
                    return_alignment=True,
                    return_stats=True,
                )
                if stats_dict["nsites_in_windows_after_filtering"] < self.min_length:
                    self.rejected_after_filtering += 1
                    continue
                if self.wex.stdout:
                    logger.debug("wrote alignment to stdout")
                    sys.stdout.write(alignment)
                else:
                    outfile = Path(stats_dict["outfile"])
                    with open(outfile, "w", encoding="utf-8") as out:
                        out.write(alignment)
                    logger.info(
                        "wrote alignment ({}, {}) to: {}",
                        stats_dict["nsamples_after_filtering"],
                        stats_dict["nsites_in_windows_after_filtering"],
                        outfile,
                    )
                stats_records.append(
                    self._build_locus_stats_record(
                        locus_index,
                        locus_label,
                        stats_dict["outfile"],
                        stats_dict,
                    )
                )
            elif self.out_format == "bpp":
                locus_alignment, stats_dict = self.wex._write_to_phy(
                    write_stats=False,
                    prefix=self._DELIM,
                    bpp_format=True,
                    return_locus=True,
                    return_stats=True,
                )
                if stats_dict["nsites_in_windows_after_filtering"] < self.min_length:
                    self.rejected_after_filtering += 1
                    continue
                locus_data.append(locus_alignment)
                bpp_outfile = "STDOUT" if self.wex.stdout else self.outfile
                stats_records.append(
                    self._build_locus_stats_record(
                        locus_index,
                        locus_label,
                        bpp_outfile,
                        stats_dict,
                    )
                )
            else:
                raise IPyradError(f"Unrecognized output format: {self.out_format}")

        if not stats_records:
            raise IPyradError(
                "No loci passed the minimum length requirement after filtering. "
                "Try reducing --min-length or relaxing the locus filters."
            )

        if self.concatenate:
            concat_names = [name for name in self.wex.snames if name in concat_present]
            concat_sequences = [concat_buffers[name] for name in concat_names]
            self._write_concatenated_alignment(concat_names, concat_sequences)
            destination = "stdout" if self.wex.stdout else self.outfile
            logger.info(
                "wrote concatenated alignment ({}, {}) from {} loci to: {}",
                len(concat_names),
                self.concatenated_sites,
                len(stats_records),
                destination,
            )
        elif self.out_format == "bpp":
            contents = "\n".join(locus_data)
            if self.wex.stdout:
                sys.stdout.write(contents)
                logger.info(f"wrote {len(locus_data)} loci to stdout")
            else:
                with open(self.outfile, "w", encoding="utf-8") as outfile:
                    outfile.write(contents)
                logger.info(f"wrote {len(locus_data)} loci to: {self.outfile}")

        if len(stats_records) < self.nloci:
            logger.warning(
                "Requested {} loci, but only {} met the minimum length requirement "
                "before and after filtering.",
                self.nloci,
                len(stats_records),
            )

        self._write_stats_summary(stats_records)


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
        _ = kwargs.pop("min_length", kwargs.pop("length", None))
        _ = kwargs.pop("random_seed", None)
        _ = kwargs.pop("concatenate", None)
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
