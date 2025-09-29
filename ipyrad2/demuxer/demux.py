#!/usr/bin/env python

"""Some utilities used in demux.py for demultiplexing.

TODO
----
- keep paired option or only auto-detect?
- Try to speed up using 1 core for reading, 1 for processing, and N for
writing/compressing, all while restricting the size of queued reads waiting
to be written, based on this approach:
https://stackoverflow.com/questions/9770027/how-to-parse-a-large-file-taking-advantage-of-threading-in-python
"""

from typing import Dict, Tuple, List, Iterator
import sys
import itertools
from pathlib import Path
from collections import Counter
from dataclasses import dataclass, field

from loguru import logger
import pandas as pd
from pandas.errors import ParserError
from ipyrad2.utils.kmers import infer_overhang
from ipyrad2.utils.parse_names import get_name_to_fastq_dict
from ipyrad2.utils.seqs import AMBIGS, BADCHARS
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.demuxer.match import (
    BarMatchingSingleInline,
    BarMatchingCombinatorialInline,
    BarMatchingI7,
)


BASES = set("ACGTN")


@dataclass
class Demux:
    fastqs: List[Path]
    """: List of Paths to fastq files, unpaired."""
    barcodes: Path
    """: Path to the barcodes file."""
    re1: str
    """: Overhang on read1 from restriction digestion + ligation. Inferred if None."""
    re2: str
    """: Overhang on read2 from restriction digestion + ligation. Inferred if None."""
    max_mismatch: int
    """: Max number of mismatches between barcodes. Checked for conflict."""
    workers: int
    """: max number of parallel workers."""
    chunksize: int
    """: max number of reads to process between writing to disk."""
    merge_technical_replicates: bool
    """: merge replicates or append -technical-replicate-X to names."""
    outdir: Path
    """: outdir/prefix is the dir where fastqs will be written."""
    i7: bool
    """: if True then demux on i7 index instead of inline barcode(s)."""
    disable_infer_re_overhangs: bool
    """: do not infer res"""
    max_reads: int
    """: subsample only the first N reads from each file (used for testing)."""

    # attrs to be filled ----------------------------------------------
    _names_to_barcodes: Dict[str, Tuple[str, str]] = None
    """: A map of barcode strings to sample names, pre-expanded by off-by-N."""
    _filenames_to_fastqs: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)
    """: Dict mapping file short names to tuples of paired fastqs."""
    _cuts1: List[str] = None
    """: List of enzyme overhang sites to match on read1s."""
    _cuts2: List[str] = None
    """: List of enzyme overhang sites to match on read2s."""
    _barcodes_to_names: Dict[str, str] = None
    """: Dict of all acceptable barcodes (e.g., off-by-1) mapped to sample names."""
    _file_stats: Dict[str, List] = None
    """: Store stats per raw data file (pair)."""
    _sample_stats: Dict[str, int] = None
    """: Dict to store n reads per sample."""
    _technical_replicates: Dict[str, List[str]] = field(default_factory=dict)

    def __post_init__(self):
        """Run subfunctions to setup object."""
        self._get_outdir()
        self._get_barcodes_path()
        self._get_filenames_to_paired_fastqs()
        self._get_names_to_barcodes()
        self._replace_bad_name_chars()
        if not self.i7:
            if not self.disable_infer_re_overhangs:
                self._check_restriction_overhangs()
            self._get_cutters_expanded()
        self._get_barcodes_to_names_map()

    def run(self):
        """Run each file (pair) on separatre engine."""
        self._demultiplex()
        self._write_stats()
        self._merge_cleanup()

    def _get_filenames_to_paired_fastqs(self) -> None:
        self._filenames_to_fastqs = get_name_to_fastq_dict(self.fastqs)

    def _get_outdir(self) -> None:
        """Require an empty outdir to write to."""
        # get full path to the outdir
        self.outdir = Path(self.outdir).expanduser().resolve()

        # if the path exists, but is empty, that is OK.
        if self.outdir.exists():
            if any(self.outdir.iterdir()):
                raise IPyradError(
                    f"outdir '{self.outdir}' exists and contains files. "
                    "To prevent overwriting or removing data you must "
                    "manually rm this dir or change the outdir arg")
        self.outdir.mkdir(exist_ok=True)

    def _get_barcodes_path(self) -> None:
        """Get barcodes path as Path object allow for regex name."""
        bars = Path(self.barcodes)
        bpath = list(bars.parent.glob(bars.name))
        if not bpath:
            raise IPyradError(f"No barcodes file found at {self.barcodes}")
        self.barcodes = Path(bpath[0]).expanduser().resolve()

    def _get_names_to_barcodes(self) -> None:
        """Fill .names_to_barcodes dict w/ info from barcodes file.

        This logs a WARNING if technical replicates are detected to
        make sure the user is aware of how they are being handled.
        """
        # parse the tabular barcodes file on whitespace. Expects
        # there to be no header. There will be >=2 columns, >2 if
        # combinatorial barcodes.
        try:
            bardata = pd.read_csv(
                self.barcodes, header=None, sep=r"\s+",
            ).dropna()
        except ParserError as err:
            raise IPyradError(
                "Failed to parse barcodes file. Check that your sample\n"
                "names do not include spaces (invalid)"
            ) from err

        # the dataframe COULD have >3 columns, in which case we will
        # discard any extra columns to keep at most 3.
        bardata = bardata.iloc[:, :3]

        # set names on barcodes dataframe
        if bardata.shape[1] == 2:
            bardata.columns = ["sample", "barcode1"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
        else:
            bardata.columns = ["sample", "barcode1", "barcode2"]
            bardata["barcode1"] = bardata["barcode1"].str.upper()
            bardata["barcode2"] = bardata["barcode2"].str.upper()

        # check for replicate sample names in the barcodes file. These
        # are allowed, since a single sample can be sequenced multiple
        # times on the same plate with different barcodes attached,
        # representing technical replicates. THere is a demux option
        # for whether to combine tech reps, or keep as diff samples.
        if bardata['sample'].value_counts().max() > 1:
            # get duplicated names
            duplicated = (bardata['sample'].value_counts() > 1).index

            # warn that dups are present AND WILL BE merged.
            if self.merge_technical_replicates:
                logger.warning(
                    "Technical replicates are present (samples with same name "
                    "in barcodes file) and will be merged into one sample. "
                    "Stats will be reported for each replicate and for the merged sample.")

            # warn that dups are present and WILL NOT be merged.
            else:
                logger.warning(
                    "Technical replicates are present (samples with same name "
                    "in barcodes file) and will have '-technical-replicate-x' "
                    "appended to their sample names")

            # either way, relabel the samples for now, and may or may not merge later.
            for dup in duplicated:
                ridxs = bardata[bardata['sample'] == dup]
                if ridxs.shape[0] > 1:
                    for idx, index in enumerate(ridxs.index):
                        newname = f"{dup}-technical-replicate-{idx}"
                        bardata.loc[index, 'sample'] = newname

        # make sure barcodes are valid characters.
        if not all(bardata["barcode1"].apply(set("RKSYWMCATG").issuperset)):
            raise IPyradError(
                "Barcodes file contains unexpected characters in the "
                "barcode sequences suggesting it is not correctly "
                "formatted. See documentation.")

        # convert bardata to a dictionary {sample: barcode}.
        # if combinatorial barcodes are present then combine them.
        if "barcode2" in bardata.columns:
            # check that data is paired
            for fname, ftuple in self._filenames_to_fastqs.items():
                if not ftuple[1]:
                    raise IPyradError(
                        "Only paired-end reads can make use of combinatorial "
                        "barcodes. The barcode table suggests multiple barcodes "
                        "but the fastq file names suggest data are not paired."
                    )
            self._names_to_barcodes = dict(zip(
                bardata["sample"], zip(bardata["barcode1"], bardata["barcode2"])
            ))
        else:
            self._names_to_barcodes = dict(zip(
                bardata["sample"], ((i, "") for i in bardata["barcode1"])
            ))
        # report to logger
        logger.debug(f"barcodes map:\n{bardata}")

    def _replace_bad_name_chars(self) -> None:
        """Replaces bad characters in names in .names_to_barcodes."""
        names = list(self._names_to_barcodes)
        for name in names:
            if any(i in name for i in BADCHARS):
                newname = name
                for badchar in BADCHARS:
                    newname = newname.replace(badchar, "_")
                # newname = "".join([i.replace(i, "_") for i in BADCHARS])
                logger.warning(f"changing name {name} to {newname} (bad characters).")
                self._names_to_barcodes[newname] = self._names_to_barcodes.pop(name)

    def _check_restriction_overhangs(self) -> None:
        """Use kmer analysis to detect restriction overhangs in sequences."""
        read1s = [i[0] for i in self._filenames_to_fastqs.values()]
        read2s = [i[1] for i in self._filenames_to_fastqs.values()]
        max_reads = int(200_000 / len(read1s))
        infer_cut1 = infer_overhang(read1s, max_len=20, max_reads=max_reads, anchored=False)
        infer_cut2 = infer_overhang(read2s, max_len=20, max_reads=max_reads, anchored=False)

        if self.re1:
            if self.re1 != infer_cut1:
                logger.warning(
                    f"user entered {self.re1} as the restriction overhang, but kmer "
                    f"analysis suggests {infer_cut1} is the most likely restriction "
                    "overhang in R1s."
                )
            else:
                logger.info(
                    f"kmer analysis confirms {self.re1} as the restriction "
                    "overhang for R1s."
                )
        else:
            self.re1 = infer_cut1
            logger.info(
                f"kmer analysis detected {self.re1} as the restriction "
                "overhang for R1s."
            )

        if self.re2:
            if self.re2 != infer_cut2:
                logger.warning(
                    f"user entered {self.re2} as the restriction overhang, but kmer "
                    f"analysis suggests {infer_cut2} is the most likely restriction "
                    "overhang in R2s."
                )
            else:
                logger.info(
                    f"kmer analysis confirms {self.re2} as the restriction "
                    "overhang for R2s."
                )
        else:
            self.re2 = infer_cut2
            if self.re2:
                logger.info(
                    f"kmer analysis detected {self.re2} as the restriction "
                    "overhang for R2s."
                )

    def _get_cutters_expanded(self) -> None:
        """Fill `.cuts1` and `.cuts2` with ordered list of resolutions.

        Sequences will be searched for cut sites starting with the
        entered value and then proceeding to allow off-by-n matches.
        The first tested sequences will be the user entered value,
        with IUPAC resolved, followed by off-by-1 matches.
        """
        cuts1 = self.re1
        if any(i in 'RKSYWM' for i in cuts1):
            res1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cuts1]
            res2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cuts1]
            cuts1 = [res1, res2]
        else:
            cuts1 = [cuts1]
        cuts1 = cuts1 + list(set(itertools.chain(*[mutate(i) for i in cuts1])))

        cuts2 = self.re2
        if any(i in 'RKSYWM' for i in cuts2):
            res1 = [AMBIGS[i][0] if i in "RKSYWM" else i for i in cuts2]
            res2 = [AMBIGS[i][1] if i in "RKSYWM" else i for i in cuts2]
            cuts2 = [res1, res2]
        else:
            cuts2 = [cuts2]
        cuts2 = cuts2 + list(set(itertools.chain(*[mutate(i) for i in cuts2])))

        # convert all str to bytes
        self._cuts1 = [i.encode() for i in cuts1]
        self._cuts2 = [i.encode() for i in cuts2]
        # logger.info(self._cuts1)
        # logger.info(self._cuts2)

    def _get_barcodes_to_names_map(self) -> None:
        """Fills .barcodes_to_names with all acceptable barcodes: name.

        This updates the .barcodes_to_names from {str: Tuple[str,str]}
        to {str: str}.
        """
        # store perfect match to barcodes
        self._barcodes_to_names = {}

        # finished if no mismatch is allowed.
        if not self.max_mismatch:
            for name, barcode in self._names_to_barcodes.items():
                # convert tuple to string with _ separator
                barc = (
                    f"{barcode[0]}" if not barcode[1] else
                    f"{barcode[0]}_{barcode[1]}"
                )
                self._barcodes_to_names[barc.encode()] = name
            return

        # iterate over barcodes: names
        warning = False
        for name, barcode in self._names_to_barcodes.items():

            # get generators of off-by-n barcodes
            if self.max_mismatch == 1:
                gen1 = mutate(barcode[0])
                gen2 = mutate(barcode[1])
            else:
                gen1 = itertools.chain(*[(mutate(i)) for i in mutate(barcode[0])])
                gen2 = itertools.chain(*[(mutate(i)) for i in mutate(barcode[1])])
            bars1 = set(gen1)
            bars2 = set(gen2)

            # if only one barcode
            if not bars2:
                barcgen = iter(bars1)
            else:
                barcgen = (f"{i}_{j}" for (i, j) in itertools.product(bars1, bars2))

            for barc in barcgen:
                barc = barc.encode()
                if barc not in self._barcodes_to_names:
                    self._barcodes_to_names[barc] = name
                else:
                    logger.warning(
                        f"\nSample: {name} ({barc}) is within "
                        f"{self.max_mismatch} "
                        f"changes of sample ({self._barcodes_to_names[barc]}).")
                    warning = True

        if warning:
            logger.warning(
                "Ambiguous barcodes that match to multiple samples "
                "will arbitrarily be assigned to the first sample.\n"
                "If you do not like this then lower the value of "
                "max_mismatch and rerun (recommended).")

    def _demultiplex(self) -> None:
        """Send fastq tuples to barmatch() function to process in parallel."""

        # barmatching performed on each fastq (r1, r2) file pair serially.
        # Some parallelization occurs within barmatch.
        jobs = {}
        for fname, fastq_tuple in self._filenames_to_fastqs.items():
            short = tuple(i.name if i else "" for i in fastq_tuple)
            logger.info(f"processing {fname} {short}")
            jobs[fname] = barmatch(fastq_tuple, self)

        # record of stats per file
        self._file_stats = jobs

        # record of stats per sample from barmatch returned objects
        self._sample_stats = Counter()
        for _, stats in self._file_stats.items():
            for sname, hits in stats[2].items():
                # store the full name stats
                self._sample_stats[sname] += hits

                # also record stats for combined technical-replicates
                if "-technical-replicate-" in sname:
                    short_name = sname.split("-technical-replicate-")[0]
                    self._sample_stats[short_name] += hits
                    if short_name not in self._technical_replicates:
                        self._technical_replicates[short_name] = [sname]
                    else:
                        self._technical_replicates[short_name].append(sname)
                # NOTE: either the technical-replicates or the merged
                # sample will be dropped after writing the statsfile.

        # report failed samples
        for sname in self._names_to_barcodes:
            if not self._sample_stats[sname]:
                logger.warning(f"Sample {sname} has 0 reads.")
                self._sample_stats[sname] = 0

    def _write_stats(self) -> None:
        """Write to {project_dir}/`s1_demultiplex_stats.txt`.

        The stats file includes the number of reads per sample as well
        as information about demultiplexing in terms of nreads per file
        and the barcodes that were found.
        """
        # open the stats file for writing.
        stats_file = self.outdir / "demultiplexing_stats.txt"
        outfile = open(stats_file, 'w', encoding="utf-8")

        # write the per-file stats
        outfile.write("# Raw file statistics\n######################\n")
        file_df = pd.DataFrame(
            index=sorted(self._file_stats),
            columns=["total_reads", "cut_found", "bar_matched"],
        )
        for key in sorted(self._file_stats):
            stats = self._file_stats[key]
            not_cut = sum(stats[0].values())
            matched = sum(stats[1].values())
            total = not_cut + matched
            file_df.loc[key, :] = total, total - not_cut, matched
        outfile.write(file_df.to_string() + "\n\n")

        # write sample nreads stats ----------------------------------
        outfile.write("# Sample demux statistics\n######################\n")
        sample_df = pd.DataFrame(
            index=sorted(self._sample_stats),
            columns=["reads_raw"],
            data=[
                # self._sample_stats[i] for i in sorted(self._names_to_barcodes)
                self._sample_stats[i] for i in sorted(self._sample_stats)
            ],
        )
        outfile.write(sample_df.to_string() + "\n\n")
        logger.info(f"demultiplexing statistics written to {stats_file}")
        # logger.info("\n" + sample_df.to_string())

        # write verbose barcode information --------------------------
        outfile.write("# Barcode detection statistics\n######################\n")

        # record matches
        data = []
        bar_obs = Counter()
        for key in self._file_stats:
            bar_obs.update(self._file_stats[key][1])
        sorted_bar_obs = sorted(bar_obs, key=lambda x: bar_obs[x], reverse=True)

        for name in sorted(self._names_to_barcodes):
            truebar = self._names_to_barcodes[name]
            for foundbar in sorted_bar_obs:
                if name == self._barcodes_to_names[foundbar]:
                    count = bar_obs[foundbar]
                    if count:
                        if b"_" in foundbar:
                            foundbar = tuple(i.decode() for i in foundbar.split(b"_"))
                        else:
                            foundbar = (foundbar.decode(), )
                        data.append([name, truebar, foundbar, count])

        # record misses
        bad_bars = Counter()
        for key in sorted(self._file_stats):
            bad_bars.update(self._file_stats[key][0])
        bad_bar_obs = sorted(bad_bars, key=lambda x: bad_bars[x], reverse=True)
        for badbar in bad_bar_obs:
            count = bad_bars[badbar]
            if b"_" in badbar:
                badbar = tuple(i.decode() for i in badbar.split(b"_"))
            else:
                badbar = badbar.decode()
            data.append(["no_match", "", badbar, count])
        barcodes_df = pd.DataFrame(
            index=[i[0] for i in data],
            columns=["true_bar", "observed_bar", "N_records"],
            data=[i[1:] for i in data],
        )
        outfile.write(barcodes_df.to_string() + "\n")
        outfile.close()

    def _merge_cleanup(self) -> None:
        """Remove keys from _sample_stats for merging."""
        if not self.merge_technical_replicates:
            for key in self._technical_replicates.keys():
                self._sample_stats.pop(key)
        else:
            for key, value in self._technical_replicates.items():
                for rep in value:
                    self._sample_stats.pop(rep)


######################################################################
######################################################################

def mutate(barcode: str) -> Iterator[str]:
    """Mutate a sequence by 1 base (ACGT). Used for barcode mismatch."""
    for pos, _ in enumerate(barcode):
        for sub in BASES:
            newbar = list(barcode)
            newbar[pos] = sub
            yield "".join(newbar)


def barmatch(fastq_tuple, demux_obj):
    """Call .run to barmatch using a class from barmatch.py."""
    kwargs = dict(
        fastqs=fastq_tuple,
        barcodes_to_names=demux_obj._barcodes_to_names,
        cuts1=demux_obj._cuts1,
        cuts2=demux_obj._cuts2,
        merge_technical_replicates=demux_obj.merge_technical_replicates,
        outdir=demux_obj.outdir,
        workers=demux_obj.workers,
        chunksize=demux_obj.chunksize,
        max_reads=demux_obj.max_reads,
    )

    if demux_obj.i7:
        logger.info("demultiplexing on i7 index")
        barmatcher = BarMatchingI7(**kwargs)
    else:
        # TODO: maybe support other options like 2BRAD here...
        if b"_" in list(demux_obj._barcodes_to_names)[0]:
            logger.info("demultiplexing on R1+R2 inline barcodes")
            barmatcher = BarMatchingCombinatorialInline(**kwargs)
        else:
            logger.info("demultiplexing on R1 inline barcodes")
            barmatcher = BarMatchingSingleInline(**kwargs)
    try:
        barmatcher.run()

    # this is not catching...
    # except MemoryError:
    #     logger.error(
    #         "Insufficient memory.\n This can be prevented by decreasing "
    #         "the 'chunksize' parameter, which will ensure data is written "
    #         "to disk more frequently.")
    #     raise
    except KeyboardInterrupt:
        logger.warning("interrupted by user. Shutting down.")
        raise
    except Exception:
        raise
    return barmatcher.barcode_misses, barmatcher.barcode_hits, barmatcher.sample_hits


def run_demuxer(**kwargs):
    """Command-line wrapper for Demux."""
    tool = Demux(**kwargs)
    tool.run()


if __name__ == "__main__":

    import shutil
    import os
    from ipyrad2.utils.logger import set_log_level
    set_log_level("DEBUG")

    DATA = Path("/home/deren/Documents/ipyrad-tests/")
    tool = Demux(
        fastqs=DATA / "iTru*.gz",
        barcodes=DATA / "barcode*.csv",
        outdir="/tmp/DEMUX",
        max_mismatch=0,
        workers=4,
        chunksize=int(1e7),
        re1="ATCGG",
        re2="CGATCC",
        i7=False,
        merge_technical_replicates=False,
        disable_infer_re_overhangs=True,
    )
    tool.run()

    sys.exit(1)

    # tool = Demux(
    #     barcodes_path="../../tests/ipsimdata/rad_example_barcodes.txt",
    #     fastq_paths="../../tests/ipsimdata/rad_example_R1*.gz",
    #     outdir="/tmp/demux_rad_example",
    #     max_barcode_mismatch=0,
    #     re1="TGCAG",
    #     re2="",
    # )
    # tool.run()

    # tool = Demux(
    #     barcodes_path="../../tests/ipsimdata/rad_example_barcodes_techreps_badchars.txt",
    #     fastq_paths="../../tests/ipsimdata/rad_example_R1*.gz",
    #     outdir="/tmp/demux_rad_example_techreps_badchars",
    #     max_barcode_mismatch=0,
    #     re1="TGCAG",
    #     re2="",
    # )
    # tool.run()

    # tool = Demux(
    #     barcodes_path="../../tests/ipsimdata/pairgbs_wmerge_example_barcodes.txt",
    #     fastq_paths="../../tests/ipsimdata/pairgbs_wmerge_example_R*.fastq.gz",
    #     outdir="/tmp/demux_pairgbs",
    #     max_barcode_mismatch=1,
    #     re1="TGCAG",
    #     re2="TGCAG",
    # )
    # tool.run()

    # tool = Demux(
    #     barcodes_path="../../tests/ipsimdata/pairddrad_example_barcodes.txt",
    #     fastq_paths="../../tests/ipsimdata/pairddrad_example_R*.gz",
    #     outdir="/tmp/demux_pairddrad_example",
    #     max_barcode_mismatch=0,
    #     re1="TGCAG",
    #     re2="CGG",
    # )
    # tool.run()

    if os.path.exists("/home/deren/Documents/tools/ipyrad2/examples/demux_2024-8-8"):
        shutil.rmtree("/home/deren/Documents/tools/ipyrad2/examples/demux_2024-8-8")
    tool = Demux(
        barcodes="../../pedtest/barcodes-fewer-plate1.csv",
        fastqs="../../pedtest/Pedicularis_plate1_R*.fastq.gz",
        outdir="../../pedtest/demux_2024-3-16",
        max_barcode_mismatch=1,
        workers=7,
        chunksize=1e6,
        # re1="ATCGG",
        # re2="CGATCC",
    )
    tool.run()

    # COMMAND LINE TOOL EXAMPLE
    # cmd = ['ipyrad', 'demux', ']

    # if os.path.exists("/tmp/radcamp_i7"):
    #     shutil.rmtree("/tmp/radcamp_i7")
    # tool = Demux(
    #     barcodes_path="../../sandbox/radcamp/SMALL_i7_barcodes.txt",
    #     # barcodes_path="../../sandbox/radcamp/SMALL_i7_barcodes_techrep_test.txt",
    #     fastq_paths="../../sandbox/radcamp/SMALL_RAW_R*.fastq",
    #     outdir="/tmp/radcamp_i7",
    #     chunksize=10_000,
    #     max_barcode_mismatch=1,
    #     merge_technical_replicates=True,  # testing w/ alt brcodes file.
    #     i7=True,
    # )
    # tool.run()

    # # TEST i7 demux.
    # DATA = ip.Assembly("TEST_i7")
    # DATA.params.raw_fastq_path = "../../sandbox/radcamp/SMALL_RAW_R*.fastq"
    # DATA.params.barcodes_path = "../../sandbox/radcamp/SMALL_i7_barcodes.txt"
    # DATA.params.project_dir = "/tmp"
    # DATA.params.max_barcode_mismatch = 1
    # DATA.hackers.demultiplex_on_i7_tags = True

    # DATA = ip.Assembly("TEST1")
    # DATA.params.raw_fastq_path =
    # DATA.params.barcodes_path =
    # DATA.params.project_dir = "/tmp"
    # DATA.params.max_barcode_mismatch = 0
    # DATA.run('1', force=True, quiet=True)
    # print(DATA.stats)

    # DATA.params.raw_fastq_path = "../../tests/ipsimdata/pairddrad_example_*.gz"
    # DATA.params.barcodes_path = "../../tests/ipsimdata/pairddrad_example_barcodes.txt"
    # DATA.params.datatype = "pairddrad"

    # # TEST i7 demux.
    # DATA = ip.Assembly("TEST_i7")
    # DATA.params.raw_fastq_path = "../../sandbox/radcamp/SMALL_RAW_R*.fastq"
    # DATA.params.barcodes_path = "../../sandbox/radcamp/SMALL_i7_barcodes.txt"
    # DATA.params.project_dir = "/tmp"
    # DATA.params.max_barcode_mismatch = 1
    # DATA.hackers.demultiplex_on_i7_tags = True
