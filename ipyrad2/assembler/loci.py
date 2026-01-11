#!/usr/bin/env python

"""
1. Get sample cov beds (where they have >mindepth cov in loci.bed)
2. Get expanded mask by adding filtered (N) sites from VCF.
3. Write consensus sequences.
3. Write database (loci)

"""

from typing import List, Tuple, Dict
import re
import sys
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from loguru import logger
from ..utils.seqs import comp
from ..utils.parallel import run_pipeline

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BCF = str(BIN / "bcftools")
BIN_BED = str(BIN / "bedtools")

AMBIGARR = np.array(list(b"RSKYWM")).astype(np.uint8)


def write_sam_faidx(tmpdir: Path) -> Path:
    """Convert loci beds (0-based) to faidx 1-based (Chr:start-end).
    """
    loci_bed = tmpdir / "beds" / "loci.bed"
    fai_path = tmpdir / "loci.faidx.txt"
    awk_prog = 'BEGIN{OFS=""}{print $1,":",$2+1,"-",$3}'
    cmd = ["awk", awk_prog, str(loci_bed)]
    run_pipeline([cmd], fai_path)
    return fai_path


def get_reference_in_loci_beds(tmpdir: Path, reference: Path) -> Path:
    """Write the reference sequence as a sample to the consensus folder
    for all loci windows.
    """
    loci = tmpdir / "loci.faidx.txt"
    consensus_dir = tmpdir / "consensus_seqs"
    out_fasta = consensus_dir / "assembly_reference_sequence.consensus.fa"

    # run pipeline
    cmd = [BIN_SAM, "faidx", str(reference), "-r", str(loci)]
    run_pipeline([cmd], out_fasta)
    return out_fasta


def get_consensus(sname: str, reference: Path, tmpdir: Path, keep_insertions: bool) -> Path:
    """Write consensus sequences for one sample.

    Create FASTA for `sample_name` only over loci in `loci_bed`,
    applying variants from `vcf_gz` and masking `zero_bed` regions to N.
    """
    # step data files
    loci = tmpdir / "loci.faidx.txt"
    vcf_gz = tmpdir / "vcfs" / "variants.resolved.vcf.gz"
    consensus_dir = tmpdir / "consensus_seqs"
    consensus_dir.mkdir(parents=True, exist_ok=True)

    # sample files
    mask_bed = tmpdir / "beds" / f"{sname}.mask.bed"
    out_fasta = consensus_dir / f"{sname}.consensus.fa"

    cmd1 = [BIN_SAM, "faidx", str(reference), "-r", str(loci)]
    cmd2 = [
        BIN_BCF, "consensus",
        "-f", "-",               # read sliced FASTA from stdin
        "-s", f"{sname}",  # sample to apply
        "-M", "N",               # write N for missing genotypes
        "--mask", str(mask_bed), # mask zero/low-coverage intervals to N
        "--mask-with", "N",
        "--mark-del", "-",
        "--mark-ins", "lc" if keep_insertions else "+",
        "--regions-overlap", "1",# apply variants overlapping slice edges
        str(vcf_gz)
    ]
    cmd3 = ['tr', '-d', "'+'"]
    run_pipeline([cmd1, cmd2, cmd3], out_fasta)

    # warn if there is no data for a sample.
    if not out_fasta.stat().st_size:
        logger.warning(f"sample {sname} has no data passed filtering and should be dropped.")
    return out_fasta


# DEPRECATED
# def get_sample_masked_beds(sname: str, bam_file: Path, min_sample_depth: int, tmpdir: Path) -> Path:
#     """Write bed files to mask <min_depth or filtered sites per sample.

#     Where is the mask used?
#     -----------------------
#     This is used to write consensus loci for the sample. The BED has
#     the sites where this sample matches the reference. Sites not
#     in this file are masked, and will appear as N for this sample.
#     If the sample is variant relative to the reference then the variant
#     will be applied during consens writing.
#     """
#     bed_dir = tmpdir / "beds"
#     loci_bed = bed_dir / "loci.bed"
#     out_path = bed_dir / f"{sname}.mask.bed"

#     # write bedgraph w/ zeros (-bga) and do NOT pair fragments (-pc)
#     # otherwise it paints coverage into the inserts.
#     cmd1 = [BIN_BED, "genomecov", "-ibam", str(bam_file), "-bga"] #, "-pc"]
#     cmd2 = ["awk", "-v", f"MIN={min_sample_depth}", 'BEGIN{OFS="\t"} $4<MIN {print $1,$2,$3}']
#     cmd3 = [BIN_BED, "intersect", "-a", "-", "-b", str(loci_bed)]
#     cmd4 = [BIN_BED, "sort", "-i", "-"]
#     cmd5 = [BIN_BED, "merge", "-i", "-"]

#     run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5], out_path)
#     return out_path


def make_lowdepth_mask(sname: str, min_sample_depth: int, tmpdir: Path):
    """Build a per-bp mask of positions inside `loci_bed` where bedGraph depth < min_depth.

    Output mask contains only the A (loci) columns and is split into minimal sub-intervals
    where coverage is below threshold (including 0-coverage gaps).
    """
    bed_dir = tmpdir / "beds"
    loci_bed = bed_dir / "loci.bed"
    ref_info = tmpdir / "REF_info.txt"
    sample_bedgraph = bed_dir / f"{sname}.fragments.bedgraph"
    good_bed = bed_dir / f"{sname}.goodcov.bed"
    out_bed = bed_dir / f"{sname}.mask.bed"

    # 1) Threshold bedGraph: keep depth >= min_depth, drop depth column for set ops
    cmd1 = [
        "awk",
        f'BEGIN{{OFS="\\t"}} $4>={min_sample_depth} {{print $1,$2,$3}}',
        str(sample_bedgraph),
    ]
    cmd2 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]
    cmd3 = [BIN_BED, "merge", "-i", "-"]
    run_pipeline([cmd1, cmd2, cmd3], good_bed)

    # subtract good loci positions from all loci positions
    # Pass in ref_info so the sort order is retained
    cmd1 = [
        BIN_BED, "subtract",
        "-a", str(loci_bed),
        "-b", str(good_bed),
        "-sorted",
        "-g", str(ref_info)
    ]
    run_pipeline([cmd1], out_bed)
    return out_bed


def iter_fasta(fasta: Path):
    """Stream a multi-FASTA and yield (header, sequence).

    Parameters
    ----------
    source : path or open file-like (text mode)
        Path to FASTA (.fa/.fasta/.gz/.bz2) or an already-open text file handle.

    Yields
    ------
    (header, sequence) : tuple[str, str]
        `header` is the full header line after '>' (stripped).
        `sequence` is the concatenated sequence for that record.
    """
    # Open if a path-like was given
    fh = open(fasta, "rt", encoding="utf-8")
    header = None
    parts: list[str] = []
    try:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(parts)
                    yield header, seq.upper()
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)

        # flush last record
        if header is not None:
            seq = "".join(parts)
            yield header, seq.upper()
    finally:
        fh.close()


def iter_build_loci(fastas: List[Path]) -> Tuple[List[str], List[str]]:
    """Read all FASTA files in `indir` and group sequences by header.

    Returns
    -------
    groups : dict
        {header: [(filename, sequence), ...]} in file order.
    file_order : list[str]
        Filenames in the order they were processed.
    """
    # do not re-sort fastas here, use the input order.
    iterators = [iter_fasta(i) for i in fastas]
    names = [i.name.rsplit(".consensus.fa")[0] for i in fastas]

    while 1:
        try:
            locus = []
            for fit in iterators:
                header, seq = next(fit)
                locus.append(seq)
            yield header, names, locus
        except StopIteration:
            break


def build_locus_fasta_database(
    name: str,
    snames: List[str],
    reference: Path,
    tmpdir: Path,
    masks: List[str],
) -> Tuple[Path, Path]:
    """..."""
    # get sorted consensus fastas with reference on top
    consensus_dir = tmpdir / "consensus_seqs"
    fastas = [consensus_dir / f"{i}.consensus.fa" for i in sorted(snames)]

    # insert reference as first sample
    reference_fa = consensus_dir / "assembly_reference_sequence.consensus.fa"
    fastas = [reference_fa] + fastas

    # get names
    snames = [i.name.rsplit(".consensus.fa")[0] for i in fastas]

    # file paths
    database = tmpdir / f"{name}.database.fa"
    bed_mask = tmpdir / f"{name}.re_mask.bed"

    # restriction site sequences to be masked
    re_masks = []
    if masks:
        for mask in masks:
            re_masks.append(re.compile(mask))
            re_masks.append(re.compile(comp(mask)[::-1]))

    # iterate over loci
    beds = []
    with open(database, "w") as out_fa, open(bed_mask, "w") as out_bed:

        # iterate over loci pulled from fasta files
        lit = iter_build_loci(fastas)
        for header, names, locus in lit:

            # filter cut-sites from locus
            hits = set()
            if masks:
                for seq in locus:
                    for search in re_masks:
                        for hit in search.finditer(seq):
                            hits.add((hit.start(), hit.end()))

                # store masks to bed
                for h in hits:
                    scaff, pos = header.split(":", 1)
                    start = int(pos.split("-")[0])
                    # print("BEDMASK", scaff, start + h[0], start + h[1])
                    beds.append((scaff, start + h[0], start + h[1]))

            # build fasta
            loc = []
            for n, seq in zip(snames, locus):
                if len(seq) > seq.count("N"):
                    # mask RE sites
                    if hits:
                        seq = list(seq)
                        for h in hits:
                            seq[h[0]:h[1]] = "N" * (h[1] - h[0])
                        seq = "".join(seq)
                    # store locus
                    loc.append(f">{header} {n}\n{seq}")

            # write locus
            out_fa.write("\n".join(loc) + "\n\n")

        # write beds
        out_bed.write("\n".join(f"{i}\t{j}\t{k}" for i, j, k in beds))
    return database, bed_mask


def iter_parse_loci(database_fasta: Path):
    """Generator of (header, {names: seqs}) from database.fa"""
    ii = iter_fasta(database_fasta)
    last_scaff_pos = None
    while 1:
        try:
            locus = {}
            for fit in ii:
                header, seq = fit
                scaff_pos, sname = header.rsplit(" ", 1)
                if scaff_pos != last_scaff_pos:
                    if last_scaff_pos:
                        yield last_scaff_pos, locus
                    locus = {}
                    last_scaff_pos = scaff_pos
                locus[sname] = seq

            # flush last record
            if locus:
                yield last_scaff_pos, locus
                break

        except StopIteration:
            break


def filter_trim_locus(
    header: str,
    locus_dict: Dict[str, str],
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
):
    """Process loci from iter_parse_loci().
    """
    # parse input locus
    scaff, pos = header.split(":")
    rstart, rend = [int(i) for i in pos.split("-")]
    snames = list(locus_dict.keys())
    seqs = [list(bytes(seq, "utf-8")) for seq in locus_dict.values()]
    seqs = np.array(seqs, dtype=np.uint8)

    # dicts to fill and return
    filters = {
        "min_length": False,
        "min_samples": False,
        "max_variant_frequency": False,
        "max_shared_hetero_frequency": False,
        "max_depth_outlier": False,
    }
    stats = {
        "locus_cov": 0,       # number of samples in locus
        "variant_sites": 0,
        "variant_phylo_informative_sites": 0,
        "nsites": 0,
        "nsites_sample_cov_greater_than_1": 0,
        "nsites_sample_cov_greater_than_2": 0,
        "nsites_sample_cov_greater_than_3": 0,
        "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 0,
        "variant_site_frequency": 0,
        "variant_site_frequency_where_sample_cov_greater_than_2": 0,
        "variant_phylo_informative_site_frequency": 0,
        "variant_phylo_informative_site_frequency_where_sample_cov_greater_than_3": 0,
    }

    # apply min_samples filter --------------------------------------
    # -1 to exclude the reference sequence from sample coverage counts
    if seqs.shape[0] - 1 < min_locus_sample_coverage:
        filters["min_samples"] = True

    # apply edge trimming ---- --------------------------------------
    # get number of bases to trim from each side where sample cov < min_trim_sample_cov
    # Exclude reference sequence from site sample coverage counts (start at idx 1)
    site_sample_covs = np.sum((seqs[1:] != 78) & (seqs[1:] != 45), axis=0)
    cov_sufficient = np.where(site_sample_covs >= min_locus_trim_sample_coverage)[0]
    try:
        trim_left = int(cov_sufficient[0])
    except IndexError:
        trim_left = 0
    try:
        trim_right = seqs.shape[1] - int(cov_sufficient[-1]) - 1
    except IndexError:
        trim_right = 0
    tseqs = seqs[:, trim_left:seqs.shape[1] - trim_right]
    tsite_sample_covs = site_sample_covs[trim_left:seqs.shape[1] - trim_right]

    # get snps array. Start at row 1 to exclude reference from the stats
    snpsarr = snp_count(tseqs, rowstart=1)
    stats["variant_sites"] = int(np.sum(snpsarr > 0))
    stats["variant_phylo_informative_sites"] = int(np.sum(snpsarr == 2))

    # do not count sites where variation is not possible (sample_cov=1)
    # -1 to exclude the reference sequence from sample coverage counts
    stats["locus_cov"] = int(tseqs.shape[0]) - 1
    stats["nsites"] = int(tseqs.shape[1])
    stats["nsites_sample_cov_greater_than_1"] = int(np.sum(tsite_sample_covs > 1))
    stats["nsites_sample_cov_greater_than_2"] = int(np.sum(tsite_sample_covs > 2))
    stats["nsites_sample_cov_greater_than_3"] = int(np.sum(tsite_sample_covs > 3))
    stats["nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage"] = int(np.sum(tsite_sample_covs >= min_locus_trim_sample_coverage))

    # calculate proportion variable from sites with enough sample cov to detect pis
    if stats["nsites_sample_cov_greater_than_2"]:
        stats["variant_site_frequency_where_sample_cov_greater_than_2"] = float(stats["variant_sites"] / stats["nsites_sample_cov_greater_than_2"])

    # In addition to filtering by min-locus-length also filter by the ------------------
    # the number of non-N sites, since small loci with non-overlapping
    # filled sites could lead to almost no info. Here we just set a hard cutoff
    if min_locus_sample_coverage >= 4:
        if stats["nsites_sample_cov_greater_than_3"] < 15:  # hard-coded kind of arbitrary
            # logger.warning(f"FILTER BY MIN LENGTH: ({stats["nsites_sample_cov_greater_than_3"]}): {header} {trim_left} {trim_right}\n{seqs}\n{tseqs}\n{site_sample_covs}\n{tsite_sample_covs}\n{cov_sufficient}")
            filters["min_length"] = True
    elif min_locus_sample_coverage == 3:
        if stats["nsites_sample_cov_greater_than_2"] < 15:  # hard-coded kind of arbitrary
            filters["min_length"] = True
    elif min_locus_sample_coverage <= 2:
        if stats["nsites_sample_cov_greater_than_1"] < 15:  # hard-coded kind of arbitrary
            filters["min_length"] = True

    # filter for max proportion polymorphic sites ---------------------------------------
    if stats["variant_site_frequency_where_sample_cov_greater_than_2"] > max_locus_variant_frequency:
        filters["max_variant_proportion"] = True

    # filter for max shared het sites ----------------------------------------------------
    if tseqs.size:
        max_shared_h = max_heteros_count(tseqs)
        # -1 to exclude the reference sequence from sample coverage counts
        max_shared_h_prop = max_shared_h / (tseqs.shape[0] - 1)
        if max_shared_h_prop > max_locus_hetero_frequency:
            filters["max_shared_hetero_frequency"] = True

    # if keeping locus revise the header for trim.
    header = f"{scaff}:{rstart + trim_left}-{rend - trim_right}"
    return header, snames, tseqs, snpsarr, filters, stats


def write_loci_and_stats_files(
    snames: List[str],
    name: str,
    outdir: Path,
    tmpdir: Path,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    # read_depth_mask: np.ndarray,
):
    """
    """
    # database file is in the tmpdir inside outdir
    database = tmpdir / f"{name}.database.fa"
    loci_file = outdir / f"{name}.loci.txt"

    # add reference to stats outputs
    refname = "assembly_reference_sequence"
    snames.append(refname)

    # get name padding for loci file
    max_len = max(len(i) for i in snames) + 2
    padded = {n: n + (" " * (max_len - len(n))) for n in snames}

    # stats
    total_locus_cov = Counter()
    total_sample_cov = {i: 0 for i in snames}
    total_filters = {
        "min_length": 0,
        "min_samples": 0,
        "max_variant_frequency": 0,
        "max_shared_hetero_frequency": 0,
        "max_depth_outlier": 0,
    }
    total_stats = {
        "variant_sites": 0,
        "variant_phylo_informative_sites": 0,
        "nsites": 0,
        "nsites_sample_cov_greater_than_1": 0,
        "nsites_sample_cov_greater_than_2": 0,
        "nsites_sample_cov_greater_than_3": 0,
        "nsites_sample_cov_greater_than_or_equal_to_min_locus_trim_sample_coverage": 0,
    }

    # build
    beds = []
    loci = []
    lidx = 0    # counter of all loci
    flidx = 0   # counter of loci that passed filters
    # locus_iter = iter_parse_loci(database)
    with open(loci_file, 'w') as out:
        # for lidx, liter in enumerate(locus_iter):
        for oheader, ldict in iter_parse_loci(database):

            # # skip if masked by max depth zscore
            # if read_depth_mask[lidx]:
            #     lidx += 1
            #     total_filters["max_depth_outlier"] += 1
            #     logger.debug(f"filtered by max_depth_outlier: locus {lidx}")
            #     continue

            # apply trim and filters to locus
            args = (
                oheader,
                ldict,
                min_locus_sample_coverage,
                min_locus_trim_sample_coverage,
                min_locus_length,
                max_locus_hetero_frequency,
                max_locus_variant_frequency,
            )
            result = filter_trim_locus(*args)
            header, tnames, tseqs, snpsarr, filters, stats = result

            # update total dicts
            for key in total_filters:
                total_filters[key] += int(result[4][key])

            # tmp debugging code
            # if sum(filters.values()):
            #     logger.debug(result[4])
            #     logger.debug(header)
            #     for sname, seq in zip(tnames, tseqs):
            #         logger.debug(f"\n{padded[sname]}{bytes(seq).decode()}")

            # store for writing if locus passed filters
            if not sum(filters.values()):
                # store locus bed
                scaff, pos = header.split(":")
                pos0, pos1 = (int(i) for i in pos.split("-"))
                beds.append((scaff, pos0 - 1, pos1, tseqs.shape[0]))

                # increment sample counters
                for sname in tnames:
                    total_sample_cov[sname] += 1
                # Do not count the reference sequence in locus coverage, so
                # decrement the # of samples at a locus (n - 1)
                total_locus_cov[len(tnames) - 1] += 1
                for stat in total_stats:
                    total_stats[stat] += stats[stat]

                # build locus with snpstring
                locus = []
                for sname, seq in zip(tnames, tseqs):
                    locus.append(f"{padded[sname]}{bytes(seq).decode()}")
                snpsarr[snpsarr == 0] = 32
                snpsarr[snpsarr == 1] = 45
                snpsarr[snpsarr == 2] = 42
                snpstring = bytes(snpsarr).decode()
                locus.append(f"//{' ' * (max_len - 2)}{snpstring}|{flidx}:{header}\n")

                # store
                loci.append("\n".join(locus))
                flidx += 1
            lidx += 1

            # write in chunks
            if not flidx % 5000:
                if loci:
                    out.write("".join(loci))
                    loci = []

        # write last chunk
        if loci:
            out.write("".join(loci))

    # write locus stats -----------------------------------------------
    with open(outdir / f"{name}.stats_counts.tsv", "w") as out:
        out.write("# Locus stats and filtering (nloci tagged and excluded for each filter; one locus can hit multile filters)\n")
        out.write(f"nloci_before_filtering\t{lidx}\n")
        for key in total_filters:
            out.write(f"{key}_filter\t{total_filters[key]}\n")
        out.write(f"nloci_after_filtering\t{flidx}\n")
        for key in total_stats:
            out.write(f"{key}\t{total_stats[key]}\n")

    # write sample coverage -------------------------------------------
    # First remove the reference sequence, if it exists
    total_sample_cov.pop("assembly_reference_sequence", [])
    sample_cov = pd.DataFrame(index=["nloci"], data={i: total_sample_cov[i] for i in total_sample_cov}).T
    sample_cov.to_string(outdir / f"{name}.stats_sample_cov.txt")

    # write locus coverage stats --------------------------------------
    # TODO: add pre-filtered bed stats here.
    locus_cov = pd.DataFrame(index=['nloci'], data={i: total_locus_cov[i] for i in range(len(snames))}).T
    locus_cov.to_string(outdir / f"{name}.stats_locus_coverage.txt")

    # report stats files to user
    logger.debug(f"wrote stats files to {outdir / f'{name}.stats_*'}")

    # write a bed file with beds of loci filtered and sites trimmed to
    # be used to mask these from the final VCF
    with open(outdir / f"{name}.bed", "w") as out:
        out.write("\n".join("\t".join(map(str, i)) for i in beds))


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


def snp_count(seqs: np.ndarray, rowstart: int = 0) -> np.ndarray:
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
        ncol = seqs[rowstart:, site]

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


if __name__ == "__main__":

    import ipyrad as ip

    DIR = "Ama"
    NAME = "COL"

    DIR = "Ped"
    NAME = "BIG2"
    data = ip.load_json(f"../../{DIR}/{NAME}.json")
    data.stepdir = data.params.project_dir / f"{NAME}_outfiles"
    data.stepdir.mkdir(exist_ok=True)

    data.files.loci_bed = data.params.project_dir / f"{NAME}_clusters_within" / "beds" / "loci.bed"
    data.files.loci_vcf = data.params.project_dir / f"{NAME}_clusters_within" / "vcfs" / "loci.multi.filtered.vcf.gz"
    data.files.loci_vcf = data.params.project_dir / f"{NAME}_clusters_within" / "vcfs" / "variants.resolved.vcf.gz"
    data.files.loci_database = data.stepdir / f"{data.name}.database.fa"

    # parse loci from database
    samples = data.samples
    write_loci_and_stats_files(data, samples)

    # ii = iter_parse_loci(data.files.loci_database)

    # header, locus = next(ii)
    # print(filter_trim_locus(data, header, locus))

    # header, locus = next(ii)
    # print(filter_trim_locus(data, header, locus))

    sys.exit(0)

