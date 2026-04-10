#!/usr/bin/env python

"""Coverage-BED helpers for the active `ipyrad2 assemble` workflow."""

from __future__ import annotations

import sys
import shutil
from pathlib import Path
import numpy as np
from loguru import logger
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline
from ..utils.parallel import stream_pipeline_lines

BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")
BIN_BED = str(BIN / "bedtools")
BIN_BCF = str(BIN / "bcftools")
CALLABLE_REFERENCE_BASES = frozenset("ACGT")


def get_name_from_bam(bam_file: Path) -> str:
    """Return the sample name recorded in a BAM header."""
    cmd = [BIN_SAM, "samples", bam_file]
    _, out, _ = run_pipeline([cmd])
    return out.decode().strip().split()[0]


def samtools_index_reference(reference: Path, threads: int) -> None:
    """Index the reference FASTA with samtools if needed downstream."""
    cmd = [BIN_SAM, "faidx", reference]
    run_pipeline([cmd])


def get_reference_sort_order(reference: Path, tmpdir: Path) -> Path:
    """Write the reference scaffold order file used by bedtools sorting."""
    out_path = tmpdir / "REF_info.txt"

    fai_path = reference.with_suffix(reference.suffix + ".fai")
    if not fai_path.exists():
        samtools_index_reference(reference, 4)

    cmd = ["cut", "-f", "1,2", str(fai_path)]
    run_pipeline([cmd], out_path)
    return out_path


def sort_bed_by_reference_order(in_bed: Path, out_bed: Path, ref_info: Path) -> Path:
    """Sort a BED file to match the canonical scaffold order in REF_info.txt."""
    cmd = [BIN_BED, "sort", "-i", str(in_bed), "-g", str(ref_info)]
    run_pipeline([cmd], out_bed)
    return out_bed


def _iter_selected_fasta_records(
    reference_fasta: Path,
    contigs: set[str],
):
    """Yield `(name, sequence)` for requested FASTA records only."""
    name: str | None = None
    chunks: list[str] | None = None
    found_any = False
    with Path(reference_fasta).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and chunks is not None:
                    found_any = True
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                if not name:
                    raise IPyradError(
                        f"Reference FASTA contains an empty header: {reference_fasta}"
                    )
                chunks = [] if name in contigs else None
                continue
            if chunks is not None:
                chunks.append(line)
    if name is not None and chunks is not None:
        found_any = True
        yield name, "".join(chunks)
    elif not found_any and contigs:
        raise IPyradError(
            f"Reference FASTA contains no sequence records: {reference_fasta}"
        )


def _iter_callable_slices(sequence: str, start: int, end: int):
    """Yield contiguous callable `(start, end)` slices inside one interval."""
    run_start: int | None = None
    for pos in range(start, end):
        base = sequence[pos].upper()
        if base in CALLABLE_REFERENCE_BASES:
            if run_start is None:
                run_start = pos
            continue
        if run_start is not None:
            yield run_start, pos
            run_start = None
    if run_start is not None:
        yield run_start, end


def write_callable_regions_bed(
    regions_bed: Path,
    reference_fasta: Path,
    out_bed: Path,
) -> Path:
    """Write BED fragments limited to A/C/G/T reference runs inside regions_bed."""
    intervals: list[tuple[int, str, int, int]] = []
    by_contig: dict[str, list[tuple[int, int, int]]] = {}
    with Path(regions_bed).open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 3:
                raise IPyradError(
                    f"BED line {line_no} in {regions_bed} does not have at least 3 columns."
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise IPyradError(
                    f"BED line {line_no} in {regions_bed} has a non-integer start/end."
                ) from exc
            if start < 0 or end < start:
                raise IPyradError(
                    f"BED line {line_no} in {regions_bed} has invalid coordinates: {chrom}:{start}-{end}"
                )
            if start == end:
                continue
            order = len(intervals)
            intervals.append((order, chrom, start, end))
            by_contig.setdefault(chrom, []).append((order, start, end))

    out_bed = Path(out_bed)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    if not intervals:
        out_bed.write_text("", encoding="utf-8")
        return out_bed

    fragments_by_order: dict[int, list[tuple[int, int]]] = {
        order: [] for order, *_ in intervals
    }
    seen_contigs: set[str] = set()
    for chrom, sequence in _iter_selected_fasta_records(
        reference_fasta, set(by_contig)
    ):
        seen_contigs.add(chrom)
        seq_len = len(sequence)
        for order, start, end in by_contig[chrom]:
            if end > seq_len:
                raise IPyradError(
                    f"BED interval exceeds reference length for {chrom}: {start}-{end} > {seq_len}"
                )
            fragments_by_order[order].extend(
                _iter_callable_slices(sequence, start, end)
            )

    missing = sorted(set(by_contig) - seen_contigs)
    if missing:
        joined = ", ".join(missing[:10])
        raise IPyradError(
            f"Loci BED references scaffolds that were not found in {reference_fasta}: {joined}"
        )

    with out_bed.open("w", encoding="utf-8") as out:
        for order, chrom, _start, _end in intervals:
            for frag_start, frag_end in fragments_by_order[order]:
                out.write(f"{chrom}\t{frag_start}\t{frag_end}\n")
    return out_bed


def get_coverage_bed_graphs(
    sname: str,
    bam_file: Path,
    is_paired: bool,
    reference: Path,
    tmpdir: Path,
    min_map_q: int,
    min_sample_depth: int,
    min_merge_distance: int,
    threads: int,
) -> Path:
    r"""Produce a fragments BED (full inserts) from a coordinate-sorted BAM.

    Shell command
    -------------
    >>> $ samtools collate -u -@ 2 -O S1.sorted.bam \
    >>>   | bedtools bamtobed -bedpe -i - \
    >>>   | awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}' \
    >>>   | bedtools sort -i - \
    >>>   | bedtools genomecov -i - -g REF.info -bg > S1.bedgraph
    """
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(parents=True, exist_ok=True)
    coll_dir = tmpdir / f"{sname}.collate"
    coll_dir.mkdir(exist_ok=True)
    out_bed_count = bed_dir / f"{sname}.fragments.bedgraph"  # has counts
    out_bed_merge = bed_dir / f"{sname}.fragments.merged.bed"  # merged, no counts
    fai_path = reference.with_suffix(reference.suffix + ".fai")

    # Collate keeps read mates adjacent before converting to BED/BEDPE. The
    # resulting bedgraph is later reused for low-depth masking and depth stats.
    cmd1 = [
        BIN_SAM,
        "collate",
        "-@",
        str(min(threads, 4)),  # doesn't benefit from >4
        "-T",
        str(coll_dir / f"{sname}"),
        "-r",
        "1000000",
        "-u",
        "-O",
        str(bam_file),
    ]
    # compute bedpe table for each pair record (skipping pairs if one was filtered out), e.g.,
    # Chr1    60908424        60908533        Chr1    60908434        60908543        LH00150:341:22HGMLLT3:3:1160:35876:20034        49      -       +
    # Chr3    5915120         5915253         Chr3    5915131         5915265         LH00150:341:22HGMLLT3:3:2236:23034:9892         57      +       -
    # Chr2    109898149       109898235       Chr2    109898287       109898367       LH00150:341:22HGMLLT3:3:1287:23015:18000        54      +       -
    cmd2 = [BIN_BED, "bamtobed"]
    if is_paired:
        cmd2.append("-bedpe")
    cmd2.extend(["-i", "-"])
    # extract only the records from this table where the mean mapq passes this steps filter
    # Chr1    60908424        60908533        Chr1    60908434        60908543        LH00150:341:22HGMLLT3:3:1160:35876:20034        49      -       +
    # Chr3    5915120         5915253         Chr3    5915131         5915265         LH00150:341:22HGMLLT3:3:2236:23034:9892         57      +       -
    # Chr2    109898149       109898235       Chr2    109898287       109898367       LH00150:341:22HGMLLT3:3:1287:23015:18000        54      +       -
    cmd3 = ["awk", "-v", f"q={min_map_q}"]
    if is_paired:
        cmd3 += [r'BEGIN{OFS="\t"} ($8+0) >= q']
    else:
        cmd3 += [r'BEGIN{OFS="\t"} ($5+0) >= q']
    # check start/end values and extract only (chrom, start, end)
    # Chr4    108577598       108577715
    # Chr4    106107228       106107344
    # Chr1    45721044        45721144
    if is_paired:
        cmd4 = [
            "awk",
            r'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e}',
        ]
    else:
        cmd4 = ["awk", r'BEGIN{OFS="\t"} {print $1,$2,$3}']
    # sort beds by genome coordinates
    # Chr1    825268  825547
    # Chr1    825268  825547
    # Chr1    825268  825547
    # Chr1    833321  833418
    # Chr1    833321  833418
    # Chr1    833321  833418
    cmd5 = [BIN_BED, "sort", "-i", "-", "-faidx", str(fai_path)]
    # get coverage in bedgraph format with counts within regions
    # Chr1    833321  833418  14
    # Chr1    833419  833520  2
    # Chr1    837052  837165  4
    # Chr1    837165  837230  18
    # Chr1    837230  837240  14
    cmd6 = [BIN_BED, "genomecov", "-i", "-", "-g", str(fai_path), "-bg"]
    # filter out sites below MIN_DEPTH coverage
    # Chr1    833321  833418  14
    # Chr1    837052  837165  4
    # Chr1    837165  837230  18
    # Chr1    837230  837240  14
    cmd7 = ["awk", "-v", f"MIN={min_sample_depth}", r"$4>=MIN", "-"]
    # pipe one stream forward and save another to file. This saved file
    # has the per-site depths that will be used later for...
    cmd8 = ["tee", str(out_bed_count)]
    # merge beds within MIN_MERGE_DISTANCE of each other.
    # -d should always be >0 to merge across small low depth or alignment
    # gaps within samples. The default value of 300 should generally be
    # sufficient to group paired reads. Larger values can also merge
    # neighboring (linked) loci.
    # bedtools merge expects plain coordinate order; denovo locus ids such as
    # locus_3_8 / locus_3_16 can be valid in reference order but fail that check.
    cmd9 = ["sort", "-k1,1", "-k2,2n", "-T", str(tmpdir)]
    cmd10 = [BIN_BED, "merge", "-d", str(min_merge_distance), "-i", "-"]
    cmd11 = [BIN_BED, "sort", "-i", "-", "-g", str(fai_path)]
    # Chr1    833321  833418
    # Chr1    837052  837240
    run_pipeline(
        [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9, cmd10, cmd11],
        out_bed_merge,
    )
    shutil.rmtree(coll_dir)
    logger.debug(f"wrote bed graph for {sname}")
    return out_bed_merge


def get_across_sample_loci_bed(
    snames: list[str],
    min_sample_coverage: int,
    min_merge_distance: int,
    min_locus_length: int,
    suffix: str,
    tmpdir: Path,
) -> Path:
    """Merge per-sample BEDs into the shared loci BED used by assemble."""
    ref_info = tmpdir / "REF_info.txt"
    bed_dir = tmpdir / "beds"
    bed_paths = [bed_dir / f"{sname}{suffix}" for sname in snames]
    if not bed_paths:
        raise IPyradError(
            "No sample BED files were provided for shared locus delimiting."
        )
    out_bed = bed_dir / "loci.bed"
    # logger.warning(bed_paths)
    # for each bed get [chrom, start, end, nsamples, samplenames, A-present, B-present, ...]
    # Chr1    1789726 1789745 1       C       0       0       1       0
    # Chr1    1790639 1790658 1       B       0       1       0       0
    # Chr1    1792068 1792357 4       A,B,C,D 1       1       1       1
    # Chr1    1792357 1792384 3       A,B,C   1       1       1       0
    # Chr1    1792384 1792386 2       B,C     0       1       1       0
    # Chr1    1799627 1799701 1       B       0       1       0       0
    # Chr1    1810262 1810282 1       D       0       0       0       1
    cmd1 = [BIN_BED, "multiinter", "-g", str(ref_info)]
    cmd1 += ["-i"] + [str(p) for p in bed_paths]
    cmd1 += ["-names"] + snames

    # require at least MIN_SAMPLES_COVERAGE in each bed and print only first 5 cols
    # Chr1    1792068 1792357 4       A,B,C,D
    # Chr1    1792357 1792384 3       A,B,C
    # Chr1    2344873 2344902 3       A,B,D
    # Chr1    2665674 2665760 3       B,C,D
    # Chr1    2824851 2824932 4       A,B,C,D
    # Chr1    3045768 3045944 3       A,B,D
    cmd2 = [
        "awk",
        f'BEGIN{{OFS="\\t"}} $4>={int(min_sample_coverage)} {{print $1,$2,$3,$4,$5}}',
    ]

    # merge sub-intervals by MIN_MERGE_DISTANCE
    # Chr1    1792068 1792384 4
    # Chr1    2344873 2344902 3
    # Chr1    2665674 2665760 3
    # Chr1    2824851 2824932 4
    # Chr1    3045768 3045944 3
    cmd3 = ["sort", "-k1,1", "-k2,2n", "-T", str(tmpdir)]

    cmd4 = [
        BIN_BED,
        "merge",
        "-i",
        "-",
        "-d",
        str(int(min_merge_distance)),
        "-c",
        "4",
        "-o",
        "min",
    ]

    # cmd5: filter intervals shorter than MIN_LOCUS_LENGTH
    # Chr1    1792068 1792384 4
    # Chr1    2665674 2665760 3
    # Chr1    2824851 2824932 4
    # Chr1    3045768 3045944 3
    cmd5 = ["awk", "-v", f"L={min_locus_length}", 'BEGIN{OFS=FS="\t"} ($3-$2) >= L']
    cmd6 = [BIN_BED, "sort", "-i", "-", "-g", str(ref_info)]

    # run pipeline
    run_pipeline([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6], out_bed)
    return out_bed


def get_sample_depth_stats_in_final_loci(
    sname: str, loci_bed: Path, tmpdir: Path
) -> dict[str, float]:
    """Return per-sample depth summaries across the final shared loci BED.

    This reuses the existing per-sample bedgraph and computes mean locus depth
    over the full locus length, so uncovered portions contribute zero depth.
    """
    cov_bed = tmpdir / "beds" / f"{sname}.fragments.bedgraph"
    cmd = [
        BIN_BED,
        "intersect",
        "-a",
        str(loci_bed),
        "-b",
        str(cov_bed),
        "-wao",
    ]

    locus_depths: list[float] = []
    current_key: tuple[str, int, int] | None = None
    current_length = 0
    current_weighted_depth = 0.0

    def flush_current() -> None:
        nonlocal current_key, current_length, current_weighted_depth
        if current_key is None:
            return
        mean_depth = current_weighted_depth / current_length if current_length else 0.0
        locus_depths.append(float(mean_depth))

    # bedtools -wao emits the A intervals in order, potentially repeated for
    # multiple bedgraph overlaps, so we can aggregate each locus on the fly.
    for line in stream_pipeline_lines([cmd]):
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        locus_key = (chrom, start, end)

        if current_key != locus_key:
            flush_current()
            current_key = locus_key
            current_length = max(0, end - start)
            current_weighted_depth = 0.0

        overlap = int(fields[-1])
        if overlap <= 0:
            continue
        cov_field = fields[-2]
        if cov_field == ".":
            continue
        current_weighted_depth += overlap * float(cov_field)

    flush_current()

    if not locus_depths:
        return {
            "shared_loci_with_nonzero_depth": 0,
            "mean_depth_shared_loci": 0.0,
            "median_depth_shared_loci": 0.0,
            "mean_depth_nonzero_shared_loci": 0.0,
            "median_depth_nonzero_shared_loci": 0.0,
        }

    covs = np.array(locus_depths, dtype=float)
    nonzero = covs[covs > 0]
    return {
        "shared_loci_with_nonzero_depth": int(nonzero.size),
        "mean_depth_shared_loci": float(np.mean(covs)),
        "median_depth_shared_loci": float(np.median(covs)),
        "mean_depth_nonzero_shared_loci": float(np.mean(nonzero))
        if nonzero.size
        else 0.0,
        "median_depth_nonzero_shared_loci": float(np.median(nonzero))
        if nonzero.size
        else 0.0,
    }
