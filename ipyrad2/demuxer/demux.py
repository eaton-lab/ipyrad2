#!/usr/bin/env python

"""Some utilities used in demux.py for demultiplexing."""

from typing import Dict, Tuple, List, Iterator
import itertools
import shutil
from pathlib import Path
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass, field

from loguru import logger
import pandas as pd
from pandas.errors import ParserError
from ipyrad2.utils.kmers import (
    InferredJunctionSet,
    get_overhangs_from_barcoded_reads,
    validate_named_motif,
    validate_named_motif_list,
)
from ipyrad2.utils.names import get_name_to_fastq_dict
from ipyrad2.utils.seqs import AMBIGS, BADCHARS
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.demuxer.demux_pipeline import _demux_spool_dir, run_demux_pipeline
from ipyrad2.demuxer.demux_report import (
    DEMUX_STATS_PREFIX,
    format_logged_motif_set,
    format_preserved_file_preview,
    warn_multi_motif_inference,
    write_demux_stats,
)
from ipyrad2.demuxer.match import DemuxRunConfig, get_demux_mode_label, run_serial_demux
from ipyrad2.demuxer.sample_names import (
    is_technical_replicate_name,
    technical_replicate_base_name,
)


BASES = set("ACGTN")
BARCODE_CHARS = set("RKSYWMCATG")


def _barcode_patterns_by_length(
    barcode_pairs: Dict[str, Tuple[str, str]],
    read_end: int,
) -> Dict[int, Tuple[str, ...]]:
    """Return unique raw barcode patterns grouped by length for one read end."""
    grouped = {}
    for barcodes in barcode_pairs.values():
        barcode = barcodes[read_end]
        if not barcode:
            continue
        grouped.setdefault(len(barcode), set()).add(barcode)
    return {
        length: tuple(sorted(values))
        for length, values in sorted(grouped.items())
    }


def _manual_junction_set(
    motifs: Tuple[str, ...],
    *,
    offset: int = 0,
) -> InferredJunctionSet:
    """Build junction-set metadata for explicit user-entered overhangs."""
    return InferredJunctionSet(
        motifs=motifs,
        motif_counts=tuple(0 for _ in motifs),
        offset=offset,
        total_support=0,
        runner_up_offset_support=0,
        candidate_offsets=(offset,),
    )


def _format_motif_tuple(motifs: Tuple[str, ...]) -> str:
    """Return motifs in the same bracketed style used by demux logging."""
    return f"[{', '.join(motifs) if motifs else '<none>'}]"


def _expand_cuts(motifs: Tuple[str, ...]) -> List[bytes]:
    """Expand motifs by IUPAC resolution and one-mismatch barcode-style mutation."""
    expanded = set()
    for motif in motifs:
        if any(i in "RKSYWM" for i in motif):
            cuts = [
                "".join(AMBIGS[i][0] if i in "RKSYWM" else i for i in motif),
                "".join(AMBIGS[i][1] if i in "RKSYWM" else i for i in motif),
            ]
        else:
            cuts = [motif]
        expanded.update(cuts)
        expanded.update(itertools.chain(*[mutate(i) for i in cuts]))
    return [i.encode() for i in sorted(expanded, key=lambda item: (-len(item), item))]


def _expand_barcode_candidates(barcode: str, max_mismatch: int) -> Dict[str, int]:
    """Return acceptable barcode candidates mapped to their minimum mismatch distance."""
    best = {barcode: 0}
    if max_mismatch <= 0:
        return best

    frontier = {barcode}
    for distance in range(1, max_mismatch + 1):
        next_frontier = set()
        for seq in frontier:
            for candidate in mutate(seq):
                current = best.get(candidate)
                if current is None or distance < current:
                    best[candidate] = distance
                    next_frontier.add(candidate)
        frontier = next_frontier
        if not frontier:
            break
    return best


def _freeze_sample_map(mapping: Dict[bytes, set[str]]) -> Dict[bytes, Tuple[str, ...]]:
    """Convert a mutable sample-set map to deterministic tuples."""
    return {
        barcode: tuple(sorted(samples))
        for barcode, samples in sorted(mapping.items(), key=lambda item: item[0])
        if samples
    }


def _barcode_candidates_by_length(mapping: Dict[bytes, Tuple[str, ...]]) -> Dict[int, frozenset[bytes]]:
    """Group runtime barcode candidates by length for exact boundary matching."""
    grouped: Dict[int, set[bytes]] = defaultdict(set)
    for barcode in mapping:
        grouped[len(barcode)].add(barcode)
    return {
        length: frozenset(sorted(values))
        for length, values in sorted(grouped.items())
        if values
    }


def _barcode_sample_map_from_names(
    names_to_barcodes: Dict[str, Tuple[str, str]],
    read_end: int,
) -> Dict[bytes, Tuple[str, ...]]:
    """Build an exact barcode-to-samples map for one read end."""
    mapping: Dict[bytes, set[str]] = defaultdict(set)
    for sample_name, barcodes in names_to_barcodes.items():
        barcode = barcodes[read_end]
        if barcode:
            mapping[barcode.encode()].add(sample_name)
    return _freeze_sample_map(mapping)


def _collect_boundary_collisions(
    barcode_to_samples: Dict[bytes, Tuple[str, ...]],
    motifs: Tuple[str, ...],
    read_end: str,
    source: str,
) -> List[Dict[str, str]]:
    """Return real cross-sample collisions caused by motif occurrences inside barcodes."""
    collisions: List[Dict[str, str]] = []
    seen = set()
    for barcode, longer_samples in sorted(barcode_to_samples.items(), key=lambda item: item[0]):
        barcode_str = barcode.decode()
        longer_sample_set = set(longer_samples)
        for motif in motifs:
            max_start = len(barcode_str) - len(motif)
            for start in range(1, max_start + 1):
                if barcode_str[start:start + len(motif)] != motif:
                    continue
                prefix = barcode[:start]
                shorter_samples = barcode_to_samples.get(prefix)
                if not shorter_samples:
                    continue
                shorter_sample_set = set(shorter_samples)
                if not (shorter_sample_set - longer_sample_set or longer_sample_set - shorter_sample_set):
                    continue
                key = (
                    read_end,
                    source,
                    prefix.decode(),
                    motif,
                    barcode_str,
                    tuple(sorted(shorter_samples)),
                    tuple(sorted(longer_samples)),
                )
                if key in seen:
                    continue
                seen.add(key)
                collisions.append(
                    {
                        "read_end": read_end,
                        "source": source,
                        "prefix_barcode": prefix.decode(),
                        "motif": motif,
                        "full_barcode": barcode_str,
                        "prefix_samples": ",".join(sorted(shorter_samples)),
                        "full_samples": ",".join(sorted(longer_samples)),
                    }
                )
    return collisions


@dataclass
class Demux:
    fastqs: List[Path]
    """: List of Paths to fastq files, unpaired."""
    barcodes: Path
    """: Path to the barcodes file."""
    cutsite_1: str | None
    """: 5' restriction-site remnant / cutsite motif at the start of R1. Inferred if None."""
    cutsite_2: str | None
    """: 5' restriction-site remnant / cutsite motif at the start of R2. Inferred if None."""
    max_mismatch: int
    """: Max number of mismatches between barcodes. Checked for conflict."""
    cores: int
    """: max number of parallel cores."""
    chunksize: int
    """: max number of reads to process between writing to disk."""
    merge_technical_replicates: bool
    """: merge replicates or append -technical-replicate-X to names."""
    outdir: Path
    """: outdir/prefix is the dir where fastqs will be written."""
    i7: bool
    """: if True then demux on i7 index instead of inline barcode(s)."""
    disable_infer_cutsite_motifs: bool
    """: Skip cutsite motif inference."""
    max_reads: int | None
    """: subsample only the first N reads from each file (used for testing)."""
    max_reads_kmer: int
    """: Total reads sampled across files for junction inference."""
    log_level: str
    barcode_boundary_slack: int = 1
    """: Max 5-prime barcode-boundary offset allowed for inline barcode matching."""
    pigz: bool = False
    force: bool = False

    # attrs to be filled ----------------------------------------------
    _names_to_barcodes: Dict[str, Tuple[str, str]] = None
    """: A map of barcode strings to sample names, pre-expanded by off-by-N."""
    _filenames_to_fastqs: Dict[str, Tuple[Path, Path | None]] = field(default_factory=dict)
    """: Dict mapping parsed input names to SE or PE FASTQ tuples."""
    _pe: bool = True
    """: Whether the parsed input FASTQs are paired-end."""
    _cuts1: List[str] = None
    """: List of enzyme overhang sites to match on read1s."""
    _cuts2: List[str] = None
    """: List of enzyme overhang sites to match on read2s."""
    _barcodes_to_names: Dict[bytes, str] = None
    """: Dict of all acceptable runtime barcode keys mapped to sample names."""
    _barcodes_to_samples: Dict[bytes, Tuple[str, ...]] = None
    """: Runtime barcode combinations mapped to all matching sample names."""
    _barcode1_to_samples: Dict[bytes, Tuple[str, ...]] = None
    """: Runtime R1 barcode candidates mapped to all matching sample names."""
    _barcode2_to_samples: Dict[bytes, Tuple[str, ...]] = None
    """: Runtime R2 barcode candidates mapped to all matching sample names."""
    _barcode1_candidates_by_length: Dict[int, frozenset[bytes]] = None
    """: Runtime R1 barcode candidates grouped by length for exact boundary matching."""
    _barcode2_candidates_by_length: Dict[int, frozenset[bytes]] = None
    """: Runtime R2 barcode candidates grouped by length for exact boundary matching."""
    _barcode1_mismatch_by_barcode: Dict[bytes, int] = None
    """: Minimum mismatch distance for each acceptable R1 barcode candidate."""
    _barcode2_mismatch_by_barcode: Dict[bytes, int] = None
    """: Minimum mismatch distance for each acceptable R2 barcode candidate."""
    _file_stats: Dict[str, List] = None
    """: Store stats per raw data file (pair)."""
    _sample_stats: Dict[str, int] = None
    """: Dict to store n reads per sample."""
    _technical_replicates: Dict[str, List[str]] = field(default_factory=dict)
    _barcode_lengths1: Tuple[int, ...] = None
    _barcode_lengths2: Tuple[int, ...] = None
    _re1_motifs: Tuple[str, ...] = field(default_factory=tuple)
    _re2_motifs: Tuple[str, ...] = field(default_factory=tuple)
    _re1_inference: InferredJunctionSet | None = None
    _re2_inference: InferredJunctionSet | None = None
    _re1_detected_inference: InferredJunctionSet | None = None
    _re2_detected_inference: InferredJunctionSet | None = None
    _re1_source: str | None = None
    _re2_source: str | None = None
    _re1_motif_decision: str | None = None
    _re2_motif_decision: str | None = None
    _barcode_boundary_collisions: List[Dict[str, str]] = field(default_factory=list)

    def __post_init__(self):
        """Run subfunctions to setup object."""
        self._prepare_input_paths()
        self._prepare_barcode_table()
        self._prepare_output_artifacts()
        self._prepare_cutsite_motifs()
        self._prepare_matching_state()

    def _prepare_input_paths(self) -> None:
        """Resolve primary inputs before any demux-specific setup."""
        self._prepare_outdir_path()
        self._resolve_barcodes_path()
        self._load_input_fastqs()
        self._validate_runtime_args()

    def _prepare_barcode_table(self) -> None:
        """Load the barcode table and normalize output sample names."""
        self._load_barcode_table()
        self._sanitize_sample_names()

    def _prepare_output_artifacts(self) -> None:
        """Preflight demux-managed outputs in the destination outdir."""
        self._check_for_existing_outputs()

    def _prepare_cutsite_motifs(self) -> None:
        """Resolve user-entered or inferred cutsite motifs."""
        if self.i7:
            return
        self._validate_user_cutsite_motifs()
        if not self.disable_infer_cutsite_motifs:
            self._resolve_cutsite_motifs()
        self._ensure_cutsite_motifs_available()
        self._expand_cut_motifs()

    def _prepare_matching_state(self) -> None:
        """Build runtime barcode maps and collision checks for demux."""
        self._build_runtime_barcode_maps()
        if not self.i7:
            self._check_barcode_boundary_collisions()

    def _warn_partial_outputs(self) -> None:
        """Warn that demux may have left partial outputs behind after a failure."""
        logger.warning(
            "demux failed; output directory '{}' may contain partial files and should be removed before rerun.",
            self.outdir,
        )

    def run(self):
        """Run each file (pair) on separate demux engine(s)."""
        try:
            self._demultiplex()
        except KeyboardInterrupt:
            self._warn_partial_outputs()
            raise
        except Exception:
            self._warn_partial_outputs()
            raise
        self._write_stats()
        self._merge_cleanup()

    def _load_input_fastqs(self) -> None:
        self._filenames_to_fastqs = get_name_to_fastq_dict(self.fastqs, None, None)
        paired_states = {
            fastq_tuple[1] is not None
            for fastq_tuple in self._filenames_to_fastqs.values()
        }
        if len(paired_states) != 1:
            raise IPyradError(
                "some but not all files have R1 and R2 pairs. Check inputs."
            )
        self._pe = paired_states.pop()
        logger.info("Found {} data", "PE" if self._pe else "SE")

    def _validate_runtime_args(self) -> None:
        """Validate demux runtime arguments."""
        if self.cores < 1:
            raise IPyradError("cores must be >= 1.")
        if self.chunksize < 1:
            raise IPyradError("chunksize must be >= 1.")
        if self.max_reads is not None and self.max_reads < 1:
            raise IPyradError("max_reads must be >= 1 when set.")
        if self.max_reads_kmer < 1:
            raise IPyradError("max_reads_kmer must be >= 1.")
        if not 0 <= self.max_mismatch <= 2:
            raise IPyradError("max_mismatch must be between 0 and 2.")
        if self.barcode_boundary_slack not in (0, 1):
            raise IPyradError("barcode_boundary_slack must be 0 or 1.")

    def _prepare_outdir_path(self) -> None:
        """Normalize the outdir path and create missing parent directories."""
        self.outdir = Path(self.outdir).expanduser().resolve()
        if self.outdir.exists() and not self.outdir.is_dir():
            raise IPyradError(f"outdir '{self.outdir}' exists and is not a directory.")
        self.outdir.mkdir(parents=True, exist_ok=True)

    def _resolve_barcodes_path(self) -> None:
        """Resolve the barcodes path from one concrete file or one glob match."""
        bars = Path(self.barcodes)
        bpath = sorted(bars.parent.glob(bars.name), key=lambda item: str(item))
        if not bpath:
            raise IPyradError(f"No barcodes file found at {self.barcodes}")
        if len(bpath) > 1:
            preview = ", ".join(path.name for path in bpath[:5])
            suffix = "" if len(bpath) <= 5 else f", ... and {len(bpath) - 5} more"
            raise IPyradError(
                f"Barcode path pattern matches multiple files. Select exactly one: {preview}{suffix}"
            )
        self.barcodes = Path(bpath[0]).expanduser().resolve()
        if not self.barcodes.is_file():
            raise IPyradError(f"Barcodes path is not a file: {self.barcodes}")

    def _load_barcode_table(self) -> None:
        """Fill .names_to_barcodes dict w/ info from barcodes file.

        This logs a WARNING if technical replicates are detected to
        make sure the user is aware of how they are being handled.
        """
        # parse the tabular barcodes file on whitespace. Expects
        # there to be no header. There will be >=2 columns, >2 if
        # combinatorial barcodes.
        try:
            bardata = pd.read_csv(
                self.barcodes,
                header=None,
                sep=r"\s+",
                skip_blank_lines=True,
                dtype="string",
            )
        except ParserError as err:
            raise IPyradError(
                "Failed to parse barcodes file. Check that your sample\n"
                "names do not include spaces (invalid)"
            ) from err

        if bardata.empty or bardata.shape[1] < 2:
            raise IPyradError(
                "Barcodes file must contain at least two whitespace-delimited columns: "
                "sample and barcode1."
            )

        n_columns = bardata.shape[1]
        if self.i7:
            if n_columns > 2:
                logger.warning(
                    "Ignoring barcode2 and any extra barcode columns because "
                    "--i7 uses only barcode1."
                )
            bardata = bardata.iloc[:, :2].copy()
        else:
            # the dataframe COULD have >3 columns, in which case we will
            # discard any extra columns to keep at most 3.
            bardata = bardata.iloc[:, :3].copy()

        invalid_rows = bardata.isna().any(axis=1)
        if invalid_rows.any():
            bad_rows = ", ".join(str(i + 1) for i in bardata.index[invalid_rows][:5])
            suffix = "" if invalid_rows.sum() <= 5 else f", ... and {invalid_rows.sum() - 5} more"
            raise IPyradError(
                "Barcodes file contains incomplete rows. Each non-blank row must define "
                "sample and barcode columns. Bad row numbers: "
                f"{bad_rows}{suffix}"
            )

        # set names on barcodes dataframe
        if bardata.shape[1] == 2:
            bardata.columns = ["sample", "barcode1"]
            bardata["sample"] = bardata["sample"].astype("string").str.strip()
            bardata["barcode1"] = bardata["barcode1"].astype("string").str.upper().str.strip()
        else:
            bardata.columns = ["sample", "barcode1", "barcode2"]
            bardata["sample"] = bardata["sample"].astype("string").str.strip()
            bardata["barcode1"] = bardata["barcode1"].astype("string").str.upper().str.strip()
            bardata["barcode2"] = bardata["barcode2"].astype("string").str.upper().str.strip()

        required_columns = ["sample", "barcode1"] + (["barcode2"] if "barcode2" in bardata.columns else [])
        empty_rows = bardata[required_columns].eq("").any(axis=1)
        if empty_rows.any():
            bad_rows = ", ".join(str(i + 1) for i in bardata.index[empty_rows][:5])
            suffix = "" if empty_rows.sum() <= 5 else f", ... and {empty_rows.sum() - 5} more"
            raise IPyradError(
                "Barcodes file contains empty sample or barcode fields. Bad row numbers: "
                f"{bad_rows}{suffix}"
            )

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

        # make sure barcodes are valid characters and not monomorphic.
        for row in bardata.itertuples(index=False):
            validate_named_motif(
                row.barcode1,
                f"barcode1 for sample '{row.sample}'",
                allowed_chars=BARCODE_CHARS,
            )
            if hasattr(row, "barcode2"):
                validate_named_motif(
                    row.barcode2,
                    f"barcode2 for sample '{row.sample}'",
                    allowed_chars=BARCODE_CHARS,
                )

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
            self._barcode_lengths2 = tuple(sorted({len(i) for i in bardata["barcode2"]}))
        else:
            self._names_to_barcodes = dict(zip(
                bardata["sample"], ((i, "") for i in bardata["barcode1"])
            ))
            self._barcode_lengths2 = ()
        self._barcode_lengths1 = tuple(sorted({len(i) for i in bardata["barcode1"]}))
        # report to logger
        logger.debug(f"barcodes map:\n{bardata}")

    def _sanitize_sample_names(self) -> None:
        """Replace unsupported characters in sample names without allowing collisions."""
        sanitized: Dict[str, Tuple[str, str]] = {}
        collisions: Dict[str, List[str]] = defaultdict(list)

        for name, barcodes in self._names_to_barcodes.items():
            newname = name
            if any(i in name for i in BADCHARS):
                for badchar in BADCHARS:
                    newname = newname.replace(badchar, "_")
                logger.warning(f"changing name {name} to {newname} (bad characters).")
            collisions[newname].append(name)
            sanitized[newname] = barcodes

        bad = {newname: names for newname, names in collisions.items() if len(names) > 1}
        if bad:
            details = "; ".join(
                f"{newname} <- {', '.join(sorted(names))}"
                for newname, names in sorted(bad.items())
            )
            raise IPyradError(
                f"Sanitized sample names would collide. Revise the barcodes file names: {details}"
            )
        self._names_to_barcodes = sanitized

    def _final_output_sample_names(self) -> List[str]:
        """Return the final output sample names and reject unsafe name collisions."""
        output_to_sources: Dict[str, List[str]] = defaultdict(list)
        for name in self._names_to_barcodes:
            output_name = (
                technical_replicate_base_name(name)
                if self.merge_technical_replicates
                else name
            )
            output_to_sources[output_name].append(name)

        bad = []
        for output_name, sources in sorted(output_to_sources.items()):
            if len(sources) == 1:
                continue
            allowed_merge = self.merge_technical_replicates and all(
                is_technical_replicate_name(source)
                and technical_replicate_base_name(source) == output_name
                for source in sources
            )
            if not allowed_merge:
                bad.append(f"{output_name} <- {', '.join(sorted(sources))}")

        if bad:
            detail = "; ".join(bad)
            raise IPyradError(
                f"Final demux output sample names would collide. Revise sample names in the barcodes file: {detail}"
            )
        return sorted(output_to_sources)

    def _managed_output_artifacts(self) -> Tuple[Path, ...]:
        """Return the exact output artifacts managed by this demux invocation."""
        artifacts: List[Path] = []
        for sample_name in self._final_output_sample_names():
            artifacts.append(self.outdir / f"{sample_name}_R1.fastq.gz")
            artifacts.append(self.outdir / f"{sample_name}_R1.fastq")
            if self._pe:
                artifacts.append(self.outdir / f"{sample_name}_R2.fastq.gz")
                artifacts.append(self.outdir / f"{sample_name}_R2.fastq")
        return tuple(artifacts)

    def _warn_preserved_outdir_files(self, managed_artifacts: Tuple[Path, ...]) -> None:
        """Warn once for preserved stats files and unrelated FASTQ outputs."""
        managed_set = set(managed_artifacts)
        existing_stats = sorted(self.outdir.glob(f"{DEMUX_STATS_PREFIX}*.txt"))
        if existing_stats:
            logger.warning(
                "existing demux stats files are present in {} and will not be overwritten ({} total): {}",
                self.outdir,
                len(existing_stats),
                format_preserved_file_preview(existing_stats),
            )

        preserved_fastqs = [
            path for path in sorted(self.outdir.glob("*.fastq.gz"))
            if path not in managed_set
        ]
        if preserved_fastqs:
            logger.warning(
                "existing FASTQ.gz files are present in {} and will not be overwritten ({} total): {}",
                self.outdir,
                len(preserved_fastqs),
                format_preserved_file_preview(preserved_fastqs),
            )

    def _check_for_existing_outputs(self) -> None:
        """Preflight this run's managed artifacts and optionally remove them."""
        managed_artifacts = self._managed_output_artifacts()
        self._warn_preserved_outdir_files(managed_artifacts)
        spool_dir = _demux_spool_dir(self.outdir)

        existing_artifacts = [path for path in managed_artifacts if path.exists()]
        if (existing_artifacts or spool_dir.exists()) and not self.force:
            if spool_dir.exists():
                raise IPyradError(
                    f"Existing demux spool directory found in outdir: {spool_dir}. Use --force to overwrite."
                )
            raise IPyradError(
                "One or more files matching the expected output names exist in outdir. "
                "Use --force to overwrite."
            )
        if existing_artifacts and self.force:
            logger.info(
                "removing {} existing demux output artifact(s) from {} because --force was set",
                len(existing_artifacts),
                self.outdir,
            )
            for artifact in existing_artifacts:
                artifact.unlink()
        if spool_dir.exists() and self.force:
            logger.info(
                "removing existing demux spool directory from {} because --force was set",
                self.outdir,
            )
            if spool_dir.is_dir():
                shutil.rmtree(spool_dir)
            else:
                spool_dir.unlink()

    def _validate_user_cutsite_motifs(self) -> None:
        """Validate explicit user-provided cutsite motifs."""
        if self.cutsite_1:
            self._re1_motifs = validate_named_motif_list(self.cutsite_1, "R1 cutsite motif")
            self._re1_inference = _manual_junction_set(self._re1_motifs)
            self._re1_source = "manual"
            self._re1_motif_decision = "inference disabled; using manual motifs"
        if self.cutsite_2:
            self._re2_motifs = validate_named_motif_list(
                self.cutsite_2,
                "R2 cutsite motif",
                allow_empty=True,
            )
            self._re2_inference = _manual_junction_set(self._re2_motifs)
            self._re2_source = "manual"
            self._re2_motif_decision = "inference disabled; using manual motifs"

    def _ensure_cutsite_motifs_available(self) -> None:
        """Ensure required cutsite motifs are available before matching."""
        if not self._re1_motifs:
            raise IPyradError(
                "Cutsite motifs are required for demux. "
                "Provide --cutsite-1 or enable cutsite motif inference."
            )
        if self._barcode_lengths2 and self._pe and not self._re2_motifs:
            raise IPyradError(
                "Cutsite motifs are required on R2 for combinatorial inline barcodes. "
                "Provide --cutsite-2 or enable cutsite motif inference."
            )

    def _infer_cutsite_motifs(
        self,
        read_end: str,
        fastqs: List[Path],
        barcode_index: int,
    ) -> InferredJunctionSet:
        """Run barcode-aware kmer inference for one read end."""
        return get_overhangs_from_barcoded_reads(
            fastqs,
            _barcode_patterns_by_length(self._names_to_barcodes, barcode_index),
            20,
            self.max_reads_kmer,
            self.cores,
            self.log_level,
            label=f"{read_end} cutsite motif inference",
            max_barcode_boundary_slack=self.barcode_boundary_slack,
        )

    def _record_manual_motif_decision(
        self,
        read_end: str,
        manual_motifs: Tuple[str, ...],
        detected: InferredJunctionSet,
    ) -> str:
        """Log whether user-entered motifs agree with kmer-detected motifs."""
        logger.info(
            "{} cutsite motif inference detected {}",
            read_end,
            format_logged_motif_set(detected),
        )
        manual_text = _format_motif_tuple(manual_motifs)
        detected_text = _format_motif_tuple(detected.motifs)
        if set(manual_motifs) == set(detected.motifs):
            logger.info(
                "{} user-defined cutsite motifs {} match detected motifs {}; using user-defined motifs.",
                read_end,
                manual_text,
                detected_text,
            )
            return "manual motifs match detected motifs; using manual motifs"
        logger.warning(
            "{} user-defined cutsite motifs {} do not match detected motifs {}; "
            "letting the user-defined motif overrule the detected motif.",
            read_end,
            manual_text,
            detected_text,
        )
        return "manual motifs override detected motifs"

    def _resolve_read_end_cutsite_motifs(
        self,
        *,
        read_end: str,
        fastqs: List[Path],
        barcode_index: int,
        manual_motifs: Tuple[str, ...],
    ) -> Tuple[Tuple[str, ...], InferredJunctionSet, str, InferredJunctionSet | None, str]:
        """Infer, compare, and select cutsite motifs for one read end."""
        try:
            detected = self._infer_cutsite_motifs(read_end, fastqs, barcode_index)
        except IPyradError as exc:
            if not manual_motifs:
                raise
            logger.warning(
                "{} cutsite motif inference failed while user-defined motifs were provided: {}. "
                "Using user-defined motifs.",
                read_end,
                exc,
            )
            return (
                manual_motifs,
                _manual_junction_set(manual_motifs),
                "manual",
                None,
                "inference failed; using manual motifs",
            )

        if manual_motifs:
            decision = self._record_manual_motif_decision(read_end, manual_motifs, detected)
            return (
                manual_motifs,
                _manual_junction_set(manual_motifs),
                "manual",
                detected,
                decision,
            )

        warn_multi_motif_inference(read_end, detected, self.max_reads_kmer)
        return detected.motifs, detected, "auto", detected, "auto-selected"

    def _resolve_cutsite_motifs(self) -> None:
        """Use kmer analysis to detect cutsite motifs in sequences."""
        (
            self._re1_motifs,
            self._re1_inference,
            self._re1_source,
            self._re1_detected_inference,
            self._re1_motif_decision,
        ) = self._resolve_read_end_cutsite_motifs(
            read_end="R1",
            fastqs=[i[0] for i in self._filenames_to_fastqs.values()],
            barcode_index=0,
            manual_motifs=self._re1_motifs,
        )

        read2s = [i[1] for i in self._filenames_to_fastqs.values() if i[1] is not None]
        if read2s and self._barcode_lengths2:
            (
                self._re2_motifs,
                self._re2_inference,
                self._re2_source,
                self._re2_detected_inference,
                self._re2_motif_decision,
            ) = self._resolve_read_end_cutsite_motifs(
                read_end="R2",
                fastqs=read2s,
                barcode_index=1,
                manual_motifs=self._re2_motifs,
            )
        else:
            self._re2_inference = _manual_junction_set(self._re2_motifs)
            self._re2_source = "manual" if self._re2_motifs else None
            if self._re2_motifs:
                self._re2_motif_decision = "R2 inference not required; using manual motifs"

        logger.info(
            "cutsite motifs set to R1={}",
            format_logged_motif_set(self._re1_inference),
        )
        logger.info(
            "cutsite motifs set to R2={}",
            format_logged_motif_set(self._re2_inference),
        )

    def _expand_cut_motifs(self) -> None:
        """Fill `.cuts1` and `.cuts2` with ordered list of resolutions.

        Sequences will be searched for cut sites starting with the
        entered value and then proceeding to allow off-by-n matches.
        The first tested sequences will be the user entered value,
        with IUPAC resolved, followed by off-by-1 matches.
        """
        self._cuts1 = _expand_cuts(self._re1_motifs)
        self._cuts2 = _expand_cuts(self._re2_motifs) if self._re2_motifs else []

    def _build_runtime_barcode_maps(self) -> None:
        """Fills .barcodes_to_names with all acceptable barcodes: name.

        This updates the .barcodes_to_names from {str: Tuple[str,str]}
        to {str: str}.
        """
        combo_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        barcode1_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        barcode2_to_samples: Dict[bytes, set[str]] = defaultdict(set)
        barcode1_mismatch_by_barcode: Dict[bytes, int] = {}
        barcode2_mismatch_by_barcode: Dict[bytes, int] = {}
        for name, barcode in self._names_to_barcodes.items():
            bars1 = _expand_barcode_candidates(barcode[0], self.max_mismatch)
            bars2 = (
                _expand_barcode_candidates(barcode[1], self.max_mismatch)
                if barcode[1]
                else {}
            )

            for bar1, distance1 in bars1.items():
                bar1_bytes = bar1.encode()
                barcode1_to_samples[bar1_bytes].add(name)
                current = barcode1_mismatch_by_barcode.get(bar1_bytes)
                if current is None or distance1 < current:
                    barcode1_mismatch_by_barcode[bar1_bytes] = distance1

            for bar2, distance2 in bars2.items():
                bar2_bytes = bar2.encode()
                barcode2_to_samples[bar2_bytes].add(name)
                current = barcode2_mismatch_by_barcode.get(bar2_bytes)
                if current is None or distance2 < current:
                    barcode2_mismatch_by_barcode[bar2_bytes] = distance2

            if not bars2:
                for bar1 in bars1:
                    barc = bar1.encode()
                    combo_to_samples[barc].add(name)
                continue

            for bar1, bar2 in itertools.product(bars1, bars2):
                barc = f"{bar1}_{bar2}".encode()
                combo_to_samples[barc].add(name)
        self._barcodes_to_samples = _freeze_sample_map(combo_to_samples)
        ambiguous = {
            barcode: samples
            for barcode, samples in self._barcodes_to_samples.items()
            if len(samples) > 1
        }
        if ambiguous:
            details = "; ".join(
                f"{barcode.decode()} -> {', '.join(samples)}"
                for barcode, samples in list(sorted(ambiguous.items()))[:6]
            )
            if len(ambiguous) > 6:
                details += f"; ... and {len(ambiguous) - 6} more"
            raise IPyradError(
                "Barcode mismatch expansion creates ambiguous barcode candidates that "
                "map to multiple samples. Lower --max-mismatch or revise the barcodes. "
                f"{details}"
            )
        self._barcodes_to_names = {
            barcode: samples[0]
            for barcode, samples in self._barcodes_to_samples.items()
        }
        self._barcode1_to_samples = _freeze_sample_map(barcode1_to_samples)
        self._barcode2_to_samples = _freeze_sample_map(barcode2_to_samples)
        self._barcode1_candidates_by_length = _barcode_candidates_by_length(self._barcode1_to_samples)
        self._barcode2_candidates_by_length = _barcode_candidates_by_length(self._barcode2_to_samples)
        self._barcode1_mismatch_by_barcode = {
            barcode: barcode1_mismatch_by_barcode[barcode]
            for barcode in sorted(barcode1_mismatch_by_barcode)
        }
        self._barcode2_mismatch_by_barcode = {
            barcode: barcode2_mismatch_by_barcode[barcode]
            for barcode in sorted(barcode2_mismatch_by_barcode)
        }

    def _check_barcode_boundary_collisions(self) -> None:
        """Detect selected motifs that create real cross-sample barcode-boundary collisions."""
        collisions: List[Dict[str, str]] = []
        exact_r1 = _barcode_sample_map_from_names(self._names_to_barcodes, 0)
        collisions.extend(_collect_boundary_collisions(exact_r1, self._re1_motifs, "R1", "exact"))
        collisions.extend(_collect_boundary_collisions(
            self._barcode1_to_samples,
            self._re1_motifs,
            "R1",
            "runtime",
        ))

        if self._barcode_lengths2:
            exact_r2 = _barcode_sample_map_from_names(self._names_to_barcodes, 1)
            collisions.extend(_collect_boundary_collisions(exact_r2, self._re2_motifs, "R2", "exact"))
            collisions.extend(_collect_boundary_collisions(
                self._barcode2_to_samples,
                self._re2_motifs,
                "R2",
                "runtime",
            ))

        self._barcode_boundary_collisions = collisions
        if not collisions:
            return

        detail = "; ".join(
            (
                f"{row['source']} {row['read_end']} {row['prefix_barcode']} + {row['motif']} "
                f"-> {row['full_barcode']} "
                f"({row['prefix_samples']} vs {row['full_samples']})"
            )
            for row in collisions[:6]
        )
        if len(collisions) > 6:
            detail += f"; ... and {len(collisions) - 6} more"

        if self._barcode_lengths2:
            logger.warning(
                "restriction motif(s) create barcode-boundary collisions that may require joint R1/R2 resolution: {}",
                detail,
            )
            return

        raise IPyradError(
            "restriction motif(s) create unrecoverable R1 barcode-boundary collisions for single-inline demux: "
            f"{detail}"
        )

    def _demultiplex(self) -> None:
        """Demultiplex each raw FASTQ tuple through the serial or pipeline runner."""
        config = self._build_run_config()

        # single reader that intermittently writes pipeline
        if self.cores == 1 and not self.pigz:
            logger.info(f"demultiplexing on {get_demux_mode_label(config)}")
            jobs = {}
            for fname, fastq_tuple in self._filenames_to_fastqs.items():
                short = tuple(i.name if i else "" for i in fastq_tuple)
                logger.info(f"processing {fname} {short}")
                jobs[fname] = run_serial_demux(fastq_tuple, config, workers=self.cores)
            self._file_stats = jobs

        # multiple readers and writers pipeline
        else:
            self._file_stats = run_demux_pipeline(
                self._filenames_to_fastqs,
                config,
                self.cores,
            )

        # record of stats per sample from barmatch returned objects
        self._collect_sample_stats()

    def _collect_sample_stats(self) -> None:
        """Aggregate per-file demux stats into sample-level counters."""
        self._sample_stats = Counter()
        for _fname, stats in self._file_stats.items():
            for sname, hits in stats[2].items():
                self._sample_stats[sname] += hits
                if is_technical_replicate_name(sname):
                    short_name = technical_replicate_base_name(sname)
                    self._sample_stats[short_name] += hits
                    replicate_names = self._technical_replicates.setdefault(short_name, [])
                    if sname not in replicate_names:
                        self._technical_replicates[short_name].append(sname)

        for sname in self._names_to_barcodes:
            if not self._sample_stats[sname]:
                logger.warning(f"Sample {sname} has 0 reads.")
                self._sample_stats[sname] = 0

        for short_name, replicate_names in self._technical_replicates.items():
            replicate_names.sort()

    def _build_run_config(self) -> DemuxRunConfig:
        """Return a serializable configuration for demux worker processes."""
        return DemuxRunConfig(
            barcodes_to_names=self._barcodes_to_names,
            barcodes_to_samples=self._barcodes_to_samples,
            barcode1_to_samples=self._barcode1_to_samples,
            barcode2_to_samples=self._barcode2_to_samples,
            barcode1_candidates_by_length=self._barcode1_candidates_by_length,
            barcode2_candidates_by_length=self._barcode2_candidates_by_length,
            barcode1_mismatch_by_barcode=self._barcode1_mismatch_by_barcode,
            barcode2_mismatch_by_barcode=self._barcode2_mismatch_by_barcode,
            barcode_lengths1=self._barcode_lengths1,
            barcode_lengths2=self._barcode_lengths2,
            barcode_boundary_slack=self.barcode_boundary_slack,
            cuts1=self._cuts1,
            cuts2=self._cuts2,
            merge_technical_replicates=self.merge_technical_replicates,
            outdir=self.outdir,
            chunksize=self.chunksize,
            max_reads=self.max_reads,
            i7=self.i7,
            log_level=self.log_level,
            pigz=self.pigz,
        )

    def _write_stats(self) -> None:
        """Write the numbered demux stats report."""
        write_demux_stats(
            outdir=self.outdir,
            file_stats=self._file_stats,
            sample_stats=self._sample_stats,
            names_to_barcodes=self._names_to_barcodes,
            barcodes_to_names=self._barcodes_to_names,
            i7=self.i7,
            re1_source=self._re1_source,
            re1_inference=self._re1_inference,
            re1_detected_inference=self._re1_detected_inference,
            re1_motif_decision=self._re1_motif_decision,
            re2_source=self._re2_source,
            re2_inference=self._re2_inference,
            re2_detected_inference=self._re2_detected_inference,
            re2_motif_decision=self._re2_motif_decision,
            barcode_boundary_collisions=self._barcode_boundary_collisions,
        )

    def _merge_cleanup(self) -> None:
        """Remove keys from _sample_stats for merging."""
        if not self.merge_technical_replicates:
            for key in self._technical_replicates.keys():
                self._sample_stats.pop(key, None)
        else:
            for key, value in self._technical_replicates.items():
                for rep in value:
                    self._sample_stats.pop(rep, None)


######################################################################
######################################################################

def mutate(barcode: str) -> Iterator[str]:
    """Mutate a sequence by 1 base (ACGT). Used for barcode mismatch."""
    for pos, _ in enumerate(barcode):
        for sub in BASES:
            newbar = list(barcode)
            newbar[pos] = sub
            yield "".join(newbar)


def run_demuxer(**kwargs):
    """Command-line wrapper for Demux."""
    tool = Demux(**kwargs)
    tool.run()
    # DATA.params.max_barcode_mismatch = 1
    # DATA.hackers.demultiplex_on_i7_tags = True
