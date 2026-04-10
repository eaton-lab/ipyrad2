import csv
import re
import random
from pathlib import Path

import pandas as pd
import pytest

from ipyrad2.denovo import denovo as denovo_module
from ipyrad2.denovo import align as align_module
from ipyrad2.denovo import graph as graph_module
from ipyrad2.denovo.graph_split import split_component as split_component_constrained
from ipyrad2.denovo import cluster as cluster_module
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.parallel import PipelineTimeoutError


def _alignment_summary(
    *,
    total_loci: int = 1,
    single_sequence_loci: int = 0,
    identical_sequence_loci: int = 0,
    mafft_required_loci: int = 0,
    mafft_threads_per_job: int = 0,
    mafft_worker_processes: int = 0,
    alignment_mode: str = "mafft",
    mafft_timeout_seconds: int = 0,
    joined_spacer_loci: int = 0,
    mixed_reconciled_spacer_loci: int = 0,
    stripped_output_loci: int = 1,
    output_spacer_length: int = 0,
) -> align_module.AlignmentRunSummary:
    return align_module.AlignmentRunSummary(
        total_loci=total_loci,
        single_sequence_loci=single_sequence_loci,
        identical_sequence_loci=identical_sequence_loci,
        mafft_required_loci=mafft_required_loci,
        mafft_threads_per_job=mafft_threads_per_job,
        mafft_worker_processes=mafft_worker_processes,
        alignment_mode=alignment_mode,
        mafft_timeout_seconds=mafft_timeout_seconds,
        joined_spacer_loci=joined_spacer_loci,
        mixed_reconciled_spacer_loci=mixed_reconciled_spacer_loci,
        stripped_output_loci=stripped_output_loci,
        output_spacer_length=output_spacer_length,
    )


def _make_executable(tmp_path: Path, name: str) -> Path:
    path = tmp_path / name
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _patch_required_binaries(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> tuple[Path, Path]:
    vsearch = _make_executable(tmp_path, "vsearch")
    mafft = _make_executable(tmp_path, "mafft")
    monkeypatch.setattr(denovo_module, "BIN_VSEARCH", str(vsearch))
    monkeypatch.setattr(denovo_module, "BIN_MAFFT", str(mafft))
    return vsearch, mafft


def _summary_header() -> str:
    return (
        "sample\tcluster_id\tseed\tlength\tcluster_length\tn_unique\tn_reads\t"
        "record_type\tcluster_sequence\tarm_boundary\n"
    )


def _summary_row(
    *,
    sample: str,
    cluster_id: int,
    seed: str,
    cluster_sequence: str,
    record_type: str,
    n_unique: int = 1,
    n_reads: int = 1,
    length: int | None = None,
    cluster_length: int | None = None,
    arm_boundary: int | None = None,
) -> str:
    cluster_sequence = cluster_sequence.upper()
    cluster_length = (
        len(cluster_sequence) if cluster_length is None else int(cluster_length)
    )
    arm_boundary = cluster_length if arm_boundary is None else int(arm_boundary)
    if length is None:
        uses_spacer = record_type == "joined" and arm_boundary < cluster_length
        length = cluster_length + (
            cluster_module.CLUSTER_JOINED_SPACER_LEN if uses_spacer else 0
        )
    return (
        f"{sample}\t{cluster_id}\t{seed}\t{length}\t{cluster_length}\t"
        f"{n_unique}\t{n_reads}\t{record_type}\t{cluster_sequence}\t{arm_boundary}\n"
    )


def _unique_test_sequence(idx: int, width: int = 8) -> str:
    """Return a deterministic per-index sequence-like token for test fixtures."""
    return f"A{idx:0{width}d}"


def _read_graph_output_tables(workdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load graph-stage mapping/stats outputs written beside one workdir."""
    return (
        pd.read_csv(workdir.parent / "denovo.loci.mapping.tsv", sep="\t"),
        pd.read_csv(workdir.parent / "denovo.loci.stats.tsv", sep="\t"),
    )


def _write_audit_summary(audit_dir: Path, rows: list[dict[str, object]]) -> None:
    """Write one components.summary.tsv fixture with the canonical columns."""
    audit_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=graph_module.AUDIT_SUMMARY_FIELDS).to_csv(
        audit_dir / "components.summary.tsv",
        sep="\t",
        index=False,
    )


def _report_has_value_line(text: str, key: str, value: str) -> bool:
    """Return True when one aligned report line matches `key  value`."""
    return (
        re.search(rf"^{re.escape(key)}\s+{re.escape(value)}$", text, re.MULTILINE)
        is not None
    )


def test_iter_status_records_splits_crlf_and_preserves_trailing_partial() -> None:
    records, pending = denovo_module._iter_status_records(
        "Reading file 100%\rSearching 1%\rSearching 2%\nCreating index 100%\rpartial"
    )

    assert records == [
        "Reading file 100%",
        "Searching 1%",
        "Searching 2%",
        "Creating index 100%",
    ]
    assert pending == "partial"


def test_extract_searching_percent_ignores_non_searching_status() -> None:
    assert denovo_module._extract_searching_percent("Searching 17%") == 17
    assert denovo_module._extract_searching_percent("Searching 100%") == 100
    assert denovo_module._extract_searching_percent("Reading file 100%") is None
    assert (
        denovo_module._extract_searching_percent(
            "Matching unique query sequences: 0 of 10 (0.00%)"
        )
        is None
    )


def test_write_stripped_clustering_fasta_strips_joined_spacer(tmp_path: Path) -> None:
    joined = tmp_path / "sample.joined.fa"
    merged = tmp_path / "sample.merged.fa"
    out_fa = tmp_path / "sample.cluster.fa"
    joined.write_text(f">s1;J1\nAAA{'N' * 24}TTT\n", encoding="utf-8")
    merged.write_text(">s1;M1\nAAATTT\n", encoding="utf-8")

    seed_to_meta = denovo_module._write_stripped_clustering_fasta(
        out_fa, [joined, merged]
    )

    assert out_fa.read_text(encoding="utf-8") == ">s1;J1\nAAATTT\n>s1;M1\nAAATTT\n"
    assert seed_to_meta == {
        "s1;J1": ("joined", 3),
        "s1;M1": ("merged", 6),
    }


def test_build_sample_summary_rehydrates_joined_consensus_from_stripped_cluster(
    tmp_path: Path,
) -> None:
    workdir = tmp_path
    (workdir / "s1.joined.fa").write_text(
        f">s1;J1\nAAA{'N' * 24}TTT\n", encoding="utf-8"
    )
    (workdir / "s1.consensus.fa").write_text(
        ">centroid=s1;J1;size=9;seqs=2\nAAATTT\n", encoding="utf-8"
    )
    (workdir / "s1.clusters.tsv").write_text(
        "S\t0\t6\t*\t*\t*\t*\t*\ts1;J1;size=9\t*\n"
        "H\t0\t6\t100.0\t+\t0\t0\t0\ts1;J2;size=4\ts1;J1;size=9\n",
        encoding="utf-8",
    )

    df = cluster_module.build_sample_summary("s1", workdir)

    assert df.loc[0, "cluster_sequence"] == "AAATTT"
    assert df.loc[0, "record_type"] == "joined"
    assert df.loc[0, "arm_boundary"] == 3
    assert df.loc[0, "length"] == 30


def test_build_sample_summary_uses_injected_metadata_without_raw_fastas(
    tmp_path: Path,
) -> None:
    workdir = tmp_path
    (workdir / "s1.consensus.fa").write_text(
        ">centroid=s1;J1;size=9;seqs=2\nAAATTT\n",
        encoding="utf-8",
    )
    (workdir / "s1.clusters.tsv").write_text(
        "S\t0\t6\t*\t*\t*\t*\t*\ts1;J1;size=9\t*\n"
        "H\t0\t6\t100.0\t+\t0\t0\t0\ts1;J2;size=4\ts1;J1;size=9\n",
        encoding="utf-8",
    )

    df = cluster_module.build_sample_summary(
        "s1",
        workdir,
        seed_to_meta={"s1;J1": ("joined", 3)},
    )

    assert df.loc[0, "record_type"] == "joined"
    assert df.loc[0, "arm_boundary"] == 3
    assert df.loc[0, "length"] == 30


def test_run_denovo_rejects_mixed_input_layout(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        denovo_module,
        "get_name_to_fastq_dict",
        lambda *args, **kwargs: {
            "sample_a": (tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz"),
            "sample_b": (tmp_path / "b.fastq.gz", None),
        },
    )

    with pytest.raises(IPyradError, match="consistently single-end or paired-end"):
        denovo_module.run_denovo(
            fastqs=[tmp_path / "reads.fastq.gz"],
            outdir=tmp_path / "out",
            within_similarity=0.95,
            across_similarity=0.85,
            min_derep_size=2,
            min_length=35,
            min_merge_overlap=20,
            max_merge_diffs=4,
            delim_str=None,
            delim_idx=1,
            allow_reverse_complement=False,
            cores=6,
            threads=3,
            no_alignment=False,
            force=False,
            imap=None,
            use_all_samples=False,
            keep_intermediates=False,
            log_level="INFO",
        )


def test_run_denovo_requires_working_binaries(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        denovo_module,
        "get_name_to_fastq_dict",
        lambda *args, **kwargs: {"sample_a": (tmp_path / "a.fastq.gz", None)},
    )
    monkeypatch.setattr(denovo_module, "BIN_VSEARCH", str(tmp_path / "missing-vsearch"))
    monkeypatch.setattr(denovo_module, "BIN_MAFFT", str(tmp_path / "missing-mafft"))

    with pytest.raises(IPyradError, match="vsearch binary is not executable"):
        denovo_module.run_denovo(
            fastqs=[tmp_path / "reads.fastq.gz"],
            outdir=tmp_path / "out",
            within_similarity=0.95,
            across_similarity=0.85,
            min_derep_size=2,
            min_length=35,
            min_merge_overlap=20,
            max_merge_diffs=4,
            delim_str=None,
            delim_idx=1,
            allow_reverse_complement=False,
            cores=6,
            threads=3,
            no_alignment=False,
            force=False,
            imap=None,
            use_all_samples=False,
            keep_intermediates=False,
            log_level="INFO",
        )


def test_run_denovo_writes_curated_outputs_and_cleans_workdir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sample_fastq = tmp_path / "sample.fastq.gz"
    sample_fastq.write_text("", encoding="utf-8")
    outdir = tmp_path / "OUT"
    calls: dict[str, object] = {}
    vsearch_binary, mafft_binary = _patch_required_binaries(monkeypatch, tmp_path)

    monkeypatch.setattr(
        denovo_module,
        "get_name_to_fastq_dict",
        lambda *args, **kwargs: {"sample_a": (sample_fastq, None)},
    )

    def fake_run_with_pool(jobs, log_level, max_workers, msg):
        calls["pool"] = {
            "keys": list(jobs),
            "log_level": log_level,
            "max_workers": max_workers,
            "msg": msg,
        }
        return {"sample_a": None}

    def fake_concat_summaries(outdir):
        (outdir / "concat.summary.tsv").write_text(
            _summary_header()
            + _summary_row(
                sample="sample_a",
                cluster_id=0,
                seed="sample_a;S1",
                cluster_sequence="ACGTACGTAA",
                record_type="single",
                n_reads=5,
                length=10,
            ),
            encoding="utf-8",
        )
        return pd.DataFrame()

    def fake_vsearch_cluster_across(outdir, summary_tsv, across_similarity, threads):
        calls["across"] = {
            "outdir": outdir,
            "summary_tsv": summary_tsv,
            "across_similarity": across_similarity,
            "threads": threads,
        }

    def fake_make_global_tables(
        outdir,
        cores,
        log_level,
        within_similarity,
    ):
        mapping = pd.DataFrame(
            [
                {
                    "locus": 1,
                    "component_id": 1,
                    "subcomponent_id": 1,
                    "locus_name": "locus_1_1",
                    "contract_group": "contract_1_1",
                    "sample": "sample_a",
                    "n_reads": 5,
                    "n_unique": 1,
                    "length": 10,
                    "cluster_length": 10,
                    "merged": 0,
                    "record_type": "single",
                    "cluster_id": 0,
                    "core": "sample_a;S1",
                }
            ]
        )
        stats = pd.DataFrame(
            [
                {
                    "locus": 1,
                    "component_id": 1,
                    "subcomponent_id": 1,
                    "locus_name": "locus_1_1",
                    "n_samples": 1,
                    "n_cores": 1,
                    "n_contracted_groups": 0,
                    "n_reads_sum": 5,
                    "n_reads_mean": 5.0,
                    "n_reads_std": 0.0,
                    "length_mean": 10.0,
                    "length_std": 0.0,
                    "merged_freq": 0.0,
                    "samples": "sample_a",
                }
            ]
        )
        mapping.to_csv(outdir.parent / "denovo.loci.mapping.tsv", sep="\t", index=False)
        stats.to_csv(outdir.parent / "denovo.loci.stats.tsv", sep="\t", index=False)
        calls["graph_cores"] = cores
        calls["graph_log_level"] = log_level
        calls["graph_within_similarity"] = within_similarity
        return graph_module.GraphTableSummary(
            loci_written=int(stats.shape[0]),
            consensus_records=int(mapping.shape[0]),
        )

    def fake_write_ordered_consensus_stream_to_file(**kwargs):
        calls["consensus"] = kwargs
        kwargs["out_fa"].write_text(">locus_1\nACGTACGTAA\n", encoding="utf-8")
        return _alignment_summary(
            total_loci=1,
            single_sequence_loci=1,
            alignment_mode="mafft",
        )

    monkeypatch.setattr(denovo_module, "run_with_pool", fake_run_with_pool)
    monkeypatch.setattr(
        denovo_module,
        "build_sample_summary",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError(
                "run_denovo should not build sample summaries in the main process"
            )
        ),
    )
    monkeypatch.setattr(denovo_module, "concat_summaries", fake_concat_summaries)
    monkeypatch.setattr(
        denovo_module, "vsearch_cluster_across", fake_vsearch_cluster_across
    )
    monkeypatch.setattr(denovo_module, "make_global_tables", fake_make_global_tables)
    monkeypatch.setattr(
        denovo_module,
        "write_ordered_consensus_stream_to_file",
        fake_write_ordered_consensus_stream_to_file,
    )

    denovo_module.run_denovo(
        fastqs=[sample_fastq],
        outdir=outdir,
        within_similarity=0.95,
        across_similarity=0.85,
        min_derep_size=2,
        min_length=35,
        min_merge_overlap=20,
        max_merge_diffs=4,
        delim_str=None,
        delim_idx=1,
        allow_reverse_complement=False,
        cores=6,
        threads=3,
        no_alignment=False,
        force=False,
        imap=None,
        use_all_samples=False,
        keep_intermediates=False,
        log_level="INFO",
    )

    assert (outdir / "denovo_reference.fa").exists()
    assert (outdir / "denovo.loci.mapping.tsv").exists()
    assert (outdir / "denovo.loci.stats.tsv").exists()
    assert (outdir / "denovo.stats.txt").exists()
    assert not (outdir / denovo_module.WORKDIR_NAME).exists()

    stats_text = (outdir / "denovo.stats.txt").read_text(encoding="utf-8")
    assert "# Inputs" in stats_text
    assert "# Clustering Parameters" in stats_text
    assert "# Denovo Summary" in stats_text
    assert "# Locus QC" in stats_text
    assert "# Component QC" in stats_text
    assert "# Component Node Summary" in stats_text
    assert "# Selected Sample Summary" in stats_text
    assert "# Locus Occupancy" in stats_text
    assert "# Runtime" in stats_text
    assert "# Outputs" in stats_text
    assert _report_has_value_line(stats_text, "selected_sample_count", "1")
    assert _report_has_value_line(stats_text, "sample_selection_mode", "all")
    assert _report_has_value_line(stats_text, "vsearch_binary", str(vsearch_binary))
    assert _report_has_value_line(stats_text, "mafft_binary", str(mafft_binary))
    assert _report_has_value_line(stats_text, "keep_intermediates", "False")
    assert _report_has_value_line(stats_text, "alignment_mode", "mafft")
    assert _report_has_value_line(stats_text, "vsearch_threads_per_job", "3")
    assert _report_has_value_line(stats_text, "across_vsearch_threads", "6")
    assert _report_has_value_line(
        stats_text, "duplicated_component_reconciliation", "same-sample graph"
    )
    assert _report_has_value_line(stats_text, "mafft_threads_per_job", "0")
    assert _report_has_value_line(stats_text, "single_sequence_loci", "1")
    assert _report_has_value_line(stats_text, "mafft_required_loci", "0")
    assert _report_has_value_line(stats_text, "stripped_output_loci", "1")
    assert "sample_names:" not in stats_text
    assert "occupancy_distribution:" not in stats_text
    assert _report_has_value_line(stats_text, "selected_sample_count", "1")
    assert re.search(r"^sample_a\s+1\s+5\s+0\s+0\s+1\s*$", stats_text, re.MULTILINE)
    assert re.search(r"^1\s+1\s+1\.000000\s*$", stats_text, re.MULTILINE)

    assert calls["pool"] == {
        "keys": ["sample_a"],
        "log_level": "INFO",
        "max_workers": 2,
        "msg": "Dereplicating and clustering",
    }
    assert calls["consensus"]["mapping_tsv"] == outdir / "denovo.loci.mapping.tsv"
    assert (
        calls["consensus"]["summary_tsv"]
        == outdir / denovo_module.WORKDIR_NAME / "concat.summary.tsv"
    )
    assert calls["consensus"]["out_fa"] == outdir / "denovo_reference.fa"
    assert calls["consensus"]["mafft_binary"] == str(mafft_binary)
    assert calls["consensus"]["cores"] == 6
    assert calls["consensus"]["alignment_mode"] == "mafft"
    assert calls["graph_cores"] == 6
    assert calls["graph_log_level"] == "INFO"
    assert (
        calls["across"]["summary_tsv"]
        == outdir / denovo_module.WORKDIR_NAME / "concat.summary.tsv"
    )
    assert calls["across"]["threads"] == 6


def test_run_denovo_keep_intermediates_preserves_workdir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sample_fastq = tmp_path / "sample.fastq.gz"
    sample_fastq.write_text("", encoding="utf-8")
    outdir = tmp_path / "OUT"
    _patch_required_binaries(monkeypatch, tmp_path)

    monkeypatch.setattr(
        denovo_module,
        "get_name_to_fastq_dict",
        lambda *args, **kwargs: {"sample_a": (sample_fastq, None)},
    )
    monkeypatch.setattr(
        denovo_module, "run_with_pool", lambda *args, **kwargs: {"sample_a": None}
    )
    monkeypatch.setattr(
        denovo_module,
        "build_sample_summary",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError(
                "run_denovo should not build sample summaries in the main process"
            )
        ),
    )
    monkeypatch.setattr(
        denovo_module,
        "concat_summaries",
        lambda outdir: (outdir / "concat.summary.tsv").write_text(
            _summary_header()
            + _summary_row(
                sample="sample_a",
                cluster_id=0,
                seed="sample_a;S1",
                cluster_sequence="ACGTACGTAA",
                record_type="single",
                n_reads=5,
                length=10,
            ),
            encoding="utf-8",
        ),
    )
    monkeypatch.setattr(
        denovo_module, "vsearch_cluster_across", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        denovo_module,
        "make_global_tables",
        lambda outdir, cores, log_level, within_similarity: (
            pd.DataFrame(
                [
                    {
                        "locus": 1,
                        "locus_name": "locus_1_1",
                        "contract_group": "contract_1_1",
                        "core": "sample_a;S1",
                    }
                ]
            ).to_csv(
                outdir.parent / "denovo.loci.mapping.tsv",
                sep="\t",
                index=False,
            ),
            pd.DataFrame(
                [
                    {
                        "locus": 1,
                        "n_samples": 1,
                        "n_cores": 1,
                        "duplicated_component": False,
                        "used_reconciliation": False,
                    }
                ]
            ).to_csv(
                outdir.parent / "denovo.loci.stats.tsv",
                sep="\t",
                index=False,
            ),
            graph_module.GraphTableSummary(
                loci_written=1,
                consensus_records=1,
            ),
        )[-1],
    )
    monkeypatch.setattr(
        denovo_module,
        "write_ordered_consensus_stream_to_file",
        lambda **kwargs: (
            kwargs["out_fa"].write_text(">locus_1_1\nACGT\n", encoding="utf-8"),
            _alignment_summary(
                total_loci=1, single_sequence_loci=1, alignment_mode="mafft"
            ),
        )[1],
    )

    denovo_module.run_denovo(
        fastqs=[sample_fastq],
        outdir=outdir,
        within_similarity=0.95,
        across_similarity=0.85,
        min_derep_size=2,
        min_length=35,
        min_merge_overlap=20,
        max_merge_diffs=4,
        delim_str=None,
        delim_idx=1,
        allow_reverse_complement=False,
        cores=6,
        threads=3,
        no_alignment=False,
        force=False,
        imap=None,
        use_all_samples=False,
        keep_intermediates=True,
        log_level="INFO",
    )

    assert (outdir / denovo_module.WORKDIR_NAME).exists()


def test_prepare_output_paths_uses_current_renamed_denovo_outputs(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "OUT"
    workdir, outputs = denovo_module._prepare_output_paths(outdir, force=False)

    assert workdir == outdir / denovo_module.WORKDIR_NAME
    assert workdir.exists()
    assert outputs["mapping"] == outdir / "denovo.loci.mapping.tsv"
    assert outputs["loci_stats"] == outdir / "denovo.loci.stats.tsv"
    assert "qc" not in outputs


def test_prepare_output_paths_force_removes_current_renamed_outputs(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "OUT"
    outdir.mkdir()
    (outdir / "denovo.loci.mapping.tsv").write_text("old mapping\n", encoding="utf-8")
    (outdir / "denovo.loci.stats.tsv").write_text("old stats\n", encoding="utf-8")

    workdir, outputs = denovo_module._prepare_output_paths(outdir, force=True)

    assert workdir.exists()
    assert outputs["mapping"] == outdir / "denovo.loci.mapping.tsv"
    assert outputs["loci_stats"] == outdir / "denovo.loci.stats.tsv"
    assert not outputs["mapping"].exists()
    assert not outputs["loci_stats"].exists()


def test_run_vsearch_with_progress_tracks_searching_status_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events: list[tuple[str, int, str] | tuple[str, int]] = []
    stream_checkpoints: list[list[tuple[str, int, str] | tuple[str, int]]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, message))

        def update(self) -> None:
            events.append(("update", self.finished))

        def close(self) -> None:
            events.append(("close", self.finished))

    class _Proc:
        def __init__(self) -> None:
            self.pid = 123
            self._rc = 0

        def wait(self) -> int:
            return self._rc

        def poll(self) -> int:
            return self._rc

    monkeypatch.setattr(denovo_module, "ProgressBar", _ProgressStub)
    monkeypatch.setattr(
        denovo_module,
        "_open_vsearch_process_with_stderr_stream",
        lambda cmd: (_Proc(), 99, True),
    )

    def _fake_stderr_chunks(stderr_fd):
        yield "Reading file /tmp/x.fa 100%\rMasking 100%\rSearching 1%\r"
        stream_checkpoints.append(list(events))
        yield "Searching 2%\rSearching 2%\rSearching 99%\r"
        stream_checkpoints.append(list(events))
        yield "Searching 100%\rMatching unique query sequences: 0 of 10 (0.00%)\r"

    monkeypatch.setattr(
        denovo_module, "_iter_vsearch_stderr_chunks", _fake_stderr_chunks
    )

    denovo_module._run_vsearch_with_progress(
        ["vsearch"], message="Across-sample clustering"
    )

    assert ("init", 100, "Across-sample clustering") in events
    assert ("update", 0) in events
    assert ("update", 1) in events
    assert ("update", 2) in events
    assert ("update", 99) in events
    assert ("update", 100) in events
    assert events.count(("update", 100)) == 1
    assert ("update", 1) in stream_checkpoints[0]
    assert ("update", 2) not in stream_checkpoints[0]
    assert ("update", 99) in stream_checkpoints[1]
    assert ("update", 100) not in stream_checkpoints[1]
    assert events[-1] == ("close", 100)


def test_vsearch_pairs_cleans_large_sample_files_after_derep_and_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outdir = tmp_path
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    r1.write_text("", encoding="utf-8")
    r2.write_text("", encoding="utf-8")

    merged = outdir / "sample.merged.fa"
    joined = outdir / "sample.joined.fa"
    unmerged_r1 = outdir / "sample.unmerged_R1.fq"
    unmerged_r2 = outdir / "sample.unmerged_R2.fq"
    cluster_fa = outdir / "sample.cluster.fa"
    derep_unsorted = outdir / "sample.derep.fa"
    derep = outdir / "sample.derep.sizesorted.fa"
    consensus = outdir / "sample.consensus.fa"
    clusters = outdir / "sample.clusters.tsv"
    summary_tsv = outdir / "sample.summary.tsv"

    def fake_run_pipeline(cmds):
        cmd = cmds[-1]
        if "--fastq_mergepairs" in cmd:
            merged.write_text(">sample;M1\nAAATTT\n", encoding="utf-8")
            unmerged_r1.write_text("@r1\nAAA\n+\nIII\n", encoding="utf-8")
            unmerged_r2.write_text("@r2\nTTT\n+\nIII\n", encoding="utf-8")
            return
        if "--fastq_join" in cmd:
            joined.write_text(f">sample;J1\nAAA{'N' * 24}TTT\n", encoding="utf-8")
            return
        if "--sortbylength" in cmd:
            assert derep_unsorted.exists()
            assert not merged.exists()
            assert not joined.exists()
            assert not unmerged_r1.exists()
            assert not unmerged_r2.exists()
            assert not cluster_fa.exists()
            derep.write_text(">sample;J1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--fastx_uniques" in cmd:
            assert cmd[2] == str(cluster_fa)
            derep_unsorted.write_text(">sample;J1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--cluster_fast" in cmd:
            assert not merged.exists()
            assert not joined.exists()
            assert not unmerged_r1.exists()
            assert not unmerged_r2.exists()
            assert not cluster_fa.exists()
            assert derep_unsorted.exists()
            assert derep.exists()
            consensus.write_text(
                ">centroid=sample;J1;size=5;seqs=1\nAAATTT\n", encoding="utf-8"
            )
            clusters.write_text(
                "S\t0\t6\t*\t*\t*\t*\t*\tsample;J1;size=5\t*\n", encoding="utf-8"
            )
            return
        raise AssertionError(cmd)

    def fake_build_sample_summary(sname, outdir, *, seed_to_meta=None, **kwargs):
        assert sname == "sample"
        assert outdir == tmp_path
        assert seed_to_meta == {
            "sample;J1": ("joined", 3),
            "sample;M1": ("merged", 6),
        }
        assert derep.exists()
        assert consensus.exists()
        assert clusters.exists()
        summary_tsv.write_text("sample\tcluster_id\nsample\t0\n", encoding="utf-8")
        return pd.DataFrame()

    monkeypatch.setattr(denovo_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        denovo_module, "build_sample_summary", fake_build_sample_summary
    )

    denovo_module.vsearch_pairs(
        sname="sample",
        r1=r1,
        r2=r2,
        outdir=outdir,
        min_derep_size=5,
        min_merge_overlap=20,
        min_length=35,
        max_merge_diffs=4,
        allow_reverse_complement=False,
        within_similarity=0.95,
        threads=1,
        keep_intermediates=False,
        paired=True,
    )

    assert summary_tsv.exists()
    assert not merged.exists()
    assert not joined.exists()
    assert not unmerged_r1.exists()
    assert not unmerged_r2.exists()
    assert not cluster_fa.exists()
    assert not derep_unsorted.exists()
    assert not derep.exists()
    assert not consensus.exists()
    assert not clusters.exists()


def test_vsearch_pairs_keeps_sample_intermediates_when_requested(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outdir = tmp_path
    r1 = tmp_path / "sample.fastq.gz"
    r1.write_text("", encoding="utf-8")

    joined = outdir / "sample.joined.fa"
    cluster_fa = outdir / "sample.cluster.fa"
    derep_unsorted = outdir / "sample.derep.fa"
    derep = outdir / "sample.derep.sizesorted.fa"
    consensus = outdir / "sample.consensus.fa"
    clusters = outdir / "sample.clusters.tsv"

    def fake_run_pipeline(cmds):
        cmd = cmds[-1]
        if "--fastx_subsample" in cmd:
            joined.write_text(">sample;S1\nAAATTT\n", encoding="utf-8")
            return
        if "--fastx_uniques" in cmd:
            derep_unsorted.write_text(">sample;S1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--sortbylength" in cmd:
            derep.write_text(">sample;S1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--cluster_fast" in cmd:
            consensus.write_text(
                ">centroid=sample;S1;size=5;seqs=1\nAAATTT\n", encoding="utf-8"
            )
            clusters.write_text(
                "S\t0\t6\t*\t*\t*\t*\t*\tsample;S1;size=5\t*\n", encoding="utf-8"
            )
            return
        raise AssertionError(cmd)

    def fake_build_sample_summary(sname, outdir, *, seed_to_meta=None, **kwargs):
        (outdir / f"{sname}.summary.tsv").write_text(
            "sample\tcluster_id\nsample\t0\n", encoding="utf-8"
        )
        return pd.DataFrame()

    monkeypatch.setattr(denovo_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        denovo_module, "build_sample_summary", fake_build_sample_summary
    )

    denovo_module.vsearch_pairs(
        sname="sample",
        r1=r1,
        r2=None,
        outdir=outdir,
        min_derep_size=5,
        min_merge_overlap=20,
        min_length=35,
        max_merge_diffs=4,
        allow_reverse_complement=False,
        within_similarity=0.95,
        threads=1,
        keep_intermediates=True,
        paired=False,
    )

    assert joined.exists()
    assert cluster_fa.exists()
    assert derep_unsorted.exists()
    assert derep.exists()
    assert consensus.exists()
    assert clusters.exists()


def test_vsearch_pairs_preserves_pre_derep_files_when_derep_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outdir = tmp_path
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    r1.write_text("", encoding="utf-8")
    r2.write_text("", encoding="utf-8")

    merged = outdir / "sample.merged.fa"
    joined = outdir / "sample.joined.fa"
    unmerged_r1 = outdir / "sample.unmerged_R1.fq"
    unmerged_r2 = outdir / "sample.unmerged_R2.fq"
    cluster_fa = outdir / "sample.cluster.fa"
    derep_unsorted = outdir / "sample.derep.fa"

    def fake_run_pipeline(cmds):
        cmd = cmds[-1]
        if "--fastq_mergepairs" in cmd:
            merged.write_text(">sample;M1\nAAATTT\n", encoding="utf-8")
            unmerged_r1.write_text("@r1\nAAA\n+\nIII\n", encoding="utf-8")
            unmerged_r2.write_text("@r2\nTTT\n+\nIII\n", encoding="utf-8")
            return
        if "--fastq_join" in cmd:
            joined.write_text(f">sample;J1\nAAA{'N' * 24}TTT\n", encoding="utf-8")
            return
        if "--fastx_uniques" in cmd:
            raise RuntimeError("derep failed")
        raise AssertionError(cmd)

    monkeypatch.setattr(denovo_module, "run_pipeline", fake_run_pipeline)

    with pytest.raises(RuntimeError, match="derep failed"):
        denovo_module.vsearch_pairs(
            sname="sample",
            r1=r1,
            r2=r2,
            outdir=outdir,
            min_derep_size=5,
            min_merge_overlap=20,
            min_length=35,
            max_merge_diffs=4,
            allow_reverse_complement=False,
            within_similarity=0.95,
            threads=1,
            keep_intermediates=False,
            paired=True,
        )

    assert merged.exists()
    assert joined.exists()
    assert unmerged_r1.exists()
    assert unmerged_r2.exists()
    assert cluster_fa.exists()
    assert not derep_unsorted.exists()


def test_vsearch_pairs_preserves_post_derep_files_when_summary_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outdir = tmp_path
    r1 = tmp_path / "sample.fastq.gz"
    r1.write_text("", encoding="utf-8")

    joined = outdir / "sample.joined.fa"
    derep_unsorted = outdir / "sample.derep.fa"
    derep = outdir / "sample.derep.sizesorted.fa"
    consensus = outdir / "sample.consensus.fa"
    clusters = outdir / "sample.clusters.tsv"

    def fake_run_pipeline(cmds):
        cmd = cmds[-1]
        if "--fastx_subsample" in cmd:
            joined.write_text(">sample;S1\nAAATTT\n", encoding="utf-8")
            return
        if "--fastx_uniques" in cmd:
            derep_unsorted.write_text(">sample;S1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--sortbylength" in cmd:
            derep.write_text(">sample;S1;size=5\nAAATTT\n", encoding="utf-8")
            return
        if "--cluster_fast" in cmd:
            consensus.write_text(
                ">centroid=sample;S1;size=5;seqs=1\nAAATTT\n", encoding="utf-8"
            )
            clusters.write_text(
                "S\t0\t6\t*\t*\t*\t*\t*\tsample;S1;size=5\t*\n", encoding="utf-8"
            )
            return
        raise AssertionError(cmd)

    monkeypatch.setattr(denovo_module, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(
        denovo_module,
        "build_sample_summary",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("summary failed")),
    )

    with pytest.raises(RuntimeError, match="summary failed"):
        denovo_module.vsearch_pairs(
            sname="sample",
            r1=r1,
            r2=None,
            outdir=outdir,
            min_derep_size=5,
            min_merge_overlap=20,
            min_length=35,
            max_merge_diffs=4,
            allow_reverse_complement=False,
            within_similarity=0.95,
            threads=1,
            keep_intermediates=False,
            paired=False,
        )

    assert not joined.exists()
    assert derep_unsorted.exists()
    assert derep.exists()
    assert consensus.exists()
    assert clusters.exists()


def test_select_denovo_samples_keeps_all_when_input_count_is_at_most_cap(
    tmp_path: Path,
) -> None:
    fastq_dict = {}
    for idx in range(10):
        path = tmp_path / f"sample_{idx}.fastq.gz"
        path.write_bytes(b"A" * (idx + 1))
        fastq_dict[f"sample_{idx}"] = (path, None)

    selected, mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=None,
        use_all_samples=False,
    )

    assert mode == "all"
    assert selected == fastq_dict


def test_select_denovo_samples_fills_from_next_largest_when_top_half_is_below_cap(
    tmp_path: Path,
) -> None:
    fastq_dict = {}
    for idx in range(11):
        path = tmp_path / f"sample_{idx}.fastq.gz"
        path.write_bytes(b"A" * (idx + 1))
        fastq_dict[f"sample_{idx}"] = (path, None)

    selected, mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=None,
        use_all_samples=False,
    )

    assert mode == "top-half-random"
    assert list(selected) == [f"sample_{idx}" for idx in range(10, 0, -1)]


def test_select_denovo_samples_prefers_top_half_by_input_size(tmp_path: Path) -> None:
    fastq_dict = {}
    for idx in range(30):
        path = tmp_path / f"sample_{idx}.fastq.gz"
        path.write_bytes(b"A" * (idx + 1))
        fastq_dict[f"sample_{idx}"] = (path, None)

    selected, mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=None,
        use_all_samples=False,
    )
    repeated, repeated_mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=None,
        use_all_samples=False,
    )
    eligible = [f"sample_{idx}" for idx in range(29, 14, -1)]
    expected = sorted(
        random.Random(0).sample(eligible, denovo_module.DEFAULT_MAX_DENOVO_SAMPLES)
    )

    assert mode == "top-half-random"
    assert repeated_mode == "top-half-random"
    assert len(selected) == denovo_module.DEFAULT_MAX_DENOVO_SAMPLES
    assert list(selected) == expected
    assert selected == repeated
    assert set(selected).issubset({f"sample_{idx}" for idx in range(15, 30)})


def test_select_denovo_samples_imap_uses_all_glob_matched_samples(tmp_path: Path) -> None:
    imap_path = tmp_path / "denovo.imap.tsv"
    imap_path.write_text(
        "a* pop_a\nb* pop_b\n",
        encoding="utf-8",
    )
    a1 = tmp_path / "a1.fastq.gz"
    a2 = tmp_path / "a2.fastq.gz"
    b1 = tmp_path / "b1.fastq.gz"
    b2 = tmp_path / "b2.fastq.gz"
    a1.write_bytes(b"A" * 10)
    a2.write_bytes(b"A" * 20)
    b1.write_bytes(b"A" * 40)
    b2.write_bytes(b"A" * 5)
    fastq_dict = {
        "a1": (a1, None),
        "a2": (a2, None),
        "b1": (b1, None),
        "b2": (b2, None),
    }

    selected, mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=imap_path,
        use_all_samples=False,
    )

    assert mode == "imap"
    assert list(selected) == ["a1", "a2", "b1", "b2"]


def test_select_denovo_samples_imap_does_not_cap_or_warn_when_more_than_default_max(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    warnings: list[tuple[object, ...]] = []

    class DummyLogger:
        def info(self, *args, **kwargs) -> None:
            pass

        def warning(self, *args, **kwargs) -> None:
            warnings.append(args)

    imap_path = tmp_path / "denovo.imap.tsv"
    fastq_dict = {}
    lines = []
    for idx in range(11):
        name = f"s{idx}"
        path = tmp_path / f"{name}.fastq.gz"
        path.write_bytes(b"A" * (idx + 1))
        fastq_dict[name] = (path, None)
        lines.append(f"{name} pop_{idx}\n")
    imap_path.write_text("".join(lines), encoding="utf-8")
    monkeypatch.setattr(denovo_module, "logger", DummyLogger())

    selected, mode = denovo_module._select_denovo_samples(
        fastq_dict,
        imap_path=imap_path,
        use_all_samples=False,
    )

    assert mode == "imap"
    assert len(selected) == 11
    assert warnings == []


def test_write_ordered_consensus_stream_to_file_flushes_in_mapping_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t10\t0\tjoined\t0\ts1;J1\n"
        "1\tlocus_1_1\tcontract_1_2\ts2\t5\t1\t10\t0\tjoined\t0\ts2;J1\n"
        "2\tlocus_2_1\tcontract_2_1\ts3\t6\t1\t11\t0\tjoined\t0\ts3;J2\n"
        "2\tlocus_2_1\tcontract_2_2\ts4\t6\t1\t11\t0\tjoined\t0\ts4;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAAAAAAAA",
            record_type="joined",
            n_reads=5,
            length=10,
            arm_boundary=10,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J1",
            cluster_sequence="AAAAAAAATA",
            record_type="joined",
            n_reads=5,
            length=10,
            arm_boundary=10,
        )
        + _summary_row(
            sample="s3",
            cluster_id=0,
            seed="s3;J2",
            cluster_sequence="CCCCCCCCCCC",
            record_type="joined",
            n_reads=6,
            length=11,
            arm_boundary=11,
        )
        + _summary_row(
            sample="s4",
            cluster_id=0,
            seed="s4;J2",
            cluster_sequence="CCCCCCCCCCA",
            record_type="joined",
            n_reads=6,
            length=11,
            arm_boundary=11,
        ),
        encoding="utf-8",
    )

    def fake_iter_threaded_alignment_results(
        jobs_iter,
        max_workers,
        heartbeat_s,
    ):
        jobs = list(jobs_iter)
        assert [key for key, _job in jobs] == [0, 1]
        assert max_workers == 2
        assert heartbeat_s == align_module.STALL_HEARTBEAT_SECONDS
        yield 1, (2, "locus_2_1", "CONS2", False)
        yield 0, (1, "locus_1_1", "CONS1", False)

    monkeypatch.setattr(
        align_module,
        "_iter_threaded_alignment_results",
        fake_iter_threaded_alignment_results,
    )
    summary = align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary="mafft",
        cores=2,
    )

    assert (
        out_fa.read_text(encoding="utf-8") == ">locus_1_1\nCONS1\n>locus_2_1\nCONS2\n"
    )
    assert summary.mafft_required_loci == 2
    assert summary.mafft_worker_processes == 2


def test_load_summary_records_supports_compact_arm_boundary_schema(
    tmp_path: Path,
) -> None:
    summary_tsv = tmp_path / "concat.summary.tsv"
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAAACCCCC",
            record_type="joined",
            n_reads=5,
            length=34,
            arm_boundary=5,
        ),
        encoding="utf-8",
    )

    records = align_module._load_summary_records(summary_tsv)

    assert records["s1;J1"].cluster_sequence == "AAAAACCCCC"
    assert records["s1;J1"].left_arm == "AAAAA"
    assert records["s1;J1"].right_arm == "CCCCC"


def test_load_summary_records_rejects_old_consensus_only_schema(tmp_path: Path) -> None:
    summary_tsv = tmp_path / "concat.summary.tsv"
    summary_tsv.write_text(
        "sample\tcluster_id\tseed\tlength\tn_unique\tn_reads\tmerged\tconsensus\n"
        "s1\t0\ts1;J1\t10\t1\t5\tFalse\tAAAAAAAAAA\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="missing required columns"):
        align_module._load_summary_records(summary_tsv)


def test_write_ordered_consensus_stream_to_file_no_alignment_shows_progress(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t10\t0\tsingle\t0\ts1;J1\n"
        "2\tlocus_2_1\tcontract_2_1\ts2\t6\t1\t11\t0\tsingle\t0\ts2;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAAAAAAAA",
            record_type="single",
            n_reads=5,
            length=10,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J2",
            cluster_sequence="CCCCCCCCCCC",
            record_type="single",
            n_reads=6,
            length=11,
        ),
        encoding="utf-8",
    )
    events: list[tuple[str, int, str] | tuple[str, int]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, message))

        def update(self) -> None:
            events.append(("update", self.finished))

        def close(self) -> None:
            events.append(("close", self.finished))

    monkeypatch.setattr(align_module, "ProgressBar", _ProgressStub)

    align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary="mafft",
        alignment_mode="none",
    )

    assert ("init", 2, "Writing loci - total jobs: 2") in events
    assert ("update", 0) in events
    assert ("update", 1) in events
    assert ("update", 2) in events
    assert events[-1] == ("close", 2)


def test_write_ordered_consensus_stream_to_file_interrupt_raises_system_exit_130(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t10\t0\tsingle\t0\ts1;J1\n"
        "1\tlocus_1_1\tcontract_1_2\ts2\t6\t1\t11\t0\tsingle\t0\ts2;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAAAAAAAA",
            record_type="single",
            n_reads=5,
            length=10,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J2",
            cluster_sequence="CCCCCCCCCCA",
            record_type="single",
            n_reads=6,
            length=11,
        ),
        encoding="utf-8",
    )
    events: list[tuple[str, int, str] | tuple[str, int]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, message))

        def update(self) -> None:
            events.append(("update", self.finished))

        def close(self) -> None:
            events.append(("close", self.finished))

    def fake_iter_threaded_alignment_results(*args, **kwargs):
        raise KeyboardInterrupt
        yield  # pragma: no cover

    monkeypatch.setattr(align_module, "ProgressBar", _ProgressStub)
    monkeypatch.setattr(
        align_module,
        "_iter_threaded_alignment_results",
        fake_iter_threaded_alignment_results,
    )

    with pytest.raises(SystemExit) as excinfo:
        align_module.write_ordered_consensus_stream_to_file(
            mapping_tsv=mapping_tsv,
            summary_tsv=summary_tsv,
            out_fa=out_fa,
            mafft_binary="mafft",
            cores=2,
        )

    assert excinfo.value.code == 130
    assert events[-1][0] == "close"


def test_worker_build_consensus_returns_single_record_without_alignment() -> None:
    locus_id, locus_name, consensus, uses_output_spacer = (
        align_module.worker_build_consensus(
            locus_id=4,
            record=[("seed", "ACGTNN")],
            mafft_binary="mafft",
            min_prop=0.5,
            threads=1,
        )
    )

    assert locus_id == 4
    assert locus_name == "locus_4"
    assert consensus == "ACGTNN"
    assert uses_output_spacer is False


def test_consensus_from_aligned_uses_majority_base_even_below_min_prop() -> None:
    consensus = align_module.consensus_from_aligned(
        [("a", "A"), ("b", "C"), ("c", "C")],
        min_prop=0.9,
    )

    assert consensus == "C"


def test_consensus_from_aligned_breaks_exact_ties_by_acgt_order() -> None:
    consensus = align_module.consensus_from_aligned(
        [("a", "T"), ("b", "G"), ("c", "C"), ("d", "A")],
        min_prop=0.9,
    )

    assert consensus == "A"


def test_worker_build_consensus_identical_sequences_skip_mafft(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_mafft(*args, **kwargs):
        raise AssertionError("mafft should not be called for identical sequences")

    monkeypatch.setattr(align_module, "mafft_align_one", fail_mafft)

    locus_id, locus_name, consensus, uses_output_spacer = (
        align_module.worker_build_consensus(
            locus_id=7,
            record=[("a", "ACGT"), ("b", "ACGT"), ("c", "ACGT")],
            mafft_binary="mafft",
            min_prop=0.5,
            threads=1,
        )
    )

    assert locus_id == 7
    assert locus_name == "locus_7"
    assert consensus == "ACGT"
    assert uses_output_spacer is False


def test_worker_build_consensus_same_length_nonidentical_still_calls_mafft(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[object] = []

    def fake_mafft(record, mafft_binary, threads, **kwargs):
        calls.append((record, mafft_binary, threads, kwargs))
        return [("a", "ACGT"), ("b", "ACGA")]

    monkeypatch.setattr(align_module, "mafft_align_one", fake_mafft)

    locus_id, locus_name, consensus, uses_output_spacer = (
        align_module.worker_build_consensus(
            locus_id=8,
            record=[("a", "ACGT"), ("b", "ACGA")],
            mafft_binary="mafft",
            min_prop=0.5,
            threads=1,
        )
    )

    assert locus_id == 8
    assert locus_name == "locus_8"
    assert consensus == "ACGA"
    assert uses_output_spacer is False
    assert len(calls) == 1


def test_mafft_align_one_reports_locus_details_on_timeout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run_pipeline(*args, **kwargs):
        raise PipelineTimeoutError("pipeline timed out")

    monkeypatch.setattr(align_module, "run_pipeline", fake_run_pipeline)

    with pytest.raises(RuntimeError, match="mafft timed out for locus_12 after 15s"):
        align_module.mafft_align_one(
            [("a", "ACGT"), ("b", "ACGTAA")],
            mafft_binary="mafft",
            threads=2,
            locus_id=12,
            timeout_s=15,
        )


def test_write_ordered_consensus_stream_to_file_no_alignment_uses_longest_stripped_sequence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t7\t0\tjoined\t0\ts1;J1\n"
        "1\tlocus_1_1\tcontract_1_2\ts2\t6\t1\t30\t0\tjoined\t0\ts2;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAAAAA",
            record_type="joined",
            n_reads=5,
            length=7,
            arm_boundary=7,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J2",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        ),
        encoding="utf-8",
    )

    def fail_pool(*args, **kwargs):
        raise AssertionError(
            "threaded aligner should not be used for --no-alignment path"
        )

    monkeypatch.setattr(align_module, "_iter_threaded_alignment_results", fail_pool)
    align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary="mafft",
        alignment_mode="none",
    )

    assert out_fa.read_text(encoding="utf-8") == ">locus_1_1\nAAAAAAA\n"


def test_write_ordered_consensus_stream_to_file_joined_only_locus_writes_output_spacer(
    tmp_path: Path,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t30\t0\tjoined\t0\ts1;J1\n"
        "1\tlocus_1_1\tcontract_1_2\ts2\t6\t1\t30\t0\tjoined\t0\ts2;J2\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J2",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        ),
        encoding="utf-8",
    )

    summary = align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary="mafft",
        alignment_mode="none",
    )

    assert out_fa.read_text(encoding="utf-8") == f">locus_1_1\nAAA{'N' * 50}TTT\n"
    assert summary.joined_spacer_loci == 1
    assert summary.mixed_reconciled_spacer_loci == 0
    assert summary.stripped_output_loci == 0
    assert summary.output_spacer_length == 50


def test_write_ordered_consensus_stream_to_file_mixed_reconciled_locus_writes_output_spacer(
    tmp_path: Path,
) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"
    out_fa = tmp_path / "denovo_reference.fa"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\treconcile_mode\treconciled_group\toutput_form\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t30\t0\tjoined\t0\ts1;J1\tmixed\tcontract_1_1\tspaced\n"
        "1\tlocus_1_1\tcontract_1_2\ts2\t6\t1\t6\t1\tmerged\t0\ts2;M1\tmixed\tcontract_1_2\tspaced\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;M1",
            cluster_sequence="AAATTT",
            record_type="merged",
            n_reads=6,
            length=6,
        ),
        encoding="utf-8",
    )

    summary = align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=mapping_tsv,
        summary_tsv=summary_tsv,
        out_fa=out_fa,
        mafft_binary="mafft",
        alignment_mode="none",
    )

    assert out_fa.read_text(encoding="utf-8") == f">locus_1_1\nAAA{'N' * 50}TTT\n"
    assert summary.joined_spacer_loci == 0
    assert summary.mixed_reconciled_spacer_loci == 1
    assert summary.stripped_output_loci == 0


def test_iter_alignment_jobs_yields_only_mafft_required_loci(tmp_path: Path) -> None:
    mapping_tsv = tmp_path / "denovo.loci.mapping.tsv"
    summary_tsv = tmp_path / "concat.summary.tsv"

    mapping_tsv.write_text(
        "locus\tlocus_name\tcontract_group\tsample\tn_reads\tn_unique\tlength\tmerged\trecord_type\tcluster_id\tcore\n"
        "1\tlocus_1_1\tcontract_1_1\ts1\t5\t1\t4\t0\tsingle\t0\ts1;J1\n"
        "2\tlocus_2_1\tcontract_2_1\ts2\t5\t1\t4\t0\tsingle\t0\ts2;J2\n"
        "2\tlocus_2_1\tcontract_2_2\ts3\t5\t1\t4\t0\tsingle\t0\ts3;J2\n"
        "3\tlocus_3_1\tcontract_3_1\ts4\t5\t1\t4\t0\tsingle\t0\ts4;J3\n"
        "3\tlocus_3_1\tcontract_3_2\ts5\t5\t1\t4\t0\tsingle\t0\ts5;J3\n",
        encoding="utf-8",
    )
    summary_tsv.write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAAA",
            record_type="single",
            n_reads=5,
            length=4,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J2",
            cluster_sequence="CCCC",
            record_type="single",
            n_reads=5,
            length=4,
        )
        + _summary_row(
            sample="s3",
            cluster_id=0,
            seed="s3;J2",
            cluster_sequence="CCCA",
            record_type="single",
            n_reads=5,
            length=4,
        )
        + _summary_row(
            sample="s4",
            cluster_id=0,
            seed="s4;J3",
            cluster_sequence="GGGG",
            record_type="single",
            n_reads=5,
            length=4,
        )
        + _summary_row(
            sample="s5",
            cluster_id=0,
            seed="s5;J3",
            cluster_sequence="GGGA",
            record_type="single",
            n_reads=5,
            length=4,
        ),
        encoding="utf-8",
    )

    jobs = list(
        align_module.iter_alignment_jobs(
            mapping_tsv=mapping_tsv,
            summary_tsv=summary_tsv,
            mafft_binary="mafft",
            min_prop=0.5,
            threads=1,
            alignment_mode="mafft",
        )
    )

    assert [key for key, _job in jobs] == [1, 2]
    assert jobs[0][1]["locus_id"] == 2
    assert jobs[1][1]["locus_id"] == 3
    assert jobs[0][1]["threads"] == 1


def test_run_denovo_no_alignment_passes_alignment_mode_none(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sample_fastq = tmp_path / "sample.fastq.gz"
    sample_fastq.write_text("", encoding="utf-8")
    outdir = tmp_path / "OUT"
    calls: dict[str, object] = {}
    _patch_required_binaries(monkeypatch, tmp_path)

    monkeypatch.setattr(
        denovo_module,
        "get_name_to_fastq_dict",
        lambda *args, **kwargs: {"sample_a": (sample_fastq, None)},
    )
    monkeypatch.setattr(
        denovo_module, "run_with_pool", lambda *args, **kwargs: {"sample_a": None}
    )
    monkeypatch.setattr(
        denovo_module,
        "build_sample_summary",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError(
                "run_denovo should not build sample summaries in the main process"
            )
        ),
    )
    monkeypatch.setattr(
        denovo_module,
        "concat_summaries",
        lambda outdir: (outdir / "concat.summary.tsv").write_text(
            _summary_header()
            + _summary_row(
                sample="sample_a",
                cluster_id=0,
                seed="sample_a;S1",
                cluster_sequence="ACGTACGTAA",
                record_type="single",
                n_reads=5,
                length=10,
            ),
            encoding="utf-8",
        ),
    )
    monkeypatch.setattr(
        denovo_module, "vsearch_cluster_across", lambda *args, **kwargs: None
    )

    def fake_make_global_tables(
        outdir,
        cores,
        log_level,
        within_similarity,
    ):
        mapping = pd.DataFrame(
            [
                {
                    "locus": 1,
                    "locus_name": "locus_1_1",
                    "contract_group": "contract_1_1",
                    "core": "sample_a;S1",
                }
            ]
        )
        stats = pd.DataFrame(
            [
                {
                    "locus": 1,
                    "n_samples": 1,
                    "n_cores": 1,
                    "duplicated_component": False,
                    "used_reconciliation": False,
                }
            ]
        )
        mapping.to_csv(outdir.parent / "denovo.loci.mapping.tsv", sep="\t", index=False)
        stats.to_csv(outdir.parent / "denovo.loci.stats.tsv", sep="\t", index=False)
        return graph_module.GraphTableSummary(
            loci_written=int(stats.shape[0]),
            consensus_records=int(mapping.shape[0]),
        )

    def fake_write_ordered_consensus_stream_to_file(**kwargs):
        calls["consensus"] = kwargs
        kwargs["out_fa"].write_text(">locus_1_1\nACGTACGTAA\n", encoding="utf-8")
        return _alignment_summary(total_loci=1, alignment_mode="none")

    monkeypatch.setattr(denovo_module, "make_global_tables", fake_make_global_tables)
    monkeypatch.setattr(
        denovo_module,
        "write_ordered_consensus_stream_to_file",
        fake_write_ordered_consensus_stream_to_file,
    )

    denovo_module.run_denovo(
        fastqs=[sample_fastq],
        outdir=outdir,
        within_similarity=0.95,
        across_similarity=0.85,
        min_derep_size=2,
        min_length=35,
        min_merge_overlap=20,
        max_merge_diffs=4,
        delim_str=None,
        delim_idx=1,
        allow_reverse_complement=False,
        cores=6,
        threads=3,
        no_alignment=True,
        force=False,
        imap=None,
        use_all_samples=False,
        keep_intermediates=False,
        log_level="INFO",
    )

    assert calls["consensus"]["alignment_mode"] == "none"
    assert calls["consensus"]["cores"] == 6
    stats_text = (outdir / "denovo.stats.txt").read_text(encoding="utf-8")
    assert _report_has_value_line(stats_text, "alignment_mode", "none")
    assert _report_has_value_line(stats_text, "mafft_worker_processes", "0")


def test_write_denovo_stats_formats_assemble_style_sections(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "OUT"
    outdir.mkdir()
    outpath = outdir / "denovo.stats.txt"
    outputs = {
        "reference": outdir / "denovo_reference.fa",
        "mapping": outdir / "denovo.loci.mapping.tsv",
        "loci_stats": outdir / "denovo.loci.stats.tsv",
        "run_stats": outpath,
        "audit_dir": outdir / "denovo.audit",
    }

    denovo_module._write_denovo_stats(
        outpath,
        all_fastq_dict={
            "sample_a": (tmp_path / "sample_a.fastq.gz", None),
            "sample_b": (tmp_path / "sample_b.fastq.gz", None),
            "sample_c": (tmp_path / "sample_c.fastq.gz", None),
        },
        selected_fastq_dict={
            "sample_a": (tmp_path / "sample_a.fastq.gz", None),
            "sample_b": (tmp_path / "sample_b.fastq.gz", None),
        },
        selection_mode="top-half-random",
        paired=False,
        within_similarity=0.95,
        across_similarity=0.85,
        min_derep_size=5,
        min_length=35,
        min_merge_overlap=20,
        max_merge_diffs=4,
        allow_reverse_complement=False,
        cores=6,
        threads=3,
        workers=2,
        alignment_summary=_alignment_summary(
            total_loci=3,
            single_sequence_loci=1,
            identical_sequence_loci=1,
            mafft_required_loci=1,
            mafft_threads_per_job=1,
            mafft_worker_processes=2,
            alignment_mode="mafft",
            mafft_timeout_seconds=900,
            joined_spacer_loci=2,
            mixed_reconciled_spacer_loci=1,
            stripped_output_loci=1,
            output_spacer_length=50,
        ),
        keep_intermediates=False,
        workdir=outdir / denovo_module.WORKDIR_NAME,
        graph_summary=graph_module.GraphTableSummary(
            loci_written=3,
            consensus_records=5,
            duplicated_components_seen=2,
            same_sample_reconciliation_attempted=2,
            components_reconciled=1,
            joined_only_reconciled_loci=1,
            mixed_reconciled_loci=1,
            mixed_reconciled_groups=2,
        ),
        qc_summary=denovo_module.DenovoQcSummary(
            selected_sample_count=2,
            total_input_sample_count=3,
            consensus_records=5,
            loci_written=3,
            singleton_loci=2,
            singleton_locus_fraction=2 / 3,
            loci_with_2plus_samples=1,
            loci_with_half_or_more_selected_samples=3,
            loci_with_all_selected_samples=1,
            mean_samples_per_locus=4 / 3,
            median_samples_per_locus=1.0,
            mean_cores_per_locus=5 / 3,
            median_cores_per_locus=2.0,
            max_samples_per_locus=2,
            max_cores_per_locus=3,
            multi_core_single_sample_loci=1,
            duplicated_component_loci=2,
            reconciled_loci=1,
            audited_components=4,
            processed_components=3,
            oversize_unsplit_components=1,
            largest_component_nodes=11,
            component_input_nodes_p50=2.0,
            component_input_nodes_p90=7.5,
            component_input_nodes_p99=10.9,
            component_input_nodes_max=11,
            component_contracted_nodes_p50=1.0,
            component_contracted_nodes_p90=4.5,
            component_contracted_nodes_p99=5.9,
            component_contracted_nodes_max=6,
            occupancy_counts=((1, 2), (2, 1)),
            selected_sample_rows=(
                denovo_module.SelectedSampleQc(
                    sample="sample_b",
                    consensus_records=12,
                    n_reads_sum=1200,
                    joined_records=4,
                    merged_records=2,
                    single_records=0,
                ),
                denovo_module.SelectedSampleQc(
                    sample="sample_a",
                    consensus_records=3,
                    n_reads_sum=45,
                    joined_records=0,
                    merged_records=1,
                    single_records=2,
                ),
            ),
        ),
        outputs=outputs,
    )

    text = outpath.read_text(encoding="utf-8")

    assert "# Selected Sample Summary" in text
    assert "# Locus Occupancy" in text
    assert "# Component Node Summary" in text
    assert re.search(
        r"^sample\s+consensus_records\s+n_reads_sum\s+joined_records\s+merged_records\s+single_records$",
        text,
        re.MULTILINE,
    )
    assert re.search(r"^sample_b\s+12\s+1,200\s+4\s+2\s+0\s*$", text, re.MULTILINE)
    assert re.search(r"^sample_a\s+3\s+45\s+0\s+1\s+2\s*$", text, re.MULTILINE)
    assert re.search(
        r"^quantile\s+input_nodes\s+contracted_nodes$",
        text,
        re.MULTILINE,
    )
    assert re.search(r"^p90\s+7\.500\s+4\.500\s*$", text, re.MULTILINE)
    assert re.search(
        r"^samples_with_data\s+loci\s+fraction_of_final_loci$",
        text,
        re.MULTILINE,
    )
    assert re.search(r"^0\s+0\s+0\.000000\s*$", text, re.MULTILINE)
    assert re.search(r"^1\s+2\s+0\.666667\s*$", text, re.MULTILINE)
    assert re.search(r"^2\s+1\s+0\.333333\s*$", text, re.MULTILINE)
    assert _report_has_value_line(text, "selected_sample_count", "2")
    assert _report_has_value_line(text, "intermediates", "cleaned on success")


def test_collect_denovo_qc_summarizes_final_outputs(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "OUT"
    workdir = outdir / denovo_module.WORKDIR_NAME
    audit_dir = outdir / "denovo.audit"
    workdir.mkdir(parents=True)
    audit_dir.mkdir()
    outputs = {
        "reference": outdir / "denovo_reference.fa",
        "mapping": outdir / "denovo.loci.mapping.tsv",
        "loci_stats": outdir / "denovo.loci.stats.tsv",
        "run_stats": outdir / "denovo.stats.txt",
        "audit_dir": audit_dir,
        "workdir": workdir,
    }
    selected_fastq_dict = {
        "s1": (tmp_path / "s1.fastq.gz", None),
        "s2": (tmp_path / "s2.fastq.gz", None),
        "s3": (tmp_path / "s3.fastq.gz", None),
        "s4": (tmp_path / "s4.fastq.gz", None),
    }

    pd.DataFrame(
        [
            {
                "n_samples": 1,
                "n_cores": 1,
                "duplicated_component": False,
                "used_reconciliation": False,
            },
            {
                "n_samples": 2,
                "n_cores": 2,
                "duplicated_component": True,
                "used_reconciliation": True,
            },
            {
                "n_samples": 2,
                "n_cores": 3,
                "duplicated_component": True,
                "used_reconciliation": False,
            },
            {
                "n_samples": 1,
                "n_cores": 2,
                "duplicated_component": False,
                "used_reconciliation": False,
            },
        ]
    ).to_csv(outputs["loci_stats"], sep="\t", index=False)
    _write_audit_summary(
        audit_dir,
        [
            {
                "component_id": 1,
                "n_input_nodes": 3,
                "n_contracted_nodes": 2,
                "status": "processed",
            },
            {
                "component_id": 2,
                "n_input_nodes": 9,
                "n_contracted_nodes": 4,
                "status": "oversize_unsplit",
            },
        ],
    )
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;A",
            cluster_sequence="AAAAAA",
            record_type="joined",
            n_reads=10,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;B",
            cluster_sequence="CCCCCC",
            record_type="merged",
            n_reads=6,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;A",
            cluster_sequence="GGGGGG",
            record_type="single",
            n_reads=8,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=1,
            seed="s2;B",
            cluster_sequence="TTTTTT",
            record_type="single",
            n_reads=3,
            length=6,
        ),
        encoding="utf-8",
    )

    qc_summary = denovo_module._collect_denovo_qc(
        selected_fastq_dict=selected_fastq_dict,
        total_input_sample_count=8,
        graph_summary=graph_module.GraphTableSummary(
            loci_written=4,
            consensus_records=4,
        ),
        outputs=outputs,
    )

    assert qc_summary.selected_sample_count == 4
    assert qc_summary.singleton_loci == 2
    assert qc_summary.loci_with_2plus_samples == 2
    assert qc_summary.loci_with_half_or_more_selected_samples == 2
    assert qc_summary.loci_with_all_selected_samples == 0
    assert qc_summary.mean_samples_per_locus == pytest.approx(1.5)
    assert qc_summary.median_samples_per_locus == pytest.approx(1.5)
    assert qc_summary.multi_core_single_sample_loci == 1
    assert qc_summary.duplicated_component_loci == 2
    assert qc_summary.reconciled_loci == 1
    assert qc_summary.oversize_unsplit_components == 1
    assert qc_summary.largest_component_nodes == 9
    assert qc_summary.occupancy_counts == ((1, 2), (2, 2))
    assert qc_summary.selected_sample_rows[0].sample == "s1"
    assert qc_summary.selected_sample_rows[0].joined_records == 1
    assert qc_summary.selected_sample_rows[0].merged_records == 1
    assert qc_summary.selected_sample_rows[1].sample == "s2"
    assert qc_summary.selected_sample_rows[2].sample == "s3"
    assert qc_summary.selected_sample_rows[2].consensus_records == 0


def test_collect_denovo_qc_handles_empty_outputs(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "OUT"
    workdir = outdir / denovo_module.WORKDIR_NAME
    audit_dir = outdir / "denovo.audit"
    workdir.mkdir(parents=True)
    audit_dir.mkdir()
    outputs = {
        "reference": outdir / "denovo_reference.fa",
        "mapping": outdir / "denovo.loci.mapping.tsv",
        "loci_stats": outdir / "denovo.loci.stats.tsv",
        "run_stats": outdir / "denovo.stats.txt",
        "audit_dir": audit_dir,
        "workdir": workdir,
    }
    selected_fastq_dict = {"sample_a": (tmp_path / "sample_a.fastq.gz", None)}

    pd.DataFrame(
        columns=["n_samples", "n_cores", "duplicated_component", "used_reconciliation"]
    ).to_csv(
        outputs["loci_stats"],
        sep="\t",
        index=False,
    )
    _write_audit_summary(audit_dir, [])
    (workdir / "concat.summary.tsv").write_text(_summary_header(), encoding="utf-8")

    qc_summary = denovo_module._collect_denovo_qc(
        selected_fastq_dict=selected_fastq_dict,
        total_input_sample_count=1,
        graph_summary=graph_module.GraphTableSummary(),
        outputs=outputs,
    )

    assert qc_summary.loci_written == 0
    assert qc_summary.singleton_loci == 0
    assert qc_summary.occupancy_counts == ()
    assert qc_summary.audited_components == 0
    assert qc_summary.selected_sample_rows == (
        denovo_module.SelectedSampleQc(
            sample="sample_a",
            consensus_records=0,
            n_reads_sum=0,
            joined_records=0,
            merged_records=0,
            single_records=0,
        ),
    )


def test_make_global_tables_contracts_joined_duplicates_and_writes_hierarchical_locus_names(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;J2",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=4,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s3",
            cluster_id=0,
            seed="s3;M1",
            cluster_sequence="GGGGGG",
            record_type="merged",
            n_reads=7,
            length=6,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;J1\ts1;J2\t100.0\t+\t100.0\t30\t30",
                "s1;J1\ts2;J1\t99.0\t+\t100.0\t30\t30",
                "s1;J2\ts2;J1\t99.0\t+\t100.0\t30\t30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert mapping_df["locus_name"].tolist()[:3] == [
        "locus_1_1",
        "locus_1_1",
        "locus_1_1",
    ]
    assert (
        mapping_df.loc[
            mapping_df["locus_name"] == "locus_1_1", "contract_group"
        ].nunique()
        == 2
    )
    assert set(mapping_df.loc[mapping_df["sample"] == "s1", "contract_group"]) == {
        "contract_1_1"
    }
    assert summary.consensus_records == 4
    assert summary.loci_written == 2
    assert stats_df.loc[0, "n_samples"] == 2
    assert stats_df.loc[0, "n_cores"] == 3
    assert stats_df.loc[0, "n_contracted_groups"] == 1
    assert "locus_2_1" in mapping_df["locus_name"].tolist()


def test_make_global_tables_reconciles_mixed_joined_and_merged_component(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;M1",
            cluster_sequence="AAATTT",
            record_type="merged",
            n_reads=4,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;J1\ts2;J1\t99.0\t+\t100.0\t30\t30",
                "s1;M1\ts2;J1\t99.0\t+\t100.0\t6\t30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    graph_summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_1_1", "locus_1_1"]
    assert graph_summary.components_reconciled == 1
    assert mapping_df["reconcile_mode"].tolist() == ["mixed", "mixed", "mixed"]
    assert mapping_df["output_form"].tolist() == ["spaced", "spaced", "spaced"]
    assert set(mapping_df.loc[mapping_df["sample"] == "s1", "contract_group"]) == {
        "contract_1_1"
    }
    assert bool(stats_df.loc[0, "used_mixed_reconciliation"]) is True
    assert stats_df.loc[0, "output_form"] == "spaced"
    assert stats_df.loc[0, "n_reconciled_groups"] == 1


def test_make_global_tables_rounds_mean_and_std_fields_to_3_decimals(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;S1",
            cluster_sequence="AAAAAAA",
            record_type="single",
            n_reads=1,
            length=7,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;M1",
            cluster_sequence="CCCCCCCC",
            record_type="merged",
            n_reads=2,
            length=8,
        )
        + _summary_row(
            sample="s3",
            cluster_id=0,
            seed="s3;S1",
            cluster_sequence="GGGGGGGGGG",
            record_type="single",
            n_reads=4,
            length=10,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;S1\ts2;M1\t99.0\t+\t100.0\t7\t8",
                "s1;S1\ts3;S1\t99.0\t+\t100.0\t7\t10",
                "s2;M1\ts3;S1\t99.0\t+\t100.0\t8\t10",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    graph_module.make_global_tables(workdir)

    with open(
        workdir.parent / "denovo.loci.stats.tsv", "rt", encoding="utf-8", newline=""
    ) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))

    assert row["n_reads_mean"] == "2.333"
    assert row["n_reads_std"] == "1.528"
    assert row["length_mean"] == "8.333"
    assert row["length_std"] == "1.528"
    assert row["merged_freq"] == "0.333"


def test_make_global_tables_reconciles_joined_only_identical_duplicates_without_direct_same_sample_edge(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;J2",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=4,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;J1\ts2;J1\t99.0\t+\t100.0\t30\t30",
                "s1;J2\ts2;J1\t99.0\t+\t100.0\t30\t30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    graph_summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_1_1", "locus_1_1"]
    assert graph_summary.joined_only_reconciled_loci == 1
    assert mapping_df["reconcile_mode"].tolist() == [
        "joined_only",
        "joined_only",
        "joined_only",
    ]
    assert mapping_df["output_form"].tolist() == ["spaced", "spaced", "spaced"]
    assert set(mapping_df.loc[mapping_df["sample"] == "s1", "contract_group"]) == {
        "contract_1_1"
    }
    assert bool(stats_df.loc[0, "used_reconciliation"]) is True
    assert bool(stats_df.loc[0, "used_joined_only_reconciliation"]) is True
    assert stats_df.loc[0, "reconcile_mode"] == "joined_only"
    assert stats_df.loc[0, "output_form"] == "spaced"

    audit_dir = tmp_path / "denovo.audit"
    assert (audit_dir / "components.summary.tsv").exists()
    assert (audit_dir / "component_1.members.tsv").exists()
    assert (audit_dir / "component_1.fa").exists()
    audit_summary = pd.read_csv(audit_dir / "components.summary.tsv", sep="\t")
    assert audit_summary.loc[0, "reconcile_mode"] == "joined_only"
    assert int(audit_summary.loc[0, "n_final_loci"]) == 1


def test_make_global_tables_parallelizes_multiple_components(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;S1",
            cluster_sequence="AAAAAA",
            record_type="single",
            n_reads=5,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;S1",
            cluster_sequence="CCCCCC",
            record_type="single",
            n_reads=6,
            length=6,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text("", encoding="utf-8")

    calls: dict[str, object] = {}

    def fake_run_with_pool_iter(
        jobs_iter,
        log_level,
        max_workers,
        max_inflight=None,
        msg=None,
        njobs=None,
        progress_increment=None,
    ):
        jobs = list(jobs_iter)
        calls["count"] = len(jobs)
        calls["log_level"] = log_level
        calls["max_workers"] = max_workers
        calls["msg"] = msg
        calls["njobs"] = njobs
        for key, (func, kwargs) in jobs:
            yield key, func(**kwargs)

    monkeypatch.setattr(graph_module, "run_with_pool_iter", fake_run_with_pool_iter)

    summary = graph_module.make_global_tables(
        workdir,
        cores=2,
        log_level="DEBUG",
    )
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert calls == {
        "count": 2,
        "log_level": "DEBUG",
        "max_workers": 2,
        "msg": "Splitting global clusters",
        "njobs": 2,
    }
    assert summary.loci_written == 2
    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_2_1"]
    assert stats_df["locus"].tolist() == [1, 2]


def test_make_global_tables_parallel_path_matches_serial(
    tmp_path: Path,
) -> None:
    def _write_fixture(workdir: Path) -> None:
        workdir.mkdir()
        (workdir / "concat.summary.tsv").write_text(
            _summary_header()
            + _summary_row(
                sample="s1",
                cluster_id=0,
                seed="s1;J1",
                cluster_sequence="AAATTT",
                record_type="joined",
                n_reads=5,
                length=30,
                arm_boundary=3,
            )
            + _summary_row(
                sample="s1",
                cluster_id=1,
                seed="s1;J2",
                cluster_sequence="AAATTT",
                record_type="joined",
                n_reads=4,
                length=30,
                arm_boundary=3,
            )
            + _summary_row(
                sample="s2",
                cluster_id=0,
                seed="s2;J1",
                cluster_sequence="AAATTT",
                record_type="joined",
                n_reads=6,
                length=30,
                arm_boundary=3,
            )
            + _summary_row(
                sample="s3",
                cluster_id=0,
                seed="s3;M1",
                cluster_sequence="GGGGGG",
                record_type="merged",
                n_reads=7,
                length=6,
            ),
            encoding="utf-8",
        )
        (workdir / "global_hits.uc.tsv").write_text(
            "\n".join(
                [
                    "s1;J1\ts1;J2\t100.0\t+\t100.0\t30\t30",
                    "s1;J1\ts2;J1\t99.0\t+\t100.0\t30\t30",
                    "s1;J2\ts2;J1\t99.0\t+\t100.0\t30\t30",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    serial_workdir = tmp_path / "serial"
    parallel_workdir = tmp_path / "parallel"
    _write_fixture(serial_workdir)
    _write_fixture(parallel_workdir)

    serial_summary = graph_module.make_global_tables(
        serial_workdir,
        cores=1,
        log_level="INFO",
    )
    parallel_summary = graph_module.make_global_tables(
        parallel_workdir,
        cores=2,
        log_level="INFO",
    )
    serial_mapping, serial_stats = _read_graph_output_tables(serial_workdir)
    parallel_mapping, parallel_stats = _read_graph_output_tables(parallel_workdir)

    pd.testing.assert_frame_equal(serial_mapping, parallel_mapping)
    pd.testing.assert_frame_equal(serial_stats, parallel_stats)
    assert serial_summary == parallel_summary


def test_make_global_tables_writes_audits_in_completion_order_but_keeps_final_outputs_ordered(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;A1",
            cluster_sequence="AAAAAA",
            record_type="single",
            n_reads=5,
            length=6,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;A2",
            cluster_sequence="AAAAAA",
            record_type="single",
            n_reads=4,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=2,
            seed="s2;A1",
            cluster_sequence="CCCCCC",
            record_type="single",
            n_reads=6,
            length=6,
        )
        + _summary_row(
            sample="s3",
            cluster_id=3,
            seed="s3;B1",
            cluster_sequence="GGGGGG",
            record_type="single",
            n_reads=7,
            length=6,
        )
        + _summary_row(
            sample="s3",
            cluster_id=4,
            seed="s3;B2",
            cluster_sequence="GGGGGG",
            record_type="single",
            n_reads=3,
            length=6,
        )
        + _summary_row(
            sample="s4",
            cluster_id=5,
            seed="s4;B1",
            cluster_sequence="TTTTTT",
            record_type="single",
            n_reads=8,
            length=6,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;A1\ts2;A1\t99.0\t+\t100.0\t6\t6",
                "s1;A2\ts2;A1\t99.0\t+\t100.0\t6\t6",
                "s3;B1\ts4;B1\t99.0\t+\t100.0\t6\t6",
                "s3;B2\ts4;B1\t99.0\t+\t100.0\t6\t6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    def fake_run_with_pool_iter(
        jobs_iter,
        log_level,
        max_workers,
        max_inflight=None,
        msg=None,
        njobs=None,
        progress_increment=None,
    ):
        del log_level, max_workers, max_inflight, msg, njobs, progress_increment
        jobs = list(jobs_iter)
        results = [(key, func(**kwargs)) for key, (func, kwargs) in jobs]
        yield results[1]
        yield results[0]

    monkeypatch.setattr(graph_module, "run_with_pool_iter", fake_run_with_pool_iter)

    summary = graph_module.make_global_tables(
        workdir,
        cores=2,
    )
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    audit_summary = pd.read_csv(
        tmp_path / "denovo.audit" / "components.summary.tsv", sep="\t"
    )
    assert audit_summary["component_id"].tolist() == [2, 1]
    assert summary.loci_written == 2
    assert mapping_df["locus_name"].tolist() == ["locus_1_1"] * 3 + ["locus_2_1"] * 3
    assert stats_df["component_id"].tolist() == [1, 2]


def test_make_global_tables_logs_when_flush_is_blocked_by_earlier_component(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;A1",
            cluster_sequence="AAAAAA",
            record_type="single",
            n_reads=5,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=1,
            seed="s2;A1",
            cluster_sequence="CCCCCC",
            record_type="single",
            n_reads=6,
            length=6,
        )
        + _summary_row(
            sample="s3",
            cluster_id=2,
            seed="s3;A1",
            cluster_sequence="GGGGGG",
            record_type="single",
            n_reads=7,
            length=6,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text("", encoding="utf-8")
    monotonic_values = iter([0.0, 1.0])
    warnings: list[str] = []

    def fake_run_with_pool_iter(
        jobs_iter,
        log_level,
        max_workers,
        max_inflight=None,
        msg=None,
        njobs=None,
        progress_increment=None,
    ):
        del log_level, max_workers, max_inflight, msg, njobs, progress_increment
        jobs = list(jobs_iter)
        results = [(key, func(**kwargs)) for key, (func, kwargs) in jobs]
        yield results[1]
        yield results[2]
        yield results[0]

    monkeypatch.setattr(graph_module, "run_with_pool_iter", fake_run_with_pool_iter)
    monkeypatch.setattr(graph_module, "FLUSH_STALL_LOG_SECONDS", 0.0)
    monkeypatch.setattr(graph_module.time, "monotonic", lambda: next(monotonic_values))
    monkeypatch.setattr(
        graph_module.logger, "warning", lambda message: warnings.append(str(message))
    )

    summary = graph_module.make_global_tables(
        workdir,
        cores=2,
    )
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert summary.loci_written == 3
    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_2_1", "locus_3_1"]
    assert stats_df["component_id"].tolist() == [1, 2, 3]
    assert any(
        "component 1" in message
        and "later component results already buffered" in message
        for message in warnings
    )


def test_make_global_tables_keeps_oversize_component_as_unsplit_locus_and_records_audit(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    header = _summary_header()
    summary_rows: list[str] = []
    component_one_seeds: list[str] = []
    cluster_id = 0
    total_samples = 6
    per_sample = 11
    for idx in range(total_samples * per_sample):
        sample = f"s{(idx % total_samples) + 1}"
        seed = f"{sample};BIG{idx}"
        component_one_seeds.append(seed)
        summary_rows.append(
            _summary_row(
                sample=sample,
                cluster_id=cluster_id,
                seed=seed,
                cluster_sequence=_unique_test_sequence(idx),
                record_type="single",
                n_reads=5,
                length=8,
            )
        )
        cluster_id += 1

    component_two = ("s1;KEEP0", "s2;KEEP0")
    summary_rows.extend(
        [
            _summary_row(
                sample="s1",
                cluster_id=cluster_id,
                seed=component_two[0],
                cluster_sequence="CCCCCC",
                record_type="single",
                n_reads=6,
                length=6,
            ),
            _summary_row(
                sample="s2",
                cluster_id=cluster_id + 1,
                seed=component_two[1],
                cluster_sequence="GGGGGG",
                record_type="single",
                n_reads=7,
                length=6,
            ),
        ]
    )
    (workdir / "concat.summary.tsv").write_text(
        header + "".join(summary_rows), encoding="utf-8"
    )

    edge_rows = [
        f"{left}\t{right}\t99.0\t+\t100.0\t6\t6"
        for left, right in zip(component_one_seeds, component_one_seeds[1:])
    ]
    edge_rows.append(f"{component_two[0]}\t{component_two[1]}\t99.0\t+\t100.0\t6\t6")
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(edge_rows) + "\n", encoding="utf-8"
    )
    summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert (
        mapping_df["locus_name"].tolist()
        == ["locus_1_1"] * total_samples + ["locus_2_1"] * 2
    )
    assert summary.raw_oversize_placeholder_components == 1
    assert stats_df["locus"].tolist() == [1, 2]
    assert mapping_df.loc[
        mapping_df["locus_name"] == "locus_1_1", "sample"
    ].tolist() == [f"s{idx}" for idx in range(1, total_samples + 1)]
    assert (
        mapping_df.loc[mapping_df["locus_name"] == "locus_1_1", "output_form"].tolist()
        == ["stripped"] * total_samples
    )
    assert stats_df.loc[0, "component_id"] == 1
    assert stats_df.loc[0, "n_samples"] == total_samples
    assert stats_df.loc[0, "n_cores"] == total_samples
    assert bool(stats_df.loc[0, "same_sample_reconciliation_attempted"]) is False
    assert bool(stats_df.loc[0, "used_reconciliation"]) is False
    assert stats_df.loc[0, "output_form"] == "stripped"
    assert stats_df.loc[1, "component_id"] == 2

    audit_dir = tmp_path / "denovo.audit"
    audit_summary = pd.read_csv(audit_dir / "components.summary.tsv", sep="\t")
    assert audit_summary["component_id"].tolist() == [1]
    assert audit_summary.loc[0, "status"] == "oversize_unsplit"
    assert int(audit_summary.loc[0, "n_contracted_nodes"]) == total_samples
    assert bool(audit_summary.loc[0, "used_oversize_rescue"]) is False
    assert bool(audit_summary.loc[0, "used_residue_cleanup"]) is False
    assert (
        pd.isna(audit_summary.loc[0, "discard_reason"])
        or audit_summary.loc[0, "discard_reason"] == ""
    )
    assert int(audit_summary.loc[0, "discard_limit_nodes"]) == 60
    assert int(audit_summary.loc[0, "n_input_nodes"]) == total_samples * per_sample
    assert int(audit_summary.loc[0, "n_final_loci"]) == 1
    assert bool(audit_summary.loc[0, "same_sample_reconciliation_attempted"]) is False
    assert not (audit_dir / "component_1.members.tsv").exists()
    assert not (audit_dir / "component_1.fa").exists()


def test_make_global_tables_keeps_raw_oversize_component_as_placeholder_without_contraction(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    summary_rows: list[str] = []
    component_seeds: list[str] = []
    cluster_id = 0
    for sample in ("s1", "s2", "s3", "s4", "s5"):
        for idx in range(11):
            seed = f"{sample};RAW{idx}"
            component_seeds.append(seed)
            summary_rows.append(
                _summary_row(
                    sample=sample,
                    cluster_id=cluster_id,
                    seed=seed,
                    cluster_sequence="AAAACC",
                    record_type="single",
                    n_reads=5 + idx,
                    length=6,
                )
            )
            cluster_id += 1
    (workdir / "concat.summary.tsv").write_text(
        _summary_header() + "".join(summary_rows), encoding="utf-8"
    )
    edge_rows = [
        f"{left}\t{right}\t99.0\t+\t100.0\t6\t6"
        for left, right in zip(component_seeds, component_seeds[1:])
    ]
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(edge_rows) + "\n", encoding="utf-8"
    )
    calls = {"count": 0}
    original_build_contracted_component = graph_module._build_contracted_component

    def wrapped_build_contracted_component(*args, **kwargs):
        calls["count"] += 1
        return original_build_contracted_component(*args, **kwargs)

    monkeypatch.setattr(
        graph_module, "_build_contracted_component", wrapped_build_contracted_component
    )

    summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert calls == {"count": 0}
    assert summary.raw_oversize_placeholder_components == 1
    assert mapping_df["locus_name"].tolist() == ["locus_1_1"] * 5
    assert stats_df.loc[0, "component_id"] == 1
    assert stats_df.loc[0, "n_samples"] == 5
    assert stats_df.loc[0, "n_cores"] == 5
    assert bool(stats_df.loc[0, "same_sample_reconciliation_attempted"]) is False
    assert bool(stats_df.loc[0, "used_reconciliation"]) is False

    audit_summary = pd.read_csv(
        tmp_path / "denovo.audit" / "components.summary.tsv", sep="\t"
    )
    assert audit_summary.loc[0, "status"] == "oversize_unsplit"
    assert int(audit_summary.loc[0, "n_input_nodes"]) == 55
    assert int(audit_summary.loc[0, "n_contracted_nodes"]) == 5
    assert bool(audit_summary.loc[0, "same_sample_reconciliation_attempted"]) is False
    assert bool(audit_summary.loc[0, "used_oversize_rescue"]) is False


def test_make_global_tables_keeps_component_at_oversize_limit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    header = _summary_header()
    summary_rows: list[str] = []
    component_seeds: list[str] = []
    cluster_id = 0
    for sample in ("s1", "s2", "s3", "s4", "s5"):
        for idx in range(10):
            seed = f"{sample};BOUND{idx}"
            component_seeds.append(seed)
            summary_rows.append(
                _summary_row(
                    sample=sample,
                    cluster_id=cluster_id,
                    seed=seed,
                    cluster_sequence="AAAACC",
                    record_type="single",
                    n_reads=5,
                    length=6,
                )
            )
            cluster_id += 1
    (workdir / "concat.summary.tsv").write_text(
        header + "".join(summary_rows), encoding="utf-8"
    )

    edge_rows = [
        f"{left}\t{right}\t99.0\t+\t100.0\t6\t6"
        for left, right in zip(component_seeds, component_seeds[1:])
    ]
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(edge_rows) + "\n", encoding="utf-8"
    )

    calls: dict[str, object] = {}

    def fake_process_component(
        component_id,
        component,
        summary_records,
        component_edges,
        seed_order,
        within_similarity,
        max_component_nodes,
    ):
        del component_edges, seed_order, within_similarity, max_component_nodes
        calls["component_id"] = component_id
        calls["component_size"] = len(component)
        seed = component[0]
        record = summary_records[seed]
        merged = int(record.record_type == "merged")
        return component_id, graph_module.ComponentResult(
            component_id=component_id,
            parts=(
                graph_module.ComponentPart(
                    component_id=component_id,
                    subcomponent_id=1,
                    mapping_rows=(
                        {
                            "component_id": component_id,
                            "subcomponent_id": 1,
                            "contract_group": "contract_1_1",
                            "sample": record.sample,
                            "n_reads": record.n_reads,
                            "n_unique": record.n_unique,
                            "length": record.length,
                            "cluster_length": record.cluster_length,
                            "merged": merged,
                            "record_type": record.record_type,
                            "cluster_id": record.cluster_id,
                            "core": seed,
                            "reconcile_mode": "none",
                            "reconciled_group": "",
                            "output_form": "stripped",
                        },
                    ),
                    stats_row={
                        "component_id": component_id,
                        "subcomponent_id": 1,
                        "n_samples": 1,
                        "n_cores": 1,
                        "n_contracted_groups": 0,
                        "n_reconciled_groups": 0,
                        "n_mixed_records": 0,
                        "n_reads_sum": record.n_reads,
                        "n_reads_mean": float(record.n_reads),
                        "n_reads_std": float("nan"),
                        "length_mean": float(record.length),
                        "length_std": float("nan"),
                        "merged_freq": float(merged),
                        "duplicated_component": True,
                        "same_sample_reconciliation_attempted": False,
                        "used_reconciliation": False,
                        "reconcile_mode": "none",
                        "used_joined_only_reconciliation": False,
                        "used_mixed_reconciliation": False,
                        "output_form": "stripped",
                        "samples": record.sample,
                    },
                ),
            ),
            audit_summary={
                "component_id": component_id,
                "n_input_nodes": len(component),
                "n_contracted_nodes": len(component),
                "n_input_samples": 5,
                "n_duplicate_samples": 5,
                "has_joined": False,
                "has_merged": False,
                "same_sample_reconciliation_attempted": False,
                "used_reconciliation": False,
                "reconcile_mode": "none",
                "used_oversize_rescue": False,
                "used_residue_cleanup": False,
                "status": "processed",
                "discard_reason": "",
                "discard_limit_nodes": 0,
                "n_final_loci": 1,
            },
            audit_rows=(),
            audit_fasta=(),
        )

    monkeypatch.setattr(graph_module, "_process_component", fake_process_component)

    summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert calls == {"component_id": 1, "component_size": 50}
    assert summary.raw_oversize_placeholder_components == 0
    assert mapping_df["locus_name"].tolist() == ["locus_1_1"]
    assert stats_df.loc[0, "component_id"] == 1


def test_oversize_placeholder_locus_is_written_to_final_reference(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    header = _summary_header()
    summary_rows: list[str] = []
    component_seeds: list[str] = []
    cluster_id = 0
    expected_consensus = "TTTTGGGGCC"
    counts = {"s1": 11, "s2": 10, "s3": 10, "s4": 10, "s5": 10}
    sample_idx = {sample: 0 for sample in counts}
    while any(counts.values()):
        for sample in tuple(counts):
            if counts[sample] <= 0:
                continue
            idx = sample_idx[sample]
            seed = f"{sample};BIG{idx}"
            component_seeds.append(seed)
            cluster_sequence = _unique_test_sequence(cluster_id)
            n_reads = 5
            cluster_length = len(cluster_sequence)
            if sample == "s5" and idx == 0:
                cluster_sequence = expected_consensus
                n_reads = 9
                cluster_length = len(cluster_sequence)
            summary_rows.append(
                _summary_row(
                    sample=sample,
                    cluster_id=cluster_id,
                    seed=seed,
                    cluster_sequence=cluster_sequence,
                    record_type="single",
                    n_reads=n_reads,
                    length=cluster_length,
                    cluster_length=cluster_length,
                )
            )
            cluster_id += 1
            sample_idx[sample] += 1
            counts[sample] -= 1
    (workdir / "concat.summary.tsv").write_text(
        header + "".join(summary_rows), encoding="utf-8"
    )
    edge_rows = [
        f"{left}\t{right}\t99.0\t+\t100.0\t6\t6"
        for left, right in zip(component_seeds, component_seeds[1:])
    ]
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(edge_rows) + "\n", encoding="utf-8"
    )

    graph_summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)
    out_fa = tmp_path / "denovo_reference.fa"
    summary = align_module.write_ordered_consensus_stream_to_file(
        mapping_tsv=tmp_path / "denovo.loci.mapping.tsv",
        summary_tsv=workdir / "concat.summary.tsv",
        out_fa=out_fa,
        mafft_binary="mafft",
        alignment_mode="none",
    )

    assert graph_summary.raw_oversize_placeholder_components == 1
    assert mapping_df["locus_name"].tolist() == ["locus_1_1"] * 5
    assert stats_df["component_id"].tolist() == [1]
    assert out_fa.read_text(encoding="utf-8") == f">locus_1_1\n{expected_consensus}\n"
    assert summary.total_loci == 1


def test_constrained_split_component_recovers_clean_residue_components() -> None:
    nodes = {"a1", "a2", "b1", "c1"}
    edges = {
        ("a1", "b1"): (0.99, 1.0),
        ("a2", "b1"): (0.98, 1.0),
        ("a2", "c1"): (0.97, 1.0),
    }
    node_samples = {
        "a1": frozenset({"s1"}),
        "a2": frozenset({"s1"}),
        "b1": frozenset({"s2"}),
        "c1": frozenset({"s3"}),
    }
    node_order = {"a1": 0, "a2": 1, "b1": 2, "c1": 3}

    result = split_component_constrained(nodes, edges, node_samples, node_order)

    assert result == [{"a1", "b1"}, {"a2", "c1"}]


def test_constrained_split_component_cleans_duplicate_rich_chain() -> None:
    nodes = {"a1", "a2", "b1", "c1", "d1"}
    edges = {
        ("a1", "b1"): (0.99, 1.0),
        ("a2", "b1"): (0.98, 1.0),
        ("a2", "c1"): (0.97, 1.0),
        ("c1", "d1"): (0.96, 1.0),
    }
    node_samples = {
        "a1": frozenset({"s1"}),
        "a2": frozenset({"s1"}),
        "b1": frozenset({"s2"}),
        "c1": frozenset({"s3"}),
        "d1": frozenset({"s4"}),
    }
    node_order = {"a1": 0, "a2": 1, "b1": 2, "c1": 3, "d1": 4}

    constrained_result = split_component_constrained(
        nodes, edges, node_samples, node_order
    )

    constrained_clean_nodes = sum(
        len(component) for component in constrained_result if len(component) > 1
    )

    assert constrained_clean_nodes == 5


def test_make_global_tables_reconciles_same_sample_direct_edge_above_within_similarity(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;M1",
            cluster_sequence="AAATTT",
            record_type="merged",
            n_reads=5,
            length=6,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;M2",
            cluster_sequence="AAATTC",
            record_type="merged",
            n_reads=4,
            length=6,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;M1",
            cluster_sequence="AAATTT",
            record_type="merged",
            n_reads=6,
            length=6,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;M1\ts1;M2\t95.0\t+\t100.0\t6\t6",
                "s1;M1\ts2;M1\t99.0\t+\t100.0\t6\t6",
                "s1;M2\ts2;M1\t95.0\t+\t100.0\t6\t6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    graph_summary = graph_module.make_global_tables(workdir, within_similarity=0.95)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_1_1", "locus_1_1"]
    assert graph_summary.components_reconciled == 1
    assert set(mapping_df.loc[mapping_df["sample"] == "s1", "contract_group"]) == {
        "contract_1_1"
    }
    assert bool(stats_df.loc[0, "used_reconciliation"]) is True
    assert bool(stats_df.loc[0, "same_sample_reconciliation_attempted"]) is True


def test_make_global_tables_does_not_require_mafft_for_split_stage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "concat.summary.tsv").write_text(
        _summary_header()
        + _summary_row(
            sample="s1",
            cluster_id=0,
            seed="s1;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=5,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s1",
            cluster_id=1,
            seed="s1;J2",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=4,
            length=30,
            arm_boundary=3,
        )
        + _summary_row(
            sample="s2",
            cluster_id=0,
            seed="s2;J1",
            cluster_sequence="AAATTT",
            record_type="joined",
            n_reads=6,
            length=30,
            arm_boundary=3,
        ),
        encoding="utf-8",
    )
    (workdir / "global_hits.uc.tsv").write_text(
        "\n".join(
            [
                "s1;J1\ts2;J1\t99.0\t+\t100.0\t30\t30",
                "s1;J2\ts2;J1\t99.0\t+\t100.0\t30\t30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        graph_module,
        "mafft_align_one",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("graph splitting should not call MAFFT")
        ),
        raising=False,
    )

    graph_summary = graph_module.make_global_tables(workdir)
    mapping_df, stats_df = _read_graph_output_tables(workdir)

    assert mapping_df["locus_name"].tolist() == ["locus_1_1", "locus_1_1", "locus_1_1"]
    assert graph_summary.components_reconciled == 1
    assert bool(stats_df.loc[0, "used_reconciliation"]) is True
