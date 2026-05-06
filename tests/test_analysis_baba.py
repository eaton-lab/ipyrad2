import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from loguru import logger

from ipyrad2.analysis.methods.baba import run_baba_method
from ipyrad2.utils.exceptions import IPyradError


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _dose_triplet(dose: int) -> np.ndarray:
    if dose == 255:
        return np.array([255, 255, ord("N")], dtype=np.uint8)
    if dose == 0:
        return np.array([0, 0, ord("A")], dtype=np.uint8)
    if dose == 1:
        return np.array([0, 1, ord("R")], dtype=np.uint8)
    return np.array([1, 1, ord("G")], dtype=np.uint8)


def _write_baba_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["a1", "a2", "b1", "b2", "c1", "d1", "o1"]
    nsites = 28
    doses = np.zeros((len(names), nsites), dtype=np.uint8)

    # 0-13: ABBA for A/B/C/O, D ancestral.
    doses[0, 0:14] = 0
    doses[1, 0:14] = 0
    doses[2, 0:14] = 2
    doses[3, 0:14] = 2
    doses[4, 0:14] = 2
    doses[5, 0:14] = 0
    doses[6, 0:14] = 0

    # 14-18: BABA for A/B/C/O, D ancestral.
    doses[0, 14:19] = 2
    doses[1, 14:19] = 2
    doses[2, 14:19] = 0
    doses[3, 14:19] = 0
    doses[4, 14:19] = 2
    doses[5, 14:19] = 0
    doses[6, 14:19] = 0

    # 19-23: BBAA for A/B/C/O, D ancestral.
    doses[0, 19:24] = 2
    doses[1, 19:24] = 2
    doses[2, 19:24] = 2
    doses[3, 19:24] = 2
    doses[4, 19:24] = 0
    doses[5, 19:24] = 0
    doses[6, 19:24] = 0

    # 24-27: ABBA again, but one A sample is missing to exercise quartet minmap.
    doses[0, 24:28] = 0
    doses[1, 24:28] = 255
    doses[2, 24:28] = 2
    doses[3, 24:28] = 2
    doses[4, 24:28] = 2
    doses[5, 24:28] = 0
    doses[6, 24:28] = 0

    genos = np.zeros((len(names), nsites, 3), dtype=np.uint8)
    for ridx in range(len(names)):
        for cidx in range(nsites):
            genos[ridx, cidx] = _dose_triplet(int(doses[ridx, cidx]))

    snpsmap = np.array(
        [
            [idx, 0, 0, 0 if idx < 14 else 1, 100 + idx]
            for idx in range(nsites)
        ],
        dtype=np.uint32,
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["nsnps"] = int(nsites)
        io5.create_dataset("genos", data=genos)
        io5["genos"].attrs["names"] = np.array(names, dtype=string_dtype)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A")] * nsites, dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(
            ["loc", "loc_idx", "loc_pos", "scaff", "pos"],
            dtype=string_dtype,
        )
    return path


def _write_imap(path: Path) -> Path:
    path.write_text(
        "\n".join(
            [
                "a1\tA",
                "a2\tA",
                "b1\tB",
                "b2\tB",
                "c1\tC",
                "d1\tD",
                "o1\tO",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def _write_minmap(path: Path) -> Path:
    path.write_text(
        "\n".join(
            [
                "A\t2",
                "B\t2",
                "C\t1",
                "D\t1",
                "O\t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def _write_rooted_order_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["p1", "p2", "p3", "o1"]
    nsites = 18
    doses = np.zeros((len(names), nsites), dtype=np.uint8)

    # 0-9: ABBA
    doses[0, 0:10] = 0
    doses[1, 0:10] = 2
    doses[2, 0:10] = 2
    doses[3, 0:10] = 0

    # 10-11: BABA
    doses[0, 10:12] = 2
    doses[1, 10:12] = 0
    doses[2, 10:12] = 2
    doses[3, 10:12] = 0

    # 12-17: BBAA
    doses[0, 12:18] = 2
    doses[1, 12:18] = 2
    doses[2, 12:18] = 0
    doses[3, 12:18] = 0

    genos = np.zeros((len(names), nsites, 3), dtype=np.uint8)
    for ridx in range(len(names)):
        for cidx in range(nsites):
            genos[ridx, cidx] = _dose_triplet(int(doses[ridx, cidx]))

    snpsmap = np.array(
        [[idx, 0, 0, 0, 100 + idx] for idx in range(nsites)],
        dtype=np.uint32,
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["nsnps"] = int(nsites)
        io5.create_dataset("genos", data=genos)
        io5["genos"].attrs["names"] = np.array(names, dtype=string_dtype)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A")] * nsites, dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(
            ["loc", "loc_idx", "loc_pos", "scaff", "pos"],
            dtype=string_dtype,
        )
    return path


def test_run_baba_method_manual_population_quartet_respects_quartet_minmap(
    tmp_path: Path,
) -> None:
    h5 = _write_baba_h5(tmp_path / "snps.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")
    minmap = _write_minmap(tmp_path / "minmap.tsv")
    tests = tmp_path / "quartets.tsv"
    tests.write_text("A B C O\n", encoding="utf-8")

    run_baba_method(
        data=h5,
        name="baba",
        outdir=tmp_path / "OUT",
        tests=tests,
        tree=None,
        imap=imap,
        minmap=minmap,
        min_sample_coverage=1,
        exclude=None,
        include_reference=False,
        resampling="none",
        bootstrap_replicates=10,
        jackknife_block_bp=1,
        jackknife_block_loci=2,
        seed=3,
        f_branch=False,
        f_branch_p_threshold=0.01,
        write_block_table=False,
        clustering_stats=False,
        cores=1,
        force=True,
        log_level="INFO",
        logged_command="ipyrad2 baba -d snps.hdf5 -o OUT --tests quartets.tsv -i imap.tsv -g minmap.tsv",
    )

    quartets = pd.read_csv(tmp_path / "OUT" / "baba.quartets.tsv", sep="\t")
    rooted = pd.read_csv(tmp_path / "OUT" / "baba.rooted.tsv", sep="\t")
    resolved = pd.read_csv(tmp_path / "OUT" / "baba.tests.resolved.tsv", sep="\t")
    manifest = (tmp_path / "OUT" / "baba.manifest.txt").read_text(encoding="utf-8")
    summary_json = _read_json(tmp_path / "OUT" / "baba.summary.json")

    assert len(quartets.index) == 1
    assert quartets.loc[0, "p1"] == "A"
    assert quartets.loc[0, "p2"] == "B"
    assert quartets.loc[0, "p3"] == "C"
    assert quartets.loc[0, "p4"] == "O"
    assert quartets.loc[0, "n_sites_tested"] == 24
    assert quartets.loc[0, "abba"] == pytest.approx(14.0)
    assert quartets.loc[0, "baba"] == pytest.approx(5.0)
    assert quartets.loc[0, "bbaa"] == pytest.approx(5.0)
    assert quartets.loc[0, "d_stat"] == pytest.approx((14.0 - 5.0) / (14.0 + 5.0))
    assert quartets.loc[0, "resampling_mode"] == "none"
    input_row = rooted[rooted["orientation"] == "input"].iloc[0]
    assert input_row["p1"] == "A"
    assert input_row["p2"] == "B"
    assert input_row["p3"] == "C"
    assert input_row["d_stat"] == pytest.approx((14.0 - 5.0) / (14.0 + 5.0))
    assert input_row["f_g"] == pytest.approx((14.0 - 5.0) / 14.0)
    assert resolved.to_dict("records") == [
        {"source": "tests", "p1": "A", "p2": "B", "p3": "C", "p4": "O"}
    ]
    assert manifest.startswith("CMD: ipyrad2 baba -d snps.hdf5 -o OUT --tests quartets.tsv -i imap.tsv -g minmap.tsv")
    assert "rooted_rows_written: 3" in manifest
    assert "quartets_resolved: 1" in manifest
    assert "Quartet Preview" in manifest
    assert summary_json["command"] == "ipyrad2 baba -d snps.hdf5 -o OUT --tests quartets.tsv -i imap.tsv -g minmap.tsv"
    assert summary_json["tool"] == "baba"
    assert summary_json["results_summary"]["quartets_rows"] == 1
    assert summary_json["results_summary"]["rooted_orientation_counts"] == {
        "bbaa": 1,
        "dmin": 1,
        "input": 1,
    }
    assert summary_json["tables"]["quartets"][0]["p1"] == "A"
    assert summary_json["tables"]["quartets"][0]["estimate_se"] is None
    assert summary_json["outputs"]["summary_json"].endswith("baba.summary.json")


def test_run_baba_method_manual_sample_quartet_without_imap_uses_sample_namespace(
    tmp_path: Path,
) -> None:
    h5 = _write_baba_h5(tmp_path / "snps.hdf5")
    tests = tmp_path / "quartets.tsv"
    tests.write_text("a1 b1 c1 o1\n", encoding="utf-8")

    run_baba_method(
        data=h5,
        name="sample",
        outdir=tmp_path / "OUT",
        tests=tests,
        tree=None,
        imap=None,
        minmap=None,
        min_sample_coverage=1,
        exclude=None,
        include_reference=False,
        resampling="none",
        bootstrap_replicates=10,
        jackknife_block_bp=1,
        jackknife_block_loci=2,
        seed=3,
        f_branch=False,
        f_branch_p_threshold=0.01,
        write_block_table=False,
        clustering_stats=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    quartets = pd.read_csv(tmp_path / "OUT" / "sample.quartets.tsv", sep="\t")
    rooted = pd.read_csv(tmp_path / "OUT" / "sample.rooted.tsv", sep="\t")
    summary_json = _read_json(tmp_path / "OUT" / "sample.summary.json")
    assert quartets.loc[0, "p1"] == "a1"
    assert quartets.loc[0, "n_sites_tested"] == 28
    assert quartets.loc[0, "abba"] == pytest.approx(18.0)
    assert quartets.loc[0, "baba"] == pytest.approx(5.0)
    assert quartets.loc[0, "bbaa"] == pytest.approx(5.0)
    assert quartets.loc[0, "d_stat"] == pytest.approx((18.0 - 5.0) / (18.0 + 5.0))
    assert set(rooted["orientation"]) == {"input", "bbaa", "dmin"}
    assert "command" not in summary_json


def test_run_baba_method_rooted_rows_match_dsuite_style_orientations(
    tmp_path: Path,
) -> None:
    h5 = _write_rooted_order_h5(tmp_path / "rooted.hdf5")
    tests = tmp_path / "quartets.tsv"
    tests.write_text("p1 p2 p3 o1\n", encoding="utf-8")

    run_baba_method(
        data=h5,
        name="rooted",
        outdir=tmp_path / "OUT",
        tests=tests,
        tree=None,
        imap=None,
        minmap=None,
        min_sample_coverage=1,
        exclude=None,
        include_reference=False,
        resampling="none",
        bootstrap_replicates=10,
        jackknife_block_bp=1,
        jackknife_block_loci=2,
        seed=5,
        f_branch=False,
        f_branch_p_threshold=0.01,
        write_block_table=False,
        clustering_stats=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    rooted = pd.read_csv(tmp_path / "OUT" / "rooted.rooted.tsv", sep="\t")

    input_row = rooted[rooted["orientation"] == "input"].iloc[0]
    bbaa_row = rooted[rooted["orientation"] == "bbaa"].iloc[0]
    dmin_row = rooted[rooted["orientation"] == "dmin"].iloc[0]

    assert input_row["p1"] == "p1"
    assert input_row["p2"] == "p2"
    assert input_row["p3"] == "p3"
    assert input_row["d_stat"] == pytest.approx((10.0 - 2.0) / (10.0 + 2.0))
    assert input_row["f_g"] == pytest.approx((10.0 - 2.0) / 10.0)

    assert bbaa_row["p1"] == "p3"
    assert bbaa_row["p2"] == "p2"
    assert bbaa_row["p3"] == "p1"
    assert bbaa_row["d_stat"] == pytest.approx((6.0 - 2.0) / (6.0 + 2.0))
    assert bbaa_row["f_g"] == pytest.approx((6.0 - 2.0) / 6.0)

    assert dmin_row["p1"] == "p1"
    assert dmin_row["p2"] == "p3"
    assert dmin_row["p3"] == "p2"
    assert dmin_row["d_stat"] == pytest.approx((10.0 - 6.0) / (10.0 + 6.0))
    assert dmin_row["f_g"] == pytest.approx((10.0 - 6.0) / 10.0)


def test_run_baba_method_tree_mode_writes_tree_branch_and_block_outputs(
    tmp_path: Path,
) -> None:
    h5 = _write_baba_h5(tmp_path / "snps.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")
    minmap = _write_minmap(tmp_path / "minmap.tsv")
    tree = tmp_path / "species.nwk"
    tree.write_text("((((A,B),C),D),O);\n", encoding="utf-8")

    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}")
    try:
        run_baba_method(
            data=h5,
            name="tree",
            outdir=tmp_path / "OUT",
            tests=None,
            tree=tree,
            imap=imap,
            minmap=minmap,
            min_sample_coverage=1,
            exclude=None,
            include_reference=False,
            resampling="auto",
            bootstrap_replicates=10,
            jackknife_block_bp=1,
            jackknife_block_loci=2,
            seed=4,
            f_branch=True,
            f_branch_p_threshold=0.01,
            write_block_table=True,
            clustering_stats=True,
            cores=1,
            force=True,
            log_level="INFO",
            logged_command="ipyrad2 baba -d snps.hdf5 -o OUT --tree species.nwk --f-branch",
        )
    finally:
        logger.remove(sink_id)

    quartets = pd.read_csv(tmp_path / "OUT" / "tree.quartets.tsv", sep="\t")
    rooted = pd.read_csv(tmp_path / "OUT" / "tree.rooted.tsv", sep="\t")
    resolved = pd.read_csv(tmp_path / "OUT" / "tree.tests.resolved.tsv", sep="\t")
    blocks = pd.read_csv(tmp_path / "OUT" / "tree.blocks.tsv", sep="\t")
    f_branch = pd.read_csv(tmp_path / "OUT" / "tree.f_branch.tsv", sep="\t")
    f_branch_matrix = pd.read_csv(tmp_path / "OUT" / "tree.f_branch.matrix.tsv", sep="\t")
    f_branch_z = pd.read_csv(tmp_path / "OUT" / "tree.f_branch.z.tsv", sep="\t")
    f_branch_p = pd.read_csv(tmp_path / "OUT" / "tree.f_branch.p.tsv", sep="\t")
    tree_used = (tmp_path / "OUT" / "tree.tree.used.nwk").read_text(encoding="utf-8")
    manifest = (tmp_path / "OUT" / "tree.manifest.txt").read_text(encoding="utf-8")
    summary_json = _read_json(tmp_path / "OUT" / "tree.summary.json")

    assert len(resolved.index) == 5
    assert len(quartets.index) == 5
    assert len(rooted.index) == 20
    assert set(["n_blocks", "mean_block_d", "sign_switches"]).issubset(quartets.columns)
    focal = quartets[
        (quartets["p1"] == "A")
        & (quartets["p2"] == "B")
        & (quartets["p3"] == "C")
        & (quartets["p4"] == "O")
    ].iloc[0]
    focal_tree = rooted[
        (rooted["orientation"] == "tree")
        & (rooted["p1"] == "A")
        & (rooted["p2"] == "B")
        & (rooted["p3"] == "C")
        & (rooted["p4"] == "O")
    ].iloc[0]
    assert focal["resampling_mode"] == "jackknife"
    assert focal["resampling_unit"] == "physical_block"
    assert focal["resampling_units"] == 24
    assert np.isfinite(focal["estimate_se"])
    assert np.isfinite(focal_tree["f_g"])
    assert len(blocks.index) >= 24
    assert not f_branch.empty
    assert {"f_branch_raw", "f_branch", "z_branch", "p_branch"}.issubset(f_branch.columns)
    assert {"branch", "branch_descendants", "A", "B", "C", "D", "O"}.issubset(f_branch_matrix.columns)
    assert {"branch", "branch_descendants"}.issubset(f_branch_z.columns)
    assert {"branch", "branch_descendants"}.issubset(f_branch_p.columns)
    assert "O" in tree_used
    assert manifest.startswith("CMD: ipyrad2 baba -d snps.hdf5 -o OUT --tree species.nwk --f-branch")
    assert "F-branch Preview" in manifest
    assert "Block Preview" in manifest
    assert "tree_skipped_balanced_quartets: 0" in manifest
    assert summary_json["command"] == "ipyrad2 baba -d snps.hdf5 -o OUT --tree species.nwk --f-branch"
    assert summary_json["results_summary"]["quartets_rows"] == 5
    assert summary_json["results_summary"]["rooted_orientation_counts"]["tree"] == 5
    assert len(summary_json["tables"]["f_branch"]) == len(f_branch.index)
    assert len(summary_json["tables"]["blocks"]) == len(blocks.index)
    assert summary_json["outputs"]["manifest"].endswith("tree.manifest.txt")
    joined = "\n".join(str(msg) for msg in messages)
    assert "resolving baba tree in population namespace" in joined
    assert "retained 7 samples, 28 SNPs, 28 loci, and 2 scaffolds after filtering" in joined
    assert "processing 5 quartets" in joined
    assert "processed quartet 1/5 A B C D" in joined
    assert "processed quartet 2/5 A B C O" in joined
    assert "aggregating tree-oriented f_branch summaries" in joined
    assert "writing baba outputs to" in joined


def test_run_baba_method_rejects_vcf_and_invalid_minmap_usage(tmp_path: Path) -> None:
    with pytest.raises(IPyradError, match="requires an SNP-capable HDF5 input"):
        run_baba_method(
            data=tmp_path / "variants.vcf.gz",
            name="bad",
            outdir=tmp_path / "OUT",
            tests=tmp_path / "quartets.tsv",
            tree=None,
            imap=None,
            minmap=None,
            min_sample_coverage=1,
            exclude=None,
            include_reference=False,
            resampling="none",
            bootstrap_replicates=10,
            jackknife_block_bp=1,
            jackknife_block_loci=2,
            seed=1,
            f_branch=False,
            f_branch_p_threshold=0.01,
            write_block_table=False,
            clustering_stats=False,
            cores=1,
            force=True,
            log_level="INFO",
        )

    h5 = _write_baba_h5(tmp_path / "snps.hdf5")
    tests = tmp_path / "quartets.tsv"
    tests.write_text("a1 b1 c1 o1\n", encoding="utf-8")
    minmap = _write_minmap(tmp_path / "minmap.tsv")
    with pytest.raises(IPyradError, match="--minmap requires --imap."):
        run_baba_method(
            data=h5,
            name="bad2",
            outdir=tmp_path / "OUT2",
            tests=tests,
            tree=None,
            imap=None,
            minmap=minmap,
            min_sample_coverage=1,
            exclude=None,
            include_reference=False,
            resampling="none",
            bootstrap_replicates=10,
            jackknife_block_bp=1,
            jackknife_block_loci=2,
            seed=1,
            f_branch=False,
            f_branch_p_threshold=0.01,
            write_block_table=False,
            clustering_stats=False,
            cores=1,
            force=True,
            log_level="INFO",
        )
