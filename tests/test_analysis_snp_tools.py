import gzip
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from loguru import logger

from ipyrad2.analysis.converters.vcf_to_hdf5 import run_vcf_to_hdf5
from ipyrad2.analysis.extracters.snps_extracter import (
    SNPSMAP_COLUMNS,
    SNPsExtracter,
    run_snps_extracter,
)
from ipyrad2.utils.exceptions import IPyradError


def _write_snps_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    genos = np.array(
        [
            [[0, 0, ord("A")], [0, 1, ord("R")], [1, 1, ord("G")]],
            [[0, 0, ord("A")], [0, 0, ord("A")], [0, 0, ord("A")]],
            [[1, 1, ord("G")], [0, 0, ord("A")], [0, 1, ord("R")]],
        ],
        dtype=np.uint8,
    )
    snpsmap = np.array(
        [
            [0, 0, 0, 0, 100],
            [0, 1, 1, 0, 120],
            [1, 0, 0, 1, 5],
        ],
        dtype=np.uint32,
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(
            ["s1", "assembly_reference_sequence", "s3"],
            dtype=string_dtype,
        )
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A"), ord("A"), ord("A")], dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(SNPSMAP_COLUMNS, dtype=string_dtype)
    return path


def _write_snps_h5_without_reference(path: Path) -> Path:
    path = _write_snps_h5(path)
    with h5py.File(path, "a") as io5:
        del io5["reference"]
    return path


def _write_missingness_snps_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    genos = np.array(
        [
            [[0, 0, ord("A")], [0, 1, ord("R")]],
            [[255, 255, ord("N")], [255, 255, ord("N")]],
            [[1, 1, ord("G")], [1, 1, ord("G")]],
        ],
        dtype=np.uint8,
    )
    snpsmap = np.array(
        [
            [0, 0, 0, 0, 10],
            [1, 0, 0, 0, 20],
        ],
        dtype=np.uint32,
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(["s1", "s2", "s3"], dtype=string_dtype)
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A"), ord("A")], dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(SNPSMAP_COLUMNS, dtype=string_dtype)
    return path


def _write_mincov_minmap_snps_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["a1", "a2", "a3", "a4", "b1", "b2"]
    doses = np.array(
        [
            [0, 0, 255],
            [1, 0, 0],
            [255, 1, 0],
            [255, 2, 1],
            [2, 0, 1],
            [2, 255, 2],
        ],
        dtype=np.uint8,
    )
    genos = np.zeros((doses.shape[0], doses.shape[1], 3), dtype=np.uint8)
    for row in range(doses.shape[0]):
        for col in range(doses.shape[1]):
            dose = int(doses[row, col])
            if dose == 255:
                genos[row, col, :2] = 255
                genos[row, col, 2] = ord("N")
            elif dose == 0:
                genos[row, col, :2] = [0, 0]
                genos[row, col, 2] = ord("A")
            elif dose == 1:
                genos[row, col, :2] = [0, 1]
                genos[row, col, 2] = ord("R")
            else:
                genos[row, col, :2] = [1, 1]
                genos[row, col, 2] = ord("G")

    snpsmap = np.array(
        [
            [0, 0, 0, 0, 10],
            [1, 0, 0, 0, 20],
            [2, 0, 0, 0, 30],
        ],
        dtype=np.uint32,
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A"), ord("A"), ord("A")], dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(SNPSMAP_COLUMNS, dtype=string_dtype)
    return path


def _write_assembly_style_snps_h5(path: Path) -> Path:
    """Write an assemble-style HDF5 where top-level names include the reference."""
    string_dtype = h5py.string_dtype(encoding="utf-8")
    genos = np.array(
        [
            [[0, 0, ord("A")], [0, 1, ord("R")], [1, 1, ord("G")]],
            [[1, 1, ord("G")], [0, 0, ord("A")], [0, 1, ord("R")]],
        ],
        dtype=np.uint8,
    )
    snpsmap = np.array(
        [
            [0, 0, 0, 0, 100],
            [0, 1, 1, 0, 120],
            [1, 0, 0, 1, 5],
        ],
        dtype=np.uint32,
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(
            ["assembly_reference_sequence", "s1", "s3"],
            dtype=string_dtype,
        )
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A"), ord("A"), ord("A")], dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(SNPSMAP_COLUMNS, dtype=string_dtype)
    return path


def _write_depth_qual_snps_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    genos = np.array(
        [
            [[0, 0, ord("A")], [0, 0, ord("A")]],
            [[0, 1, ord("R")], [0, 1, ord("R")]],
        ],
        dtype=np.uint8,
    )
    snpsmap = np.array(
        [
            [0, 0, 0, 0, 10],
            [1, 0, 0, 0, 20],
        ],
        dtype=np.uint32,
    )
    sample_dp = np.array(
        [
            [5, 1],
            [5, 5],
        ],
        dtype=np.uint32,
    )
    site_qual = np.array([60.0, 10.0], dtype=np.float32)
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(["s1", "s2"], dtype=string_dtype)
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A"), ord("A")], dtype=np.uint8),
        )
        io5.create_dataset("sample_dp", data=sample_dp)
        io5.create_dataset("site_qual", data=site_qual)
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(SNPSMAP_COLUMNS, dtype=string_dtype)
    return path


def _write_imap_files(tmp_path: Path, *, include_reference: bool) -> tuple[Path, Path]:
    imap_path = tmp_path / "imap.tsv"
    if include_reference:
        imap_path.write_text(
            "s1\tpop1\nassembly_reference_sequence\tpop1\ns3\tpop2\n",
            encoding="utf-8",
        )
    else:
        imap_path.write_text(
            "s1\tpop1\ns3\tpop2\n",
            encoding="utf-8",
        )
    minmap_path = tmp_path / "minmap.tsv"
    minmap_path.write_text("pop1\t1\npop2\t1\n", encoding="utf-8")
    return imap_path, minmap_path


def _write_vcf(path: Path) -> Path:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=testvcf",
                "##reference=testref.fa",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
                "chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT:DP\t0/0:8\t0/1:7",
                "chr1\t5\t.\tC\tT\t15\tPASS\t.\tDP:GT\t6:1/1\t9:0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def test_run_snps_extracter_defaults_to_unlinked_output_with_seed_and_excludes_reference_by_default(
    tmp_path: Path,
) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")
    expected = tool.get_view(subsample=True, random_seed=7, log_level="INFO")

    run_snps_extracter(
        data=h5,
        name="snpset",
        outdir=tmp_path / "OUT",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        cores=1,
        force=True,
        log_level="INFO",
        random_seed=7,
    )

    outdir = tmp_path / "OUT"
    samples = (outdir / "snpset.samples.txt").read_text(encoding="utf-8").splitlines()
    sample_summary = pd.read_csv(outdir / "snpset.sample_data_summary.tsv", sep="\t")
    genos = np.load(outdir / "snpset.genos.npy")
    snps = np.load(outdir / "snpset.snps.npy")
    snpsmap = pd.read_csv(outdir / "snpset.snpsmap.tsv", sep="\t")
    stats_text = (outdir / "snpset.stats.txt").read_text(encoding="utf-8")

    assert samples == ["s1", "s3"]
    assert np.array_equal(genos, expected.genos)
    assert np.array_equal(snps, expected.snps)
    assert list(snpsmap.columns) == SNPSMAP_COLUMNS
    assert snpsmap.shape == (2, 5)
    assert snpsmap["loc"].nunique() == 2
    assert list(sample_summary.columns) == [
        "sample",
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputation_algorithm",
        "imputed_genotype_fraction",
    ]
    assert sample_summary["imputation_algorithm"].eq("not-imputed").all()
    assert np.allclose(
        sample_summary["missing_fraction"],
        sample_summary["post_imputation_missing_fraction"],
    )
    assert np.allclose(sample_summary["imputed_genotype_fraction"], 0.0)
    assert "include_reference: False" in stats_text
    assert "subsample: True" in stats_text
    assert "random_seed: 7" in stats_text
    assert "impute_method: none" in stats_text
    assert "imputation_algorithm: not-imputed" in stats_text
    assert "written_formats: genos, snps, snpsmap, samples, sample_data_summary" in stats_text
    assert "linked_post_filter_snps: 3" in stats_text
    assert "exported_snps: 2" in stats_text


def test_run_snps_extracter_no_subsample_writes_linked_outputs(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")

    run_snps_extracter(
        data=h5,
        name="snpset",
        outdir=tmp_path / "OUT",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        cores=1,
        force=True,
        log_level="INFO",
        subsample=False,
    )

    genos = np.load(tmp_path / "OUT" / "snpset.genos.npy")
    snpsmap = pd.read_csv(tmp_path / "OUT" / "snpset.snpsmap.tsv", sep="\t")
    stats_text = (tmp_path / "OUT" / "snpset.stats.txt").read_text(encoding="utf-8")

    assert genos.shape == (2, 3)
    assert snpsmap.shape == (3, 5)
    assert "subsample: False" in stats_text


def test_run_snps_extracter_global_imputation_updates_all_written_outputs(tmp_path: Path) -> None:
    h5 = _write_missingness_snps_h5(tmp_path / "snps.hdf5")

    run_snps_extracter(
        data=h5,
        name="snpset",
        outdir=tmp_path / "OUT",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        cores=1,
        force=True,
        log_level="INFO",
        subsample=False,
        random_seed=3,
        impute_method="zero",
        write_phylip=True,
        write_nexus=True,
        write_fasta=True,
    )

    outdir = tmp_path / "OUT"
    genos = np.load(outdir / "snpset.genos.npy")
    snps = np.load(outdir / "snpset.snps.npy")
    sample_summary = pd.read_csv(outdir / "snpset.sample_data_summary.tsv", sep="\t")
    stats_text = (outdir / "snpset.stats.txt").read_text(encoding="utf-8")
    phy = (outdir / "snpset.phy").read_text(encoding="utf-8")
    nex = (outdir / "snpset.nex").read_text(encoding="utf-8")
    fa = (outdir / "snpset.fa").read_text(encoding="utf-8")

    assert not np.any(genos == 255)
    assert not np.any(snps == ord("N"))
    assert sample_summary["imputation_algorithm"].eq("zero-fill").all()
    assert np.allclose(sample_summary["post_imputation_missing_fraction"], 0.0)
    assert sample_summary.loc[sample_summary["sample"] == "s2", "imputed_genotype_fraction"].item() == 1.0
    assert "impute_method: zero-fill" in stats_text
    assert "imputation_algorithm: zero-fill" in stats_text
    assert "written_formats: genos, snps, snpsmap, samples, sample_data_summary, phylip, nexus, fasta" in stats_text
    assert "3 2" in phy
    assert "datatype=dna" in nex
    assert ">s1" in fa


def test_run_snps_extracter_writes_treemix_and_eems_outputs(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    imap_path, minmap_path = _write_imap_files(tmp_path, include_reference=False)

    run_snps_extracter(
        data=h5,
        name="snpset",
        outdir=tmp_path / "OUT",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap_path,
        minmap=minmap_path,
        exclude=None,
        include_reference=False,
        cores=1,
        force=True,
        log_level="INFO",
        subsample=False,
        write_treemix=True,
        write_eems=True,
    )

    with gzip.open(tmp_path / "OUT" / "snpset.treemix.gz", "rt", encoding="utf-8") as infile:
        treemix_lines = infile.read().splitlines()
    eems = pd.read_csv(tmp_path / "OUT" / "snpset.eems", sep="\t", header=None).to_numpy()

    assert treemix_lines[0] == "pop1 pop2"
    assert len(treemix_lines) == 4
    assert eems.shape == (2, 2)
    assert np.allclose(np.diag(eems), 0.0)
    assert np.allclose(eems, eems.T)


def test_snps_extracter_logs_default_minmap_message_without_claiming_override(
    tmp_path: Path,
) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap={"pop1": ["s1"], "pop2": ["s3"]},
        minmap=None,
        include_reference=False,
        cores=1,
    )
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="INFO")
    try:
        tool._get_imap_minmap(tool.imap, {})
    finally:
        logger.remove(sink_id)

    assert tool.minmap == {"pop1": 0, "pop2": 0}
    assert any(
        "global `-m` filter still applies" in msg
        and "defaulting per-population minimums to 0" in msg
        and "`-g` has no effect" in msg
        for msg in messages
    )


def test_run_snps_extracter_logs_extraction_then_no_imputation_then_export_summary(
    tmp_path: Path,
) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="INFO")
    try:
        run_snps_extracter(
            data=h5,
            name="snpset",
            outdir=tmp_path / "OUT",
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            cores=1,
            force=True,
            log_level="INFO",
            subsample=True,
            random_seed=7,
            impute_method=None,
        )
    finally:
        logger.remove(sink_id)

    extraction_idx = next(
        idx for idx, msg in enumerate(messages) if "SNP extraction summary" in msg
    )
    imputation_idx = next(
        idx
        for idx, msg in enumerate(messages)
        if "snpex SNP imputation: no imputation performed" in msg
    )
    export_idx = next(
        idx
        for idx, msg in enumerate(messages)
        if "snpex exported SNP summary:" in msg
    )

    assert extraction_idx < imputation_idx < export_idx
    assert any(
        "snpex SNP imputation: no imputation performed" in msg
        and "prepared_matrix_scope=subsampled_unlinked" in msg
        for msg in messages
    )
    assert any(
        "snpex exported SNP summary:" in msg
        and "prepared_matrix_scope=subsampled_unlinked" in msg
        and "linked_post_filter_snps=3" in msg
        and "exported_snps=2" in msg
        for msg in messages
    )
    assert not any("subsampled " in msg for msg in messages)


def test_snps_extracter_applies_global_mincov_and_population_minmap_filters(
    tmp_path: Path,
) -> None:
    h5 = _write_mincov_minmap_snps_h5(tmp_path / "snps.hdf5")

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=5,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap={"popA": ["a1", "a2", "a3", "a4"], "popB": ["b1", "b2"]},
        minmap={"popA": 2, "popB": 2},
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert int(tool.stats["filter_by_mincov"]) == 1
    assert int(tool.stats["filter_by_minmap"]) == 1
    assert int(tool.stats["post_filter_snps"]) == 1
    assert tool.snpsmap.shape[0] == 1


def test_snps_extracter_defaults_missing_minmap_to_zero_effect_per_population(
    tmp_path: Path,
) -> None:
    h5 = _write_mincov_minmap_snps_h5(tmp_path / "snps.hdf5")

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=5,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap={"popA": ["a1", "a2", "a3", "a4"], "popB": ["b1", "b2"]},
        minmap=None,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert tool.minmap == {"popA": 0, "popB": 0}
    assert int(tool.stats["filter_by_mincov"]) == 1
    assert int(tool.stats["filter_by_minmap"]) == 0
    assert int(tool.stats["post_filter_snps"]) == 2
    assert tool.snpsmap.shape[0] == 2


def test_snps_extracter_masks_low_depth_genotypes_before_site_filters(tmp_path: Path) -> None:
    h5 = _write_depth_qual_snps_h5(tmp_path / "snps.hdf5")

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        min_genotype_depth=2,
        min_site_qual=0.0,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert tool.genos.shape == (2, 1)
    assert int(tool.stats["masked_genotypes_by_min_depth"]) == 1
    assert int(tool.stats["filter_by_mincov"]) == 1
    assert int(tool.stats["filter_by_min_site_qual"]) == 0


def test_snps_extracter_filters_low_site_qual_sites(tmp_path: Path) -> None:
    h5 = _write_depth_qual_snps_h5(tmp_path / "snps.hdf5")

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        min_genotype_depth=0,
        min_site_qual=20.0,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert tool.genos.shape == (2, 1)
    assert int(tool.stats["filter_by_min_site_qual"]) == 1
    assert int(tool.stats["masked_genotypes_by_min_depth"]) == 0


def test_snps_extracter_rejects_new_filters_on_legacy_hdf5(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")

    with pytest.raises(IPyradError, match="sample_dp"):
        SNPsExtracter(
            data=h5,
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            min_genotype_depth=2,
            min_site_qual=0.0,
            include_reference=False,
            cores=1,
        )

    with pytest.raises(IPyradError, match="site_qual"):
        SNPsExtracter(
            data=h5,
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            min_genotype_depth=0,
            min_site_qual=20.0,
            include_reference=False,
            cores=1,
        )


def test_global_snpex_imputation_requires_reference_dataset(tmp_path: Path) -> None:
    h5 = _write_snps_h5_without_reference(tmp_path / "snps.hdf5")

    with pytest.raises(IPyradError, match="Global SNP imputation for snpex requires the HDF5 `reference` dataset"):
        run_snps_extracter(
            data=h5,
            name="snpset",
            outdir=tmp_path / "OUT",
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            cores=1,
            force=True,
            log_level="INFO",
            impute_method="sample",
        )


def test_snps_extracter_imap_file_includes_reference_without_flag(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    imap_path, minmap_path = _write_imap_files(tmp_path, include_reference=True)

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap_path,
        minmap=minmap_path,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert tool.snames == ["s1", "assembly_reference_sequence", "s3"]
    assert tool.imap == {
        "pop1": ["s1", "assembly_reference_sequence"],
        "pop2": ["s3"],
    }


def test_snps_extracter_expands_glob_imap_entries_from_file(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    imap_path = tmp_path / "imap.tsv"
    imap_path.write_text(
        "s*\tpop1\n"
        "assembly_reference_sequence\tpop2\n",
        encoding="utf-8",
    )
    minmap_path = tmp_path / "minmap.tsv"
    minmap_path.write_text("pop1\t1\npop2\t1\n", encoding="utf-8")

    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap_path,
        minmap=minmap_path,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")

    assert tool.snames == ["s1", "assembly_reference_sequence", "s3"]
    assert tool.imap == {
        "pop1": ["s1", "s3"],
        "pop2": ["assembly_reference_sequence"],
    }
    assert tool.minmap == {"pop1": 1, "pop2": 1}


def test_snps_extracter_include_reference_with_imap_requires_reference_assignment(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    imap_path, minmap_path = _write_imap_files(tmp_path, include_reference=False)

    with pytest.raises(
        IPyradError,
        match="assembly_reference_sequence was requested with -R, but it must also be assigned to an IMAP group.",
    ):
        SNPsExtracter(
            data=h5,
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=imap_path,
            minmap=minmap_path,
            include_reference=True,
            cores=1,
        )


def test_snps_extracter_drops_high_missing_samples_and_reruns_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_missingness_snps_h5(tmp_path / "snps.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=2,
        max_sample_missing=0.5,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=False,
        cores=1,
    )
    calls = 0
    original = tool._run_filter_pass

    def _wrapped(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(tool, "_run_filter_pass", _wrapped)

    tool.run(log_level="INFO")

    assert calls == 2
    assert tool.dropped_samples_by_missing == ["s2"]
    assert tool.snames == ["s1", "s3"]
    assert tool.imap == {"all": ["s1", "s3"]}


def test_snps_extracter_handles_assembly_style_names_without_reference_genotypes(
    tmp_path: Path,
) -> None:
    h5 = _write_assembly_style_snps_h5(tmp_path / "assembly_style.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=False,
        cores=1,
    )

    tool.run(log_level="INFO")

    assert tool.dbnames == ["assembly_reference_sequence", "s1", "s3"]
    assert tool.snames == ["s1", "s3"]
    assert tool.genos.shape == (2, 3)


def test_snps_extracter_include_reference_synthesizes_reference_row_for_assembly_style_hdf5(
    tmp_path: Path,
) -> None:
    h5 = _write_assembly_style_snps_h5(tmp_path / "assembly_style.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=True,
        cores=1,
    )

    tool.run(log_level="INFO")

    assert tool.snames == ["assembly_reference_sequence", "s1", "s3"]
    assert tool.genos.shape == (3, 3)
    assert np.array_equal(tool.genos[0], np.zeros(3, dtype=np.uint8))
    assert tool.sample_missing["assembly_reference_sequence"] == 0.0


def test_run_snps_extracter_writes_plink_files_for_exported_snp_view(tmp_path: Path) -> None:
    h5 = _write_snps_h5(tmp_path / "snps.hdf5")
    tool = SNPsExtracter(
        data=h5,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")
    expected = tool.get_view(subsample=True, random_seed=11, log_level="INFO")

    run_snps_extracter(
        data=h5,
        name="snpset",
        outdir=tmp_path / "OUT",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        cores=1,
        force=True,
        log_level="INFO",
        random_seed=11,
        write_plink=True,
    )

    fam_lines = (tmp_path / "OUT" / "snpset.fam").read_text(encoding="utf-8").splitlines()
    bim_lines = (tmp_path / "OUT" / "snpset.bim").read_text(encoding="utf-8").splitlines()
    bed_bytes = (tmp_path / "OUT" / "snpset.bed").read_bytes()

    assert fam_lines == [
        "s1\ts1\t0\t0\t0\t-9",
        "s3\ts3\t0\t0\t0\t-9",
    ]
    assert len(bim_lines) == expected.snpsmap.shape[0]
    assert bed_bytes[:3] == bytes([0x6C, 0x1B, 0x01])


def test_plink_export_requires_reference_dataset(tmp_path: Path) -> None:
    h5 = _write_snps_h5_without_reference(tmp_path / "snps.hdf5")

    with pytest.raises(IPyradError, match="PLINK export requires the HDF5 `reference` dataset"):
        run_snps_extracter(
            data=h5,
            name="snpset",
            outdir=tmp_path / "OUT",
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            cores=1,
            force=True,
            log_level="INFO",
            write_plink=True,
        )


def test_run_vcf_to_hdf5_writes_snp_compatible_database(tmp_path: Path) -> None:
    vcf = _write_vcf(tmp_path / "variants.vcf")

    outpath = run_vcf_to_hdf5(
        data=vcf,
        name="snps",
        outdir=tmp_path / "OUT",
        ld_block_size=10,
        force=True,
    )

    assert outpath == tmp_path / "OUT" / "snps.hdf5"
    assert outpath.exists()

    with h5py.File(outpath, "r") as io5:
        assert int(io5.attrs["nsnps"]) == 2
        assert list(io5.attrs["names"]) == ["s1", "s2"]
        assert io5["genos"].shape == (2, 2, 3)
        assert io5["snpsmap"].shape == (2, 5)
        assert np.array_equal(
            io5["reference"][:],
            np.array([ord("A"), ord("C")], dtype=np.uint8),
        )
        np.testing.assert_array_equal(
            io5["sample_dp"][:],
            np.array([[8, 6], [7, 9]], dtype=np.uint32),
        )
        np.testing.assert_allclose(
            io5["site_qual"][:],
            np.array([60.0, 15.0], dtype=np.float32),
        )

    tool = SNPsExtracter(
        data=outpath,
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        include_reference=False,
        cores=1,
    )
    tool.run(log_level="INFO")
    assert tool.genos.shape == (2, 2)
    assert tool.snps.shape == (2, 2)


def test_run_vcf_to_hdf5_requires_force_to_overwrite(tmp_path: Path) -> None:
    vcf = _write_vcf(tmp_path / "variants.vcf")

    run_vcf_to_hdf5(
        data=vcf,
        name="snps",
        outdir=tmp_path / "OUT",
        ld_block_size=10,
        force=True,
    )

    with pytest.raises(IPyradError, match="HDF5 file already exists"):
        run_vcf_to_hdf5(
            data=vcf,
            name="snps",
            outdir=tmp_path / "OUT",
            ld_block_size=10,
            force=False,
        )
