from types import SimpleNamespace

import numpy as np
import pytest

import ipyrad2.analysis as analysis_mod
import ipyrad2.analysis.methods.common as common_mod
from ipyrad2.analysis.extracters.snps_extracter import _MISSING_GENO
from ipyrad2.analysis.methods.snps_imputer import SNPsImputer
from ipyrad2.utils.exceptions import IPyradError


def test_snps_imputer_public_export_is_available() -> None:
    assert analysis_mod.SNPsImputer is SNPsImputer


def test_snps_imputer_sample_mode_respects_population_groups() -> None:
    data = np.array(
        [
            [2, _MISSING_GENO],
            [2, 0],
            [0, _MISSING_GENO],
            [0, 2],
        ],
        dtype=np.uint8,
    )
    imap = {"pop1": ["a", "b"], "pop2": ["c", "d"]}

    result = SNPsImputer(data, ["a", "b", "c", "d"], imap=imap, impute_method="sample", quiet=True).run()

    expected = np.array(
        [
            [2, 0],
            [2, 0],
            [0, 2],
            [0, 2],
        ],
        dtype=np.uint8,
    )
    np.testing.assert_array_equal(result, expected)
    assert result.dtype == data.dtype
    assert result.shape == data.shape
    np.testing.assert_array_equal(
        data,
        np.array(
            [
                [2, _MISSING_GENO],
                [2, 0],
                [0, _MISSING_GENO],
                [0, 2],
            ],
            dtype=np.uint8,
        ),
    )


def test_snps_imputer_sample_mode_defaults_to_global_group_when_imap_is_none() -> None:
    data = np.array(
        [
            [2, _MISSING_GENO],
            [2, 0],
        ],
        dtype=np.uint8,
    )

    result = SNPsImputer(data, ["a", "b"], imap=None, impute_method="sample", quiet=True).run()

    np.testing.assert_array_equal(result, np.array([[2, 0], [2, 0]], dtype=np.uint8))


@pytest.mark.parametrize("method", [None, False, "none", "zero", "zero-fill"])
def test_snps_imputer_null_modes_fill_missing_with_zero(method) -> None:
    data = np.array(
        [
            [0, _MISSING_GENO],
            [_MISSING_GENO, 2],
        ],
        dtype=np.uint8,
    )

    result = SNPsImputer(data, ["a", "b"], impute_method=method, quiet=True).run()

    np.testing.assert_array_equal(result, np.array([[0, 0], [0, 2]], dtype=np.uint8))


def test_snps_imputer_rejects_unsupported_methods() -> None:
    with pytest.raises(IPyradError, match="Unsupported SNPsImputer impute_method"):
        SNPsImputer(np.zeros((2, 2), dtype=np.uint8), ["a", "b"], impute_method="random")


@pytest.mark.parametrize(
    "imap, match",
    [
        ({"pop1": ["missing"]}, "not present in the genotype matrix"),
        ({"pop1": ["a"], "pop2": ["a"]}, "multiple groups"),
        ({"pop1": ["a", "a"]}, "contains duplicate sample names"),
        ({"pop1": []}, "is empty"),
    ],
)
def test_snps_imputer_validates_imap(imap, match) -> None:
    data = np.zeros((2, 2), dtype=np.uint8)

    with pytest.raises(IPyradError, match=match):
        SNPsImputer(data, ["a", "b"], imap=imap)


def test_snps_imputer_handles_all_missing_sites_without_nan_probabilities() -> None:
    data = np.array(
        [
            [_MISSING_GENO, _MISSING_GENO],
            [_MISSING_GENO, 2],
        ],
        dtype=np.uint8,
    )

    result = SNPsImputer(data, ["a", "b"], imap={"pop1": ["a", "b"]}, impute_method="sample", quiet=True).run()

    np.testing.assert_array_equal(result, np.array([[0, 2], [0, 2]], dtype=np.uint8))


def test_snps_imputer_no_missing_data_returns_unchanged_copy() -> None:
    data = np.array([[0, 1], [2, 0]], dtype=np.uint8)

    result = SNPsImputer(data, ["a", "b"], impute_method="sample", quiet=True).run()

    np.testing.assert_array_equal(result, data)
    assert result.dtype == data.dtype
    assert result is not data


def test_snps_imputer_validates_data_shape_and_name_count() -> None:
    with pytest.raises(IPyradError, match="2D genotype matrix"):
        SNPsImputer(np.array([0, 1, 2], dtype=np.uint8), ["a"])
    with pytest.raises(IPyradError, match="names length must match"):
        SNPsImputer(np.zeros((2, 2), dtype=np.uint8), ["a"])


def test_shared_numerical_input_uses_snps_imputer_and_reports_imputation(monkeypatch) -> None:
    observed = {}

    class DummyImputer:
        def __init__(self, data, names, imap, impute_method, quiet):
            observed["data"] = data.copy()
            observed["names"] = list(names)
            observed["imap"] = imap
            observed["impute_method"] = impute_method
            observed["quiet"] = quiet

        def run(self):
            return np.array([[0, 1], [2, 2]], dtype=np.uint8)

    class DummyExtracter:
        snames = ["a", "b"]
        imap = {"pop1": ["a", "b"]}

        def get_view(self, *, subsample, random_seed, log_level):
            assert subsample is True
            assert random_seed == 13
            assert log_level == "DEBUG"
            return SimpleNamespace(
                genos=np.array([[0, _MISSING_GENO], [2, 2]], dtype=np.uint8),
                snpsmap=np.array([[0, 0, 0, 0, 100], [1, 0, 0, 0, 200]], dtype=np.uint32),
            )

    monkeypatch.setattr(common_mod, "SNPsImputer", DummyImputer)

    prepared = common_mod.get_numerical_input(
        DummyExtracter(),
        subsample=True,
        random_seed=13,
        impute_method="sample",
        log_level="INFO",
    )

    np.testing.assert_array_equal(prepared.matrix, np.array([[0, 1], [2, 2]], dtype=np.uint8))
    np.testing.assert_array_equal(
        observed["data"],
        np.array([[0, _MISSING_GENO], [2, 2]], dtype=np.uint8),
    )
    assert observed["names"] == ["a", "b"]
    assert observed["imap"] == {"pop1": ["a", "b"]}
    assert observed["impute_method"] == "sample"
    assert observed["quiet"] is True
    assert prepared.imputation.algorithm == "sample"
    assert prepared.imputation.imputed_snp_count == 1
    assert prepared.imputation.total_snps == 2
    assert prepared.imputation.imputed_genotype_count == 1
    assert prepared.imputation.total_genotypes == 4


def test_shared_numerical_input_normalizes_none_to_zero_fill() -> None:
    class DummyExtracter:
        snames = ["a", "b"]
        imap = {"pop1": ["a", "b"]}

        def get_view(self, *, subsample, random_seed, log_level):
            return SimpleNamespace(
                genos=np.array([[0, _MISSING_GENO], [_MISSING_GENO, 2]], dtype=np.uint8),
                snpsmap=np.array([[0, 0, 0, 0, 100], [1, 0, 0, 0, 200]], dtype=np.uint32),
            )

    prepared = common_mod.get_numerical_input(
        DummyExtracter(),
        subsample=False,
        random_seed=3,
        impute_method="none",
        log_level="INFO",
    )

    np.testing.assert_array_equal(prepared.matrix, np.array([[0, 0], [0, 2]], dtype=np.uint8))
    assert prepared.imputation.algorithm == "zero-fill"
    assert prepared.imputation.imputed_snp_fraction == pytest.approx(1.0)
    assert prepared.imputation.imputed_genotype_fraction == pytest.approx(0.5)
