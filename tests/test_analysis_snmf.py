import numpy as np
import pytest
from loguru import logger

import ipyrad2.analysis.methods.snmf as snmf_mod


def test_encode_genotypes_disjunctive_expands_three_state_blocks() -> None:
    matrix = np.array(
        [
            [0, 1, 2],
            [2, 1, 0],
        ],
        dtype=np.uint8,
    )

    encoded = snmf_mod._encode_genotypes_disjunctive(matrix)

    expected = np.array(
        [
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0],
        ]
    )
    assert np.array_equal(encoded, expected)


def test_membership_normalization_handles_zero_rows() -> None:
    weights = np.array(
        [
            [2.0, 2.0],
            [0.0, 0.0],
            [3.0, 1.0],
        ]
    )

    membership = snmf_mod._normalize_membership(weights)

    assert np.allclose(membership.sum(axis=1), 1.0)
    assert np.allclose(membership[0], [0.5, 0.5])
    assert np.allclose(membership[1], [0.5, 0.5])
    assert np.allclose(membership[2], [0.75, 0.25])


def test_genotype_frequency_normalization_and_allele_frequencies() -> None:
    components = np.array(
        [
            [1.0, 1.0, 2.0, 0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0, 3.0, 3.0, 0.0],
        ]
    )

    genotype_frequencies = snmf_mod._normalize_genotype_frequencies(components)
    allele_frequencies = snmf_mod._derive_allele_frequencies(genotype_frequencies)

    assert np.allclose(genotype_frequencies.sum(axis=2), 1.0)
    assert np.allclose(genotype_frequencies[0, 0], [0.25, 0.25, 0.5])
    assert np.allclose(genotype_frequencies[0, 1], [1 / 3, 1 / 3, 1 / 3])
    assert np.allclose(allele_frequencies[0], [0.625, 0.5])


def test_mean_cross_entropy_scores_true_genotypes() -> None:
    truth = np.array([0, 1, 2], dtype=np.uint8)
    predicted = np.array(
        [
            [0.9, 0.05, 0.05],
            [0.1, 0.8, 0.1],
            [0.05, 0.15, 0.8],
        ]
    )

    score = snmf_mod._mean_cross_entropy(truth, predicted)

    assert score == pytest.approx(
        -np.mean(np.log([0.9, 0.8, 0.8])),
        rel=1e-9,
    )


def test_fit_snmf_selects_best_initialization(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []

    def fake_fit_once(encoded, *, k, init, seed, alpha_w, alpha_h, l1_ratio):
        calls.append((encoded.shape, k, init, seed, alpha_w, alpha_h, l1_ratio))
        reconstruction_err = 5.0 if init == "nndsvda" else 1.0
        membership = np.full((encoded.shape[0], k), 1.0 / k)
        genotype_frequencies = np.full((k, encoded.shape[1] // 3, 3), 1.0 / 3.0)
        return snmf_mod.SNMFFit(
            membership=membership,
            genotype_frequencies=genotype_frequencies,
            allele_frequencies=snmf_mod._derive_allele_frequencies(genotype_frequencies),
            reconstruction_err=reconstruction_err,
            n_iter=10 if init == "nndsvda" else 7,
            hit_max_iter=False,
        )

    monkeypatch.setattr(snmf_mod, "_fit_snmf_once", fake_fit_once)

    fit = snmf_mod._fit_snmf(
        np.array(
            [
                [0, 1, 2],
                [2, 1, 0],
            ],
            dtype=np.uint8,
        ),
        k=2,
        seed=4,
        alpha_w=1e-4,
        alpha_h="same",
        l1_ratio=1.0,
        n_init=3,
    )

    assert len(calls) == 3
    assert calls[0][2] == "nndsvda"
    assert fit.reconstruction_err == 1.0
    assert fit.n_iter == 7


def test_log_k_convergence_warning_emits_one_summary_warning() -> None:
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        snmf_mod._log_k_convergence_warning(
            k=3,
            capped_fit_count=4,
            total_fit_count=6,
            selected_fit_hit_max_iter=True,
        )
    finally:
        logger.remove(sink_id)

    assert [msg.strip() for msg in messages] == [
        "sNMF K=3 convergence warning: 4/6 fits reached max_iter=3000 (selected fit also capped)."
    ]
