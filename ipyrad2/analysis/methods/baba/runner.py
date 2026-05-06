#!/usr/bin/env python

"""CLI runner for ABBA/BABA admixture statistics."""

from __future__ import annotations

import json
from itertools import combinations
from math import erfc, sqrt
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from ...extracters.snps_extracter import SNPsExtracter
from ..common import ensure_output_paths
from ..common import require_hdf5_input
from ....utils.exceptions import IPyradError
from .models import BabaRequest
from .models import BabaResult
from .models import QuartetDefinition


BABA_FLOAT_FORMAT = "%.8f"


def _coerce_positive_int(value, label: str) -> int:
    """Return one positive integer CLI option."""
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise IPyradError(f"{label} must be an integer.") from exc
    if parsed < 1:
        raise IPyradError(f"{label} must be >= 1.")
    return parsed


def _safe_ratio(numer: float, denom: float) -> float:
    """Return one guarded float ratio."""
    if np.isclose(denom, 0.0):
        return float("nan")
    return float(numer / denom)


def _normal_two_tailed_pvalue(z_score: float) -> float:
    """Return a two-tailed normal-approximation p-value."""
    if not np.isfinite(z_score):
        return float("nan")
    return float(erfc(abs(z_score) / sqrt(2.0)))


def _decode_tree_tips(node) -> list[str]:
    """Return descendant tip labels below one toytree node."""
    return [leaf.name for leaf in node.iter_leaves()]


def _read_tree_text(path: Path) -> str:
    """Return non-empty Newick text from one file."""
    path = Path(path).expanduser().absolute()
    if not path.exists():
        raise IPyradError(f"guide tree file does not exist: {path}")
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise IPyradError("guide tree text is empty.")
    if not text.endswith(";"):
        raise IPyradError("guide tree must end with ';'.")
    return text


def _parse_tree_or_error(tree_text: str):
    """Parse one rooted binary toytree object."""
    try:
        import toytree
    except ImportError as exc:
        raise IPyradError(
            "`ipyrad2 baba --tree` requires `toytree`. Install it first."
        ) from exc

    try:
        tree = toytree.tree(tree_text)
    except Exception as exc:  # pragma: no cover - parser-specific failure modes vary
        raise IPyradError(f"failed to parse guide tree: {exc}") from exc

    if not tree.is_rooted():
        raise IPyradError("guide tree must be rooted.")

    for node in tree.treenode.traverse():
        nchildren = len(node.children)
        if nchildren not in (0, 2):
            raise IPyradError("guide tree must be strictly bifurcating.")
    return tree


def _read_quartet_tests(path: Path) -> list[QuartetDefinition]:
    """Parse one manual quartet file."""
    path = Path(path).expanduser().absolute()
    if not path.exists():
        raise IPyradError(f"quartet file does not exist: {path}")

    quartets: list[QuartetDefinition] = []
    with open(path, encoding="utf-8") as infile:
        for lineno, line in enumerate(infile, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 4:
                raise IPyradError(
                    f"quartet file malformed at line {lineno}: expected `P1 P2 P3 P4`."
                )
            if len(set(parts)) != 4:
                raise IPyradError(
                    f"quartet file malformed at line {lineno}: quartet labels must be distinct."
                )
            quartets.append(
                QuartetDefinition(
                    source="tests",
                    p1=parts[0],
                    p2=parts[1],
                    p3=parts[2],
                    p4=parts[3],
                )
            )
    if not quartets:
        raise IPyradError("quartet file is empty.")
    return quartets


def _validate_manual_quartets(
    quartets: list[QuartetDefinition],
    allowed_labels: set[str],
) -> list[QuartetDefinition]:
    """Validate that all manual quartet labels resolve in the active namespace."""
    bad = sorted(
        {
            label
            for quartet in quartets
            for label in (quartet.p1, quartet.p2, quartet.p3, quartet.p4)
            if label not in allowed_labels
        }
    )
    if bad:
        raise IPyradError(
            "quartet labels were not found in the active sample namespace: "
            + ", ".join(bad[:10])
        )
    return quartets


def _resolve_rooted_quartet(pruned_tree) -> tuple[str, str, str, str] | None:
    """Return `P1,P2,P3,P4` from one rooted 4-tip tree or None if unresolved."""
    root = pruned_tree.treenode
    if len(root.children) != 2:
        return None

    left, right = root.children
    if left.is_leaf() == right.is_leaf():
        return None

    outgroup = left if left.is_leaf() else right
    ingroup = right if left.is_leaf() else left
    if len(ingroup.children) != 2:
        return None

    child_a, child_b = ingroup.children
    if child_a.is_leaf() == child_b.is_leaf():
        return None

    p3_node = child_a if child_a.is_leaf() else child_b
    pair_node = child_b if child_a.is_leaf() else child_a
    pair = sorted(_decode_tree_tips(pair_node))
    if len(pair) != 2:
        return None
    return pair[0], pair[1], p3_node.name, outgroup.name


def _expand_tree_quartets(tree) -> tuple[list[QuartetDefinition], int]:
    """Expand rooted valid quartets from one pruned guide tree."""
    labels = sorted(tree.get_tip_labels())
    quartets: list[QuartetDefinition] = []
    skipped = 0
    for combo in combinations(labels, 4):
        pruned = tree.mod.prune(*combo)
        resolved = _resolve_rooted_quartet(pruned)
        if resolved is None:
            skipped += 1
            continue
        quartets.append(
            QuartetDefinition(
                source="tree",
                p1=resolved[0],
                p2=resolved[1],
                p3=resolved[2],
                p4=resolved[3],
            )
        )
    if not quartets:
        raise IPyradError(
            "guide tree did not yield any rooted quartets with a unique outgroup."
        )
    return quartets, skipped


def _make_preview_extracter(
    *,
    data: Path,
    min_sample_coverage: int,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    cores: int,
):
    """Create a non-running SNP extracter to normalize sample selection inputs."""
    return SNPsExtracter(
        data=data,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=minmap,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
    )


def _resolve_quartets_and_selection(
    request: BabaRequest,
) -> tuple[
    list[QuartetDefinition],
    dict[str, list[str]] | None,
    list[str],
    dict[str, int] | None,
    str | None,
    dict[str, Any],
]:
    """Resolve quartets and return extraction inputs scoped to used labels only."""
    preview = _make_preview_extracter(
        data=request.data,
        min_sample_coverage=request.min_sample_coverage,
        imap=request.imap,
        minmap=request.minmap,
        exclude=list(request.exclude),
        include_reference=request.include_reference,
        cores=request.cores,
    )

    if preview.user_imap:
        available_labels = sorted(preview.imap)
    else:
        available_labels = sorted(preview.snames)
    available_set = set(available_labels)
    namespace = "population" if preview.user_imap else "sample"
    input_mode = "tests" if request.tests is not None else "tree"
    logger.info(
        "resolving baba {} in {} namespace with {} available labels",
        input_mode,
        namespace,
        len(available_labels),
    )

    resolved_tree_text = None
    tree_meta: dict[str, Any] = {
        "namespace": namespace,
        "input_mode": input_mode,
        "available_label_count": len(available_labels),
        "dropped_labels_not_in_tree": [],
        "skipped_balanced_quartets": 0,
        "tree_tip_count": 0,
        "used_label_count": 0,
    }
    if request.tests is not None:
        quartets = _validate_manual_quartets(
            _read_quartet_tests(request.tests),
            available_set,
        )
        logger.info(
            "resolved {} manual quartets against {} labels",
            len(quartets),
            namespace,
        )
    else:
        raw_tree = _parse_tree_or_error(_read_tree_text(request.tree))
        tree_tips = set(raw_tree.get_tip_labels())
        missing_from_data = sorted(tree_tips.difference(available_set))
        if missing_from_data:
            raise IPyradError(
                "guide tree contains names not found in the active sample namespace: "
                + ", ".join(missing_from_data[:10])
            )
        dropped = sorted(available_set.difference(tree_tips))
        tree_meta["dropped_labels_not_in_tree"] = dropped
        used_tree = raw_tree.mod.prune(*sorted(tree_tips.intersection(available_set)))
        if len(used_tree.get_tip_labels()) < 4:
            raise IPyradError("guide tree must contain at least 4 usable labels.")
        quartets, skipped = _expand_tree_quartets(used_tree)
        tree_meta["skipped_balanced_quartets"] = skipped
        tree_meta["tree_tip_count"] = len(used_tree.get_tip_labels())
        resolved_tree_text = used_tree.write()
        logger.info(
            "expanded {} rooted quartets from {} retained tree tips (skipped {}; dropped {})",
            len(quartets),
            tree_meta["tree_tip_count"],
            skipped,
            len(dropped),
        )

    used_labels = sorted(
        {
            label
            for quartet in quartets
            for label in (quartet.p1, quartet.p2, quartet.p3, quartet.p4)
        }
    )
    used_label_set = set(used_labels)
    tree_meta["used_label_count"] = len(used_labels)

    if preview.user_imap:
        run_imap = {label: preview.imap[label] for label in used_labels}
        run_exclude = sorted(set(preview.exclude))
        run_minmap = None if request.minmap is None else {label: int(preview.minmap[label]) for label in used_labels}
    else:
        run_imap = None
        run_exclude = sorted(set(preview.exclude).union(available_set.difference(used_label_set)))
        run_minmap = None

    return quartets, run_imap, run_exclude, run_minmap, resolved_tree_text, tree_meta


def _compute_group_arrays(
    genos: np.ndarray,
    sample_names: list[str],
    label_to_samples: dict[str, list[str]],
    random_seed: int | None,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, np.ndarray],
]:
    """Return per-label allele frequencies, coverage, and split frequencies."""
    freqs: dict[str, np.ndarray] = {}
    coverage: dict[str, np.ndarray] = {}
    split1: dict[str, np.ndarray] = {}
    split2: dict[str, np.ndarray] = {}
    name_to_index = {name: idx for idx, name in enumerate(sample_names)}
    seed_base = (
        (0 if random_seed is None else int(random_seed) * 1009)
        + int(genos.shape[1])
        + int(genos.shape[0]) * 13
    )

    for label, members in label_to_samples.items():
        indices = [name_to_index[name] for name in members if name in name_to_index]
        if not indices:
            raise IPyradError(f"label `{label}` has no retained samples after filtering.")
        subset = genos[indices]
        calls = subset != 255
        sample_cov = np.sum(calls, axis=0).astype(np.int32, copy=False)
        alt_sum = np.sum(np.where(calls, subset, 0), axis=0, dtype=np.int64)
        denom = sample_cov.astype(float) * 2.0
        freq = np.divide(
            alt_sum.astype(float),
            denom,
            out=np.full(subset.shape[1], np.nan, dtype=float),
            where=denom > 0,
        )
        sample_afs = np.where(calls, subset.astype(float) / 2.0, np.nan)
        split1_vals = np.full(subset.shape[1], np.nan, dtype=float)
        split2_vals = np.full(subset.shape[1], np.nan, dtype=float)
        label_seed = seed_base + sum(ord(char) for char in label) * 17 + len(indices) * 101
        rng = np.random.default_rng(label_seed)
        for site_idx in range(subset.shape[1]):
            observed = sample_afs[:, site_idx]
            observed = observed[np.isfinite(observed)]
            if observed.size == 0:
                continue
            if observed.size == 1:
                split1_vals[site_idx] = observed[0]
                split2_vals[site_idx] = observed[0]
                continue
            draw1 = observed[rng.integers(0, observed.size, size=observed.size)]
            draw2 = observed[rng.integers(0, observed.size, size=observed.size)]
            split1_vals[site_idx] = float(np.mean(draw1))
            split2_vals[site_idx] = float(np.mean(draw2))
        freqs[label] = freq
        coverage[label] = sample_cov
        split1[label] = split1_vals
        split2[label] = split2_vals
    return freqs, coverage, split1, split2


def _site_pattern_arrays(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    p4: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return site-wise ABBA/BABA/BBAA weights from derived-allele frequencies."""
    abba = ((1.0 - p1) * p2 * p3 * (1.0 - p4)) + (p1 * (1.0 - p2) * (1.0 - p3) * p4)
    baba = (p1 * (1.0 - p2) * p3 * (1.0 - p4)) + ((1.0 - p1) * p2 * (1.0 - p3) * p4)
    bbaa = (p1 * p2 * (1.0 - p3) * (1.0 - p4)) + ((1.0 - p1) * (1.0 - p2) * p3 * p4)
    return abba, baba, bbaa


def _f_g_denom_per_variant(
    p1: np.ndarray,
    p3a: np.ndarray,
    p3b: np.ndarray,
    p4: np.ndarray,
) -> np.ndarray:
    """Return one Dsuite-style per-site f_G denominator."""
    return ((1.0 - p1) * p3a * p3b * (1.0 - p4)) - (p1 * (1.0 - p3a) * p3b * (1.0 - p4))


def _locus_weights(
    locs: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return unique locus IDs and per-locus aggregated values."""
    unique_locs, inverse = np.unique(locs, return_inverse=True)
    out = np.zeros(unique_locs.shape[0], dtype=float)
    np.add.at(out, inverse, values)
    return unique_locs, out


def _physical_block_partition(
    snpsmap: np.ndarray,
    block_bp: int,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Return per-site physical-block IDs and ordered metadata."""
    labels = np.empty(snpsmap.shape[0], dtype=np.int32)
    metadata: list[dict[str, Any]] = []
    key_to_idx: dict[tuple[int, int], int] = {}

    for idx, row in enumerate(snpsmap):
        scaff = int(row[3])
        block_idx = int(row[4]) // block_bp
        key = (scaff, block_idx)
        if key not in key_to_idx:
            key_to_idx[key] = len(metadata)
            metadata.append(
                {
                    "block_label": f"scaff{scaff}:{block_idx}",
                    "scaff": scaff,
                    "block_start": block_idx * block_bp,
                    "block_end": (block_idx + 1) * block_bp,
                }
            )
        labels[idx] = key_to_idx[key]
    return labels, metadata


def _locus_block_partition(
    snpsmap: np.ndarray,
    loci_per_block: int,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Return per-site consecutive-locus block IDs and ordered metadata."""
    locs = snpsmap[:, 0].astype(np.int64, copy=False)
    unique_locs = np.unique(locs)
    block_by_loc = {
        int(loc): idx // loci_per_block for idx, loc in enumerate(unique_locs.tolist())
    }
    labels = np.array([block_by_loc[int(loc)] for loc in locs], dtype=np.int32)
    metadata: list[dict[str, Any]] = []
    nblocks = int(labels.max()) + 1 if labels.size else 0
    for block_idx in range(nblocks):
        block_locs = unique_locs[block_idx * loci_per_block:(block_idx + 1) * loci_per_block]
        metadata.append(
            {
                "block_label": f"locus_block:{block_idx}",
                "scaff": -1,
                "block_start": int(block_locs[0]),
                "block_end": int(block_locs[-1]),
            }
        )
    return labels, metadata


def _ratio_sum(
    numer: np.ndarray,
    denom: np.ndarray,
) -> float:
    """Return one ratio of summed sitewise numerator and denominator values."""
    return _safe_ratio(float(np.sum(numer)), float(np.sum(denom)))


def _jackknife_ratio_stats(
    numer: np.ndarray,
    denom: np.ndarray,
    block_ids: np.ndarray,
) -> tuple[float, dict[str, Any]]:
    """Return one full estimate plus delete-one-block jackknife summaries."""
    unique_blocks = np.unique(block_ids)
    full_estimate = _ratio_sum(numer, denom)
    if unique_blocks.size < 2:
        return full_estimate, {
            "resampling_units": int(unique_blocks.size),
            "estimate_se": float("nan"),
            "z_score": float("nan"),
            "p_value": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }

    estimates = []
    for block in unique_blocks.tolist():
        keep = block_ids != block
        estimates.append(_ratio_sum(numer[keep], denom[keep]))
    est = np.asarray(estimates, dtype=float)
    est = est[np.isfinite(est)]
    if est.size < 2:
        return full_estimate, {
            "resampling_units": int(unique_blocks.size),
            "estimate_se": float("nan"),
            "z_score": float("nan"),
            "p_value": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }
    mean_est = float(np.mean(est))
    se = float(sqrt(((est.size - 1) / est.size) * np.sum((est - mean_est) ** 2)))
    z_score = full_estimate / se if np.isfinite(full_estimate) and se > 0 else float("nan")
    return full_estimate, {
        "resampling_units": int(unique_blocks.size),
        "estimate_se": se,
        "z_score": z_score,
        "p_value": _normal_two_tailed_pvalue(z_score),
        "ci_low": full_estimate - 1.96 * se if np.isfinite(full_estimate) and np.isfinite(se) else float("nan"),
        "ci_high": full_estimate + 1.96 * se if np.isfinite(full_estimate) and np.isfinite(se) else float("nan"),
    }


def _bootstrap_ratio_stats(
    numer: np.ndarray,
    denom: np.ndarray,
    snpsmap: np.ndarray,
    replicates: int,
    rng: np.random.Generator,
) -> tuple[float, dict[str, Any]]:
    """Return one full estimate plus locus-bootstrap summaries."""
    locs = snpsmap[:, 0].astype(np.int64, copy=False)
    unique_locs, numer_by_loc = _locus_weights(locs, numer)
    _, denom_by_loc = _locus_weights(locs, denom)
    nloci = unique_locs.size
    full_estimate = _ratio_sum(numer, denom)
    if nloci == 0 or replicates < 1:
        return full_estimate, {
            "resampling_units": int(nloci),
            "resampling_replicates": int(replicates),
            "resample_mean": float("nan"),
            "estimate_se": float("nan"),
            "z_score": float("nan"),
            "p_value": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }

    estimates = np.empty(replicates, dtype=float)
    for rep in range(replicates):
        draw = rng.integers(0, nloci, size=nloci)
        counts = np.bincount(draw, minlength=nloci).astype(float)
        numer_sum = float(np.dot(numer_by_loc, counts))
        denom_sum = float(np.dot(denom_by_loc, counts))
        estimates[rep] = _safe_ratio(numer_sum, denom_sum)

    valid = estimates[np.isfinite(estimates)]
    if valid.size == 0:
        return full_estimate, {
            "resampling_units": int(nloci),
            "resampling_replicates": int(replicates),
            "resample_mean": float("nan"),
            "estimate_se": float("nan"),
            "z_score": float("nan"),
            "p_value": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
        }
    se = float(np.std(valid, ddof=1)) if valid.size > 1 else float("nan")
    z_score = full_estimate / se if np.isfinite(full_estimate) and np.isfinite(se) and se > 0 else float("nan")
    return full_estimate, {
        "resampling_units": int(nloci),
        "resampling_replicates": int(replicates),
        "resample_mean": float(np.mean(valid)),
        "estimate_se": se,
        "z_score": z_score,
        "p_value": _normal_two_tailed_pvalue(z_score),
        "ci_low": float(np.percentile(valid, 2.5)),
        "ci_high": float(np.percentile(valid, 97.5)),
    }


def _resolve_resampling_layout(
    used_map: np.ndarray,
    request: BabaRequest,
) -> tuple[str, str, np.ndarray, list[dict[str, Any]], np.ndarray, list[dict[str, Any]]]:
    """Return resampling mode plus chosen and physical block layouts."""
    if used_map.shape[0]:
        physical_ids, physical_meta = _physical_block_partition(
            used_map,
            request.jackknife_block_bp,
        )
        physical_n = int(np.unique(physical_ids).size)
    else:
        physical_ids = np.empty(0, dtype=np.int32)
        physical_meta = []
        physical_n = 0

    resampling_mode = request.resampling
    if resampling_mode == "auto":
        resampling_mode = "jackknife" if physical_n >= 20 else "bootstrap"

    if resampling_mode == "jackknife":
        if physical_n >= 20:
            return (
                "jackknife",
                "physical_block",
                physical_ids,
                physical_meta,
                physical_ids,
                physical_meta,
            )
        block_ids, block_meta = _locus_block_partition(
            used_map,
            request.jackknife_block_loci,
        )
        return (
            "jackknife",
            "locus_block",
            block_ids,
            block_meta,
            physical_ids,
            physical_meta,
        )

    if resampling_mode == "bootstrap":
        return (
            "bootstrap",
            "locus",
            np.empty(0, dtype=np.int32),
            [],
            physical_ids,
            physical_meta,
        )

    return (
        "none",
        "none",
        np.empty(0, dtype=np.int32),
        [],
        physical_ids,
        physical_meta,
    )


def _apply_ratio_stats(
    numer: np.ndarray,
    denom: np.ndarray,
    used_map: np.ndarray,
    request: BabaRequest,
    mode: str,
    unit: str,
    block_ids: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, dict[str, Any]]:
    """Return one ratio estimate plus resampling summary for the requested mode."""
    if mode == "jackknife":
        estimate, stats = _jackknife_ratio_stats(numer, denom, block_ids)
        stats["resampling_mode"] = "jackknife"
        stats["resampling_unit"] = unit
        return estimate, stats
    if mode == "bootstrap":
        estimate, stats = _bootstrap_ratio_stats(
            numer,
            denom,
            used_map,
            request.bootstrap_replicates,
            rng,
        )
        stats["resampling_mode"] = "bootstrap"
        stats["resampling_unit"] = unit
        return estimate, stats
    return _ratio_sum(numer, denom), {
        "resampling_mode": "none",
        "resampling_unit": "none",
        "resampling_units": 0,
        "estimate_se": float("nan"),
        "z_score": float("nan"),
        "p_value": float("nan"),
        "ci_low": float("nan"),
        "ci_high": float("nan"),
    }


def _rooted_orientation_payload(
    quartet: QuartetDefinition,
    mode: str,
    abba: np.ndarray,
    baba: np.ndarray,
    bbaa: np.ndarray,
    fg1: np.ndarray,
    fg1_reversed: np.ndarray,
    fg2: np.ndarray,
    fg2_reversed: np.ndarray,
    fg3: np.ndarray,
    fg3_reversed: np.ndarray,
) -> dict[str, Any]:
    """Return ordered labels and sitewise arrays for one rooted trio orientation."""
    trio = (quartet.p1, quartet.p2, quartet.p3, quartet.p4)
    if mode == "d1":
        positive = float(np.sum(abba) - np.sum(baba)) >= 0.0
        if positive:
            return {
                "p1": trio[0],
                "p2": trio[1],
                "p3": trio[2],
                "p4": trio[3],
                "abba": abba,
                "baba": baba,
                "bbaa": bbaa,
                "numer": abba - baba,
                "denom": abba + baba,
                "fg_denom": fg1,
                "mode": "d1",
            }
        return {
            "p1": trio[1],
            "p2": trio[0],
            "p3": trio[2],
            "p4": trio[3],
            "abba": baba,
            "baba": abba,
            "bbaa": bbaa,
            "numer": baba - abba,
            "denom": abba + baba,
            "fg_denom": fg1_reversed,
            "mode": "d1",
        }

    if mode == "d2":
        positive = float(np.sum(abba) - np.sum(bbaa)) >= 0.0
        if positive:
            return {
                "p1": trio[0],
                "p2": trio[2],
                "p3": trio[1],
                "p4": trio[3],
                "abba": abba,
                "baba": bbaa,
                "bbaa": baba,
                "numer": abba - bbaa,
                "denom": abba + bbaa,
                "fg_denom": fg2,
                "mode": "d2",
            }
        return {
            "p1": trio[2],
            "p2": trio[0],
            "p3": trio[1],
            "p4": trio[3],
            "abba": bbaa,
            "baba": abba,
            "bbaa": baba,
            "numer": bbaa - abba,
            "denom": abba + bbaa,
            "fg_denom": fg2_reversed,
            "mode": "d2",
        }

    positive = float(np.sum(bbaa) - np.sum(baba)) >= 0.0
    if positive:
        return {
            "p1": trio[2],
            "p2": trio[1],
            "p3": trio[0],
            "p4": trio[3],
            "abba": bbaa,
            "baba": baba,
            "bbaa": abba,
            "numer": bbaa - baba,
            "denom": bbaa + baba,
            "fg_denom": fg3,
            "mode": "d3",
        }
    return {
        "p1": trio[1],
        "p2": trio[2],
        "p3": trio[0],
        "p4": trio[3],
        "abba": baba,
        "baba": bbaa,
        "bbaa": abba,
        "numer": baba - bbaa,
        "denom": bbaa + baba,
        "fg_denom": fg3_reversed,
        "mode": "d3",
    }


def _select_bbaa_mode(
    abba_sum: float,
    baba_sum: float,
    bbaa_sum: float,
) -> str:
    """Return the Dsuite trio mode implied by the dominant concordant topology."""
    if bbaa_sum >= baba_sum and bbaa_sum >= abba_sum:
        return "d1"
    if baba_sum >= bbaa_sum and baba_sum >= abba_sum:
        return "d2"
    return "d3"


def _select_dmin_mode(
    d1: float,
    d2: float,
    d3: float,
) -> str:
    """Return the Dsuite trio mode with the smallest absolute D."""
    values = {"d1": abs(d1), "d2": abs(d2), "d3": abs(d3)}
    return min(values, key=values.get)


def _build_block_rows(
    quartet_idx: int,
    quartet: QuartetDefinition,
    abba: np.ndarray,
    baba: np.ndarray,
    bbaa: np.ndarray,
    block_ids: np.ndarray,
    metadata: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return per-block site-pattern summary rows."""
    rows: list[dict[str, Any]] = []
    for block_idx in np.unique(block_ids).tolist():
        keep = block_ids == block_idx
        abba_sum = float(np.sum(abba[keep]))
        baba_sum = float(np.sum(baba[keep]))
        bbaa_sum = float(np.sum(bbaa[keep]))
        meta = metadata[int(block_idx)]
        rows.append(
            {
                "quartet_id": quartet_idx,
                "p1": quartet.p1,
                "p2": quartet.p2,
                "p3": quartet.p3,
                "p4": quartet.p4,
                "block_index": int(block_idx),
                "block_label": meta["block_label"],
                "scaff": meta["scaff"],
                "block_start": meta["block_start"],
                "block_end": meta["block_end"],
                "n_sites": int(np.sum(keep)),
                "abba": abba_sum,
                "baba": baba_sum,
                "bbaa": bbaa_sum,
                "abba_minus_baba": abba_sum - baba_sum,
                "d_stat": _safe_ratio(abba_sum - baba_sum, abba_sum + baba_sum),
            }
        )
    return rows


def _add_clustering_stats(row: dict[str, Any], block_rows: list[dict[str, Any]]) -> None:
    """Attach quartet-level block clustering summaries."""
    if not block_rows:
        row.update(
            {
                "n_blocks": 0,
                "mean_block_d": float("nan"),
                "sd_block_d": float("nan"),
                "fraction_blocks_abba_gt_baba": float("nan"),
                "longest_same_sign_run": 0,
                "sign_switches": 0,
            }
        )
        return

    ordered = sorted(block_rows, key=lambda item: (item["scaff"], item["block_start"], item["block_end"]))
    dvals = np.asarray([item["d_stat"] for item in ordered], dtype=float)
    valid = dvals[np.isfinite(dvals)]
    signs = [0 if not np.isfinite(value) or np.isclose(value, 0.0) else (1 if value > 0 else -1) for value in dvals.tolist()]
    longest = 0
    switches = 0
    prev = 0
    run = 0
    for sign in signs:
        if sign == 0:
            run = 0
            continue
        if prev == 0 or sign == prev:
            run += 1
        else:
            switches += 1
            run = 1
        prev = sign
        longest = max(longest, run)

    row.update(
        {
            "n_blocks": len(block_rows),
            "mean_block_d": float(np.mean(valid)) if valid.size else float("nan"),
            "sd_block_d": float(np.std(valid, ddof=1)) if valid.size > 1 else float("nan"),
            "fraction_blocks_abba_gt_baba": float(
                np.mean([item["abba"] > item["baba"] for item in ordered])
            ),
            "longest_same_sign_run": int(longest),
            "sign_switches": int(switches),
        }
    )


def _summarize_quartet(
    quartet_idx: int,
    quartet: QuartetDefinition,
    freqs: dict[str, np.ndarray],
    coverage: dict[str, np.ndarray],
    split1: dict[str, np.ndarray],
    split2: dict[str, np.ndarray],
    snpsmap: np.ndarray,
    request: BabaRequest,
    quartet_minmap: dict[str, int] | None,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    """Return one quartet summary row, rooted rows, and optional block rows."""
    p1 = freqs[quartet.p1]
    p2 = freqs[quartet.p2]
    p3 = freqs[quartet.p3]
    p4 = freqs[quartet.p4]
    finite_mask = np.isfinite(p1) & np.isfinite(p2) & np.isfinite(p3) & np.isfinite(p4)

    thresholds = {label: 1 for label in (quartet.p1, quartet.p2, quartet.p3, quartet.p4)}
    if quartet_minmap is not None:
        for label in thresholds:
            thresholds[label] = max(1, int(quartet_minmap[label]))
    cov_mask = (
        (coverage[quartet.p1] >= thresholds[quartet.p1])
        & (coverage[quartet.p2] >= thresholds[quartet.p2])
        & (coverage[quartet.p3] >= thresholds[quartet.p3])
        & (coverage[quartet.p4] >= thresholds[quartet.p4])
    )
    mask = finite_mask & cov_mask

    abba_all, baba_all, bbaa_all = _site_pattern_arrays(p1, p2, p3, p4)
    abba = abba_all[mask]
    baba = baba_all[mask]
    bbaa = bbaa_all[mask]
    used_map = snpsmap[mask]
    p1s1 = split1[quartet.p1][mask]
    p1s2 = split2[quartet.p1][mask]
    p2s1 = split1[quartet.p2][mask]
    p2s2 = split2[quartet.p2][mask]
    p3s1 = split1[quartet.p3][mask]
    p3s2 = split2[quartet.p3][mask]
    mp1 = p1[mask]
    mp2 = p2[mask]
    mp3 = p3[mask]
    mp4 = p4[mask]

    fg1 = _f_g_denom_per_variant(mp1, p3s1, p3s2, mp4)
    fg1_reversed = _f_g_denom_per_variant(mp2, p3s1, p3s2, mp4)
    fg2 = _f_g_denom_per_variant(mp1, p2s1, p2s2, mp4)
    fg2_reversed = _f_g_denom_per_variant(mp3, p2s1, p2s2, mp4)
    fg3 = _f_g_denom_per_variant(mp3, p1s1, p1s2, mp4)
    fg3_reversed = _f_g_denom_per_variant(mp2, p1s1, p1s2, mp4)
    if mp4.size:
        fg1 += _f_g_denom_per_variant(1.0 - mp1, 1.0 - p3s1, 1.0 - p3s2, 1.0 - mp4)
        fg1_reversed += _f_g_denom_per_variant(1.0 - mp2, 1.0 - p3s1, 1.0 - p3s2, 1.0 - mp4)
        fg2 += _f_g_denom_per_variant(1.0 - mp1, 1.0 - p2s1, 1.0 - p2s2, 1.0 - mp4)
        fg2_reversed += _f_g_denom_per_variant(1.0 - mp3, 1.0 - p2s1, 1.0 - p2s2, 1.0 - mp4)
        fg3 += _f_g_denom_per_variant(1.0 - mp3, 1.0 - p1s1, 1.0 - p1s2, 1.0 - mp4)
        fg3_reversed += _f_g_denom_per_variant(1.0 - mp2, 1.0 - p1s1, 1.0 - p1s2, 1.0 - mp4)

    abba_sum = float(np.sum(abba))
    baba_sum = float(np.sum(baba))
    bbaa_sum = float(np.sum(bbaa))
    resampling_mode, resampling_unit, block_ids, _block_meta, physical_ids, physical_meta = _resolve_resampling_layout(
        used_map,
        request,
    )
    d_estimate, d_stats = _apply_ratio_stats(
        abba - baba,
        abba + baba,
        used_map,
        request,
        resampling_mode,
        resampling_unit,
        block_ids,
        rng,
    )
    row: dict[str, Any] = {
        "quartet_id": quartet_idx,
        "source": quartet.source,
        "p1": quartet.p1,
        "p2": quartet.p2,
        "p3": quartet.p3,
        "p4": quartet.p4,
        "n_sites_tested": int(mask.sum()),
        "n_loci_tested": int(np.unique(used_map[:, 0]).size if used_map.size else 0),
        "abba": abba_sum,
        "baba": baba_sum,
        "bbaa": bbaa_sum,
        "abba_minus_baba": abba_sum - baba_sum,
        "d_stat": d_estimate,
        "abba_over_bbaa": _safe_ratio(abba_sum, bbaa_sum),
        "baba_over_bbaa": _safe_ratio(baba_sum, bbaa_sum),
        "discordant_over_bbaa": _safe_ratio(abba_sum + baba_sum, bbaa_sum),
        "resampling_mode": d_stats["resampling_mode"],
        "resampling_unit": d_stats["resampling_unit"],
        "resampling_units": d_stats["resampling_units"],
        "estimate_se": d_stats["estimate_se"],
        "z_score": d_stats["z_score"],
        "p_value": d_stats["p_value"],
        "ci_low": d_stats["ci_low"],
        "ci_high": d_stats["ci_high"],
    }
    if "resampling_replicates" in d_stats:
        row["resampling_replicates"] = d_stats["resampling_replicates"]

    block_rows: list[dict[str, Any]] = []
    if request.write_block_table or request.clustering_stats:
        if physical_meta and physical_ids.size:
            block_rows = _build_block_rows(
                quartet_idx,
                quartet,
                abba,
                baba,
                bbaa,
                physical_ids,
                physical_meta,
            )
        if request.clustering_stats:
            _add_clustering_stats(row, block_rows)

    d1 = _ratio_sum(abba - baba, abba + baba)
    d2 = _ratio_sum(abba - bbaa, abba + bbaa)
    d3 = _ratio_sum(bbaa - baba, bbaa + baba)
    rooted_modes = [
        ("input", "d1"),
        ("bbaa", _select_bbaa_mode(abba_sum, baba_sum, bbaa_sum)),
        ("dmin", _select_dmin_mode(d1, d2, d3)),
    ]
    if quartet.source == "tree":
        rooted_modes.append(("tree", "d1"))

    rooted_rows: list[dict[str, Any]] = []
    for orientation, mode in rooted_modes:
        payload = _rooted_orientation_payload(
            quartet,
            mode,
            abba,
            baba,
            bbaa,
            fg1,
            fg1_reversed,
            fg2,
            fg2_reversed,
            fg3,
            fg3_reversed,
        )
        rooted_d, rooted_d_stats = _apply_ratio_stats(
            payload["numer"],
            payload["denom"],
            used_map,
            request,
            resampling_mode,
            resampling_unit,
            block_ids,
            rng,
        )
        rooted_fg, rooted_fg_stats = _apply_ratio_stats(
            payload["numer"],
            payload["fg_denom"],
            used_map,
            request,
            resampling_mode,
            resampling_unit,
            block_ids,
            rng,
        )
        if np.isfinite(rooted_fg):
            rooted_fg = min(1.0, max(0.0, rooted_fg))
        fg_ci_low = rooted_fg_stats["ci_low"]
        fg_ci_high = rooted_fg_stats["ci_high"]
        if np.isfinite(fg_ci_low):
            fg_ci_low = min(1.0, max(0.0, fg_ci_low))
        if np.isfinite(fg_ci_high):
            fg_ci_high = min(1.0, max(0.0, fg_ci_high))
        rooted_row = {
            "quartet_id": quartet_idx,
            "source": quartet.source,
            "orientation": orientation,
            "arrangement": payload["mode"],
            "p1": payload["p1"],
            "p2": payload["p2"],
            "p3": payload["p3"],
            "p4": payload["p4"],
            "n_sites_tested": int(mask.sum()),
            "n_loci_tested": int(np.unique(used_map[:, 0]).size if used_map.size else 0),
            "bbaa": float(np.sum(payload["bbaa"])),
            "abba": float(np.sum(payload["abba"])),
            "baba": float(np.sum(payload["baba"])),
            "d_stat": rooted_d,
            "d_se": rooted_d_stats["estimate_se"],
            "d_z_score": rooted_d_stats["z_score"],
            "d_p_value": rooted_d_stats["p_value"],
            "d_ci_low": rooted_d_stats["ci_low"],
            "d_ci_high": rooted_d_stats["ci_high"],
            "f_g": rooted_fg,
            "f_g_se": rooted_fg_stats["estimate_se"],
            "f_g_ci_low": fg_ci_low,
            "f_g_ci_high": fg_ci_high,
            "resampling_mode": rooted_d_stats["resampling_mode"],
            "resampling_unit": rooted_d_stats["resampling_unit"],
            "resampling_units": rooted_d_stats["resampling_units"],
        }
        if "resampling_replicates" in rooted_d_stats:
            rooted_row["resampling_replicates"] = rooted_d_stats["resampling_replicates"]
        rooted_rows.append(rooted_row)

    return row, rooted_rows, block_rows


def _write_table(path: Path, table: pd.DataFrame) -> None:
    """Write one TSV with a shared readable float format."""
    table.to_csv(path, sep="\t", index=False, float_format=BABA_FLOAT_FORMAT)


def _result_tables(result: BabaResult) -> dict[str, pd.DataFrame | None]:
    """Return ordered result tables, including optional outputs."""
    return {
        "quartets": result.quartets,
        "rooted": result.rooted,
        "resolved_tests": result.resolved_tests,
        "f_branch": result.f_branch,
        "f_branch_matrix": result.f_branch_matrix,
        "f_branch_z": result.f_branch_z,
        "f_branch_p": result.f_branch_p,
        "blocks": result.blocks,
    }


def _build_output_paths(result: BabaResult, request: BabaRequest) -> dict[str, Path]:
    """Return output paths for one baba result."""
    paths = {
        "quartets": request.outdir / f"{request.name}.quartets.tsv",
        "rooted": request.outdir / f"{request.name}.rooted.tsv",
        "resolved_tests": request.outdir / f"{request.name}.tests.resolved.tsv",
        "manifest": request.outdir / f"{request.name}.manifest.txt",
        "summary_json": request.outdir / f"{request.name}.summary.json",
    }
    if result.tree_text is not None:
        paths["tree"] = request.outdir / f"{request.name}.tree.used.nwk"
    optional_suffixes = {
        "f_branch": "f_branch.tsv",
        "f_branch_matrix": "f_branch.matrix.tsv",
        "f_branch_z": "f_branch.z.tsv",
        "f_branch_p": "f_branch.p.tsv",
        "blocks": "blocks.tsv",
    }
    for key, table in _result_tables(result).items():
        if key not in optional_suffixes or table is None or table.empty:
            continue
        paths[key] = request.outdir / f"{request.name}.{optional_suffixes[key]}"
    return paths


def _render_table_preview(
    table: pd.DataFrame | None,
    columns: list[str],
    limit: int = 10,
) -> str:
    """Return one compact preview table for the manifest."""
    if table is None or table.empty:
        return "no rows written\n"
    preview = table.loc[:, [col for col in columns if col in table.columns]].head(limit)
    return preview.to_string(index=False, float_format=lambda value: f"{value:.8f}") + "\n"


def _normalize_json_value(value: Any) -> Any:
    """Convert numpy/pandas values into JSON-safe scalars."""
    if isinstance(value, dict):
        return {str(key): _normalize_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalize_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    return value


def _table_to_json_records(table: pd.DataFrame | None) -> list[dict[str, Any]]:
    """Return one table as JSON-safe records."""
    if table is None or table.empty:
        return []
    normalized = table.astype(object).where(pd.notnull(table), None)
    return [
        {str(key): _normalize_json_value(value) for key, value in record.items()}
        for record in normalized.to_dict(orient="records")
    ]


def _build_results_summary(result: BabaResult) -> dict[str, Any]:
    """Return summary counts for written result tables."""
    summary = {
        "quartets_rows": int(len(result.quartets.index)),
        "rooted_rows": int(len(result.rooted.index)),
        "resolved_tests_rows": int(len(result.resolved_tests.index)),
        "f_branch_rows": int(len(result.f_branch.index)) if result.f_branch is not None else 0,
        "f_branch_matrix_rows": int(len(result.f_branch_matrix.index)) if result.f_branch_matrix is not None else 0,
        "f_branch_z_rows": int(len(result.f_branch_z.index)) if result.f_branch_z is not None else 0,
        "f_branch_p_rows": int(len(result.f_branch_p.index)) if result.f_branch_p is not None else 0,
        "blocks_rows": int(len(result.blocks.index)) if result.blocks is not None else 0,
    }
    if result.rooted is None or result.rooted.empty or "orientation" not in result.rooted:
        summary["rooted_orientation_counts"] = {}
    else:
        counts = result.rooted["orientation"].value_counts().sort_index()
        summary["rooted_orientation_counts"] = {
            str(key): int(value) for key, value in counts.items()
        }
    return summary


def _build_data_summary(
    request: BabaRequest,
    extracter: SNPsExtracter,
    quartets: list[QuartetDefinition],
    tree_meta: dict[str, Any],
    result: BabaResult,
) -> dict[str, Any]:
    """Return one shared data summary payload."""
    scaffolds = (
        int(np.unique(extracter.snpsmap[:, 3]).size)
        if getattr(extracter, "snpsmap", None) is not None and extracter.snpsmap.size
        else 0
    )
    return {
        "namespace": tree_meta["namespace"],
        "input_mode": tree_meta["input_mode"],
        "samples_retained": int(len(extracter.snames)),
        "labels_retained": int(tree_meta["used_label_count"]),
        "available_labels": int(tree_meta["available_label_count"]),
        "post_filter_snps": int(extracter.stats["post_filter_snps"]),
        "post_filter_loci": int(extracter.stats["post_filter_snp_containing_linkage_blocks"]),
        "scaffolds_retained": scaffolds,
        "quartets_resolved": int(len(quartets)),
        "rooted_rows_written": int(len(result.rooted.index)),
        "tree_tip_count": int(tree_meta["tree_tip_count"]),
        "tree_dropped_labels_not_in_tree": tree_meta["dropped_labels_not_in_tree"] or [],
        "tree_skipped_balanced_quartets": int(tree_meta["skipped_balanced_quartets"]),
    }


def _build_summary_json(
    request: BabaRequest,
    extracter: SNPsExtracter,
    quartets: list[QuartetDefinition],
    tree_meta: dict[str, Any],
    result: BabaResult,
    output_paths: dict[str, Path],
) -> str:
    """Return one bundled machine-readable baba summary."""
    data_summary = _build_data_summary(request, extracter, quartets, tree_meta, result)
    results_summary = _build_results_summary(result)
    payload: dict[str, Any] = {
        "tool": "baba",
        "inputs": {
            "data": request.data,
            "tests": request.tests,
            "tree": request.tree,
            "imap": request.imap,
            "minmap": request.minmap,
            "exclude": list(request.exclude),
            "include_reference": request.include_reference,
        },
        "run": {
            "name": request.name,
            "outdir": request.outdir,
            "min_sample_coverage": request.min_sample_coverage,
            "resampling": request.resampling,
            "bootstrap_replicates": request.bootstrap_replicates,
            "jackknife_block_bp": request.jackknife_block_bp,
            "jackknife_block_loci": request.jackknife_block_loci,
            "seed": request.seed,
            "f_branch": request.f_branch,
            "f_branch_p_threshold": request.f_branch_p_threshold,
            "write_block_table": request.write_block_table,
            "clustering_stats": request.clustering_stats,
            "cores": request.cores,
            "log_level": request.log_level,
        },
        "data_summary": data_summary,
        "results_summary": results_summary,
        "tables": {
            key: _table_to_json_records(table)
            for key, table in _result_tables(result).items()
        },
        "outputs": {key: str(path) for key, path in output_paths.items()},
    }
    if request.logged_command:
        payload["command"] = request.logged_command
    return json.dumps(_normalize_json_value(payload), indent=2) + "\n"


def _build_f_branch_matrix(
    rows: list[dict[str, Any]],
    tips: list[str],
    value_key: str,
) -> pd.DataFrame:
    """Return one branch-by-recipient matrix from long-form branch rows."""
    grouped: dict[str, dict[str, Any]] = {}
    for row in rows:
        branch_id = str(row["branch_id"])
        if branch_id not in grouped:
            grouped[branch_id] = {
                "branch": branch_id,
                "branch_descendants": row["branch_descendants"],
            }
            for tip in tips:
                grouped[branch_id][tip] = float("nan")
        grouped[branch_id][str(row["recipient"])] = row[value_key]
    return pd.DataFrame(list(grouped.values()))


def _build_f_branch_outputs(
    tree_text: str,
    rooted: pd.DataFrame,
    p_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return long-form and matrix f-branch outputs from rooted tree rows."""
    tree = _parse_tree_or_error(tree_text)
    tree_rows = rooted[rooted["orientation"] == "tree"].copy()
    if tree_rows.empty:
        raise IPyradError("`--f-branch` requires rooted tree-orientation rows.")
    tips = list(tree.get_tip_labels())
    tip_set = set(tips)
    raw_map: dict[tuple[str, str], list[tuple[str, float, float]]] = {}
    thresholded_map: dict[tuple[str, str], list[tuple[str, float, float]]] = {}
    for record in tree_rows.itertuples(index=False):
        if record.p1 not in tip_set or record.p2 not in tip_set or record.p3 not in tip_set:
            continue
        raw_map.setdefault((record.p1, record.p3), []).append(
            (record.p2, float(record.f_g), float(record.d_z_score))
        )
        thresholded_value = float(record.f_g) if float(record.d_p_value) < p_threshold else 0.0
        thresholded_map.setdefault((record.p1, record.p3), []).append(
            (record.p2, thresholded_value, float(record.d_z_score))
        )
        raw_map.setdefault((record.p2, record.p3), []).append((record.p1, 0.0, 0.0))
        thresholded_map.setdefault((record.p2, record.p3), []).append((record.p1, 0.0, 0.0))

    rows: list[dict[str, Any]] = []
    branch_counter = 0

    for node in tree.treenode.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        if parent is None:
            continue
        sibling = next(child for child in parent.children if child is not node)
        desc = sorted(_decode_tree_tips(node))
        desc_set = set(desc)
        sister = sorted(_decode_tree_tips(sibling))
        branch_counter += 1
        for recipient in tips:
            raw_bmins: list[float] = []
            thresholded_bmins: list[float] = []
            z_bmins: list[float] = []
            for sister_taxon in sister:
                raw_entries = raw_map.get((sister_taxon, recipient), [])
                thresholded_entries = thresholded_map.get((sister_taxon, recipient), [])
                raw_vals = [value for b_taxon, value, _z in raw_entries if b_taxon in desc_set]
                thresholded_vals = [
                    value for b_taxon, value, _z in thresholded_entries if b_taxon in desc_set
                ]
                z_vals = [z_val for b_taxon, _value, z_val in thresholded_entries if b_taxon in desc_set]
                if raw_vals:
                    raw_bmins.append(float(np.min(raw_vals)))
                if thresholded_vals:
                    thresholded_bmins.append(float(np.min(thresholded_vals)))
                if z_vals:
                    z_bmins.append(float(np.min(z_vals)))
            raw_value = float(np.median(raw_bmins)) if raw_bmins else float("nan")
            thresholded_value = (
                float(np.median(thresholded_bmins)) if thresholded_bmins else float("nan")
            )
            z_value = float(np.median(z_bmins)) if z_bmins else float("nan")
            rows.append(
                {
                    "branch_id": f"b{branch_counter}",
                    "branch_descendants": ",".join(desc),
                    "sister_descendants": ",".join(sister),
                    "recipient": recipient,
                    "n_sister_taxa_used": len(raw_bmins),
                    "pthresh": p_threshold,
                    "f_branch_raw": raw_value,
                    "f_branch": thresholded_value,
                    "z_branch": z_value,
                    "p_branch": _normal_two_tailed_pvalue(z_value),
                }
            )

    matrix = _build_f_branch_matrix(rows, tips, "f_branch")
    z_matrix = _build_f_branch_matrix(rows, tips, "z_branch")
    p_matrix = _build_f_branch_matrix(rows, tips, "p_branch")
    return pd.DataFrame(rows), matrix, z_matrix, p_matrix


def _build_manifest(
    request: BabaRequest,
    extracter: SNPsExtracter,
    quartets: list[QuartetDefinition],
    tree_meta: dict[str, Any],
    result: BabaResult,
    output_paths: dict[str, Path],
) -> str:
    """Render one human-readable manifest."""
    data_summary = _build_data_summary(request, extracter, quartets, tree_meta, result)
    results_summary = _build_results_summary(result)
    lines = [
    ]
    if request.logged_command:
        lines.extend([f"CMD: {request.logged_command}", ""])
    lines.extend(
        [
            "Inputs",
            "------",
            "tool: baba",
            f"data: {request.data}",
            f"tests: {request.tests}",
            f"tree: {request.tree}",
            f"imap: {request.imap}",
            f"minmap: {request.minmap}",
            f"exclude: {list(request.exclude)}",
            f"include_reference: {request.include_reference}",
            "",
            "Run",
            "---",
            f"name: {request.name}",
            f"outdir: {request.outdir}",
            f"namespace: {tree_meta['namespace']}",
            f"input_mode: {tree_meta['input_mode']}",
            f"min_sample_coverage: {request.min_sample_coverage}",
            f"resampling: {request.resampling}",
            f"bootstrap_replicates: {request.bootstrap_replicates}",
            f"jackknife_block_bp: {request.jackknife_block_bp}",
            f"jackknife_block_loci: {request.jackknife_block_loci}",
            f"seed: {request.seed}",
            f"f_branch: {request.f_branch}",
            f"f_branch_p_threshold: {request.f_branch_p_threshold}",
            f"write_block_table: {request.write_block_table}",
            f"clustering_stats: {request.clustering_stats}",
            f"cores: {request.cores}",
            f"log_level: {request.log_level}",
            "",
            "Data Summary",
            "------------",
        ]
    )
    for key, value in data_summary.items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "",
            "Results Summary",
            "---------------",
        ]
    )
    for key, value in results_summary.items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "",
            "Quartet Preview",
            "--------------",
            _render_table_preview(
                result.quartets,
                ["p1", "p2", "p3", "p4", "n_sites_tested", "abba", "baba", "bbaa", "d_stat", "p_value"],
            ).rstrip(),
            "",
            "Rooted Preview",
            "-------------",
            _render_table_preview(
                result.rooted,
                ["orientation", "p1", "p2", "p3", "p4", "d_stat", "f_g", "d_p_value"],
            ).rstrip(),
        ]
    )
    if result.f_branch is not None:
        lines.extend(
            [
                "",
                "F-branch Preview",
                "----------------",
                _render_table_preview(
                    result.f_branch,
                    ["branch_id", "branch_descendants", "recipient", "f_branch_raw", "f_branch", "z_branch", "p_branch"],
                ).rstrip(),
            ]
        )
    if result.blocks is not None:
        lines.extend(
            [
                "",
                "Block Preview",
                "-------------",
                _render_table_preview(
                    result.blocks,
                    ["quartet_id", "block_label", "n_sites", "abba", "baba", "bbaa", "d_stat"],
                ).rstrip(),
            ]
        )
    lines.extend(
        [
            "",
            "Outputs",
            "-------",
        ]
    )
    for key, path in output_paths.items():
        lines.append(f"{key}: {path}")
    if request.f_branch:
        lines.extend(
            [
                "",
                "Notes",
                "-----",
                "f_branch uses tree-orientation f_G values and zeros non-significant rows before thresholded branch aggregation.",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _run_baba(request: BabaRequest) -> BabaResult:
    """Execute one baba run and return typed outputs."""
    if request.minmap is not None and request.imap is None:
        raise IPyradError("--minmap requires --imap.")
    if request.f_branch and request.tree is None:
        raise IPyradError("--f-branch requires --tree.")
    if request.f_branch_p_threshold <= 0 or request.f_branch_p_threshold > 1:
        raise IPyradError("--f-branch-p-threshold must be > 0 and <= 1.")
    logger.info(
        "starting baba run with {} input, resampling={}, and f_branch={}",
        "manual quartets" if request.tests is not None else "guide tree",
        request.resampling,
        request.f_branch,
    )

    quartets, run_imap, run_exclude, quartet_minmap, tree_text, tree_meta = (
        _resolve_quartets_and_selection(request)
    )
    logger.info(
        "extracting SNP data for {} resolved quartet labels in {} namespace",
        tree_meta["used_label_count"],
        tree_meta["namespace"],
    )

    extracter = SNPsExtracter(
        data=request.data,
        min_sample_coverage=request.min_sample_coverage,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=run_imap,
        minmap=None,
        exclude=run_exclude,
        include_reference=request.include_reference,
        cores=request.cores,
    )
    logger.info("running SNP extraction and site filtering from {}", request.data)
    extracter.run(log_level=request.log_level)
    logger.info(
        "retained {} samples, {} SNPs, {} loci, and {} scaffolds after filtering",
        len(extracter.snames),
        int(extracter.stats["post_filter_snps"]),
        int(extracter.stats["post_filter_snp_containing_linkage_blocks"]),
        int(np.unique(extracter.snpsmap[:, 3]).size) if extracter.snpsmap.size else 0,
    )

    if run_imap is not None:
        label_to_samples = {
            label: [sample for sample in members if sample in extracter.snames]
            for label, members in run_imap.items()
        }
    else:
        label_to_samples = {name: [name] for name in extracter.snames}
    logger.info(
        "precomputing allele-frequency arrays for {} {} labels",
        len(label_to_samples),
        tree_meta["namespace"],
    )

    freqs, coverage, split1, split2 = _compute_group_arrays(
        extracter.genos,
        extracter.snames,
        label_to_samples,
        request.seed,
    )

    rng = np.random.default_rng(request.seed)
    quartet_rows: list[dict[str, Any]] = []
    rooted_rows: list[dict[str, Any]] = []
    block_rows: list[dict[str, Any]] = []
    total_quartets = len(quartets)
    progress_every = max(1, total_quartets // 10)
    logger.info("processing {} quartets", total_quartets)
    for quartet_idx, quartet in enumerate(quartets, start=1):
        row, rooted, blocks = _summarize_quartet(
            quartet_idx,
            quartet,
            freqs,
            coverage,
            split1,
            split2,
            extracter.snpsmap,
            request,
            quartet_minmap,
            rng,
        )
        quartet_rows.append(row)
        rooted_rows.extend(rooted)
        block_rows.extend(blocks)
        if (
            quartet_idx == 1
            or quartet_idx == total_quartets
            or quartet_idx % progress_every == 0
        ):
            logger.info(
                "processed quartet {}/{} {} {} {} {}: sites={} loci={} resampling={}/{}",
                quartet_idx,
                total_quartets,
                quartet.p1,
                quartet.p2,
                quartet.p3,
                quartet.p4,
                row["n_sites_tested"],
                row["n_loci_tested"],
                row["resampling_mode"],
                row["resampling_unit"],
            )

    quartets_df = pd.DataFrame(quartet_rows)
    rooted_df = pd.DataFrame(rooted_rows)
    resolved_df = pd.DataFrame([quartet.as_dict() for quartet in quartets])
    blocks_df = pd.DataFrame(block_rows) if block_rows else None
    f_branch_df = None
    f_branch_matrix = None
    f_branch_z = None
    f_branch_p = None
    if request.f_branch and tree_text is not None:
        logger.info("aggregating tree-oriented f_branch summaries")
        f_branch_df, f_branch_matrix, f_branch_z, f_branch_p = _build_f_branch_outputs(
            tree_text,
            rooted_df,
            request.f_branch_p_threshold,
        )

    result = BabaResult(
        quartets=quartets_df,
        rooted=rooted_df,
        resolved_tests=resolved_df,
        manifest="",
        summary_json="",
        f_branch=f_branch_df,
        f_branch_matrix=f_branch_matrix,
        f_branch_z=f_branch_z,
        f_branch_p=f_branch_p,
        blocks=blocks_df,
        tree_text=tree_text,
    )
    output_paths = _build_output_paths(result, request)
    result.manifest = _build_manifest(
        request,
        extracter,
        quartets,
        tree_meta,
        result,
        output_paths,
    )
    result.summary_json = _build_summary_json(
        request,
        extracter,
        quartets,
        tree_meta,
        result,
        output_paths,
    )
    return result


def run_baba_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    tests: Path | None,
    tree: Path | None,
    imap: Path | None,
    minmap: Path | None,
    min_sample_coverage: int,
    exclude,
    include_reference: bool,
    resampling: str,
    bootstrap_replicates: int,
    jackknife_block_bp: int,
    jackknife_block_loci: int,
    seed: int | None,
    f_branch: bool,
    f_branch_p_threshold: float,
    write_block_table: bool,
    clustering_stats: bool,
    cores: int,
    force: bool,
    log_level: str,
    logged_command: str | None = None,
) -> None:
    """Run the ABBA/BABA analysis workflow."""
    request = BabaRequest(
        data=require_hdf5_input(data, "baba").expanduser().absolute(),
        name=str(name),
        outdir=Path(outdir).expanduser().absolute(),
        tests=None if tests is None else Path(tests).expanduser().absolute(),
        tree=None if tree is None else Path(tree).expanduser().absolute(),
        imap=None if imap is None else Path(imap).expanduser().absolute(),
        minmap=None if minmap is None else Path(minmap).expanduser().absolute(),
        min_sample_coverage=_coerce_positive_int(min_sample_coverage, "min_sample_coverage"),
        exclude=tuple() if exclude is None else tuple(exclude),
        include_reference=bool(include_reference),
        resampling=resampling,
        bootstrap_replicates=_coerce_positive_int(
            bootstrap_replicates,
            "bootstrap_replicates",
        ),
        jackknife_block_bp=_coerce_positive_int(
            jackknife_block_bp,
            "jackknife_block_bp",
        ),
        jackknife_block_loci=_coerce_positive_int(
            jackknife_block_loci,
            "jackknife_block_loci",
        ),
        seed=None if seed is None else int(seed),
        f_branch=bool(f_branch),
        f_branch_p_threshold=float(f_branch_p_threshold),
        write_block_table=bool(write_block_table),
        clustering_stats=bool(clustering_stats),
        cores=_coerce_positive_int(cores, "cores"),
        force=bool(force),
        log_level=str(log_level),
        logged_command=logged_command,
    )

    request.outdir.mkdir(parents=True, exist_ok=True)
    result = _run_baba(request)
    output_paths = _build_output_paths(result, request)
    ensure_output_paths(output_paths.values(), request.force)

    logger.info("writing baba outputs to {}", request.outdir)
    _write_table(output_paths["quartets"], result.quartets)
    _write_table(output_paths["rooted"], result.rooted)
    _write_table(output_paths["resolved_tests"], result.resolved_tests)
    output_paths["manifest"].write_text(result.manifest, encoding="utf-8")
    output_paths["summary_json"].write_text(result.summary_json, encoding="utf-8")
    if "tree" in output_paths and result.tree_text is not None:
        output_paths["tree"].write_text(result.tree_text, encoding="utf-8")
    if "f_branch" in output_paths and result.f_branch is not None:
        _write_table(output_paths["f_branch"], result.f_branch)
    if "f_branch_matrix" in output_paths and result.f_branch_matrix is not None:
        _write_table(output_paths["f_branch_matrix"], result.f_branch_matrix)
    if "f_branch_z" in output_paths and result.f_branch_z is not None:
        _write_table(output_paths["f_branch_z"], result.f_branch_z)
    if "f_branch_p" in output_paths and result.f_branch_p is not None:
        _write_table(output_paths["f_branch_p"], result.f_branch_p)
    if "blocks" in output_paths and result.blocks is not None:
        _write_table(output_paths["blocks"], result.blocks)

    logger.info(
        "completed baba run: quartets={} rooted_rows={} output_dir={}",
        len(result.quartets.index),
        len(result.rooted.index),
        request.outdir,
    )
    for key, path in output_paths.items():
        logger.info("wrote baba {} to {}", key, path)
