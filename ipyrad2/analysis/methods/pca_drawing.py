#!/usr/bin/env python

"""Optional toyplot rendering helpers for PCA-family outputs."""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ...utils.exceptions import IPyradError
from .common import require_toyplot

if TYPE_CHECKING:
    from .pca import PCAFamilyResult


def ensure_pca_plotting_available() -> None:
    """Fail fast if the optional PCA plotting dependency is unavailable."""
    require_toyplot()


def _require_pca_plot_axes(result: "PCAFamilyResult") -> None:
    """Require a PCA result with at least two coordinate axes."""
    if result.method != "pca":
        raise IPyradError("PCA plotting currently supports only PCA results.")
    first_rep = min(result.coords_by_replicate)
    coords = result.coords_by_replicate[first_rep]
    if coords.shape[1] < 2:
        raise IPyradError("PCA plotting requires at least two PCA axes.")


def _align_replicate_coords(result: "PCAFamilyResult") -> dict[int, np.ndarray]:
    """Return the first two PCA axes for each replicate with signs aligned."""
    reps = sorted(result.coords_by_replicate)
    aligned = {}
    base = result.coords_by_replicate[reps[0]][:, :2].astype(np.float64, copy=True)
    aligned[reps[0]] = base

    for rep in reps[1:]:
        current = result.coords_by_replicate[rep][:, :2].astype(np.float64, copy=True)
        for axis in (0, 1):
            if np.dot(base[:, axis], current[:, axis]) < 0:
                current[:, axis] *= -1.0
        aligned[rep] = current
    return aligned


def _mean_variances(result: "PCAFamilyResult") -> np.ndarray:
    """Return the mean explained-variance ratio across PCA replicates."""
    values = np.array(
        [result.variance_by_replicate[rep] for rep in sorted(result.variance_by_replicate)],
        dtype=np.float64,
    )
    return values.mean(axis=0)


def _sample_to_group(result: "PCAFamilyResult") -> dict[str, str]:
    """Return the plotting group for each sample."""
    mapping = {}
    for group, names in result.extracter.imap.items():
        for name in names:
            mapping[name] = group
    missing = [name for name in result.samples if name not in mapping]
    if missing:
        raise IPyradError(
            "PCA plotting requires every retained sample to belong to one group."
        )
    return mapping


def _build_marker_styles(
    toyplot,
    *,
    groups: list[str],
    nreplicates: int,
    size: int = 10,
) -> tuple[dict[str, object], dict[str, object], list[tuple[str, object]]]:
    """Return centroid and replicate marker styles plus legend items."""
    cycle = max(1, min(8, len(groups)))
    colors = itertools.cycle(
        toyplot.color.broadcast(
            toyplot.color.brewer.map("Spectral"),
            shape=cycle,
        )
    )
    shapes = itertools.cycle(
        np.concatenate(
            [
                np.tile("o", cycle),
                np.tile("s", cycle),
                np.tile("^", cycle),
                np.tile("d", cycle),
                np.tile("v", cycle),
                np.tile("<", cycle),
                np.tile("x", cycle),
            ]
        )
    )

    centroid_styles = {}
    replicate_styles = {}
    legend_items = []
    for group in groups:
        color = next(colors)
        shape = next(shapes)
        try:
            css_color = toyplot.color.to_css(color)
        except Exception:
            css_color = color
        centroid_styles[group] = toyplot.marker.create(
            size=size,
            shape=shape,
            mstyle={
                "fill": css_color,
                "stroke": "#262626",
                "stroke-opacity": 1.0,
                "stroke-width": 1.5,
                "fill-opacity": 0.75,
            },
        )
        replicate_styles[group] = toyplot.marker.create(
            size=size,
            shape=shape,
            mstyle={
                "fill": css_color,
                "stroke": "none",
                "fill-opacity": 0.9 / max(nreplicates, 1),
            },
        )
        legend_items.append((group, centroid_styles[group]))
    return centroid_styles, replicate_styles, legend_items


def _stroke_width(style: dict[str, object]) -> float:
    """Return a numeric stroke width from a toyplot style dict."""
    value = style.get("stroke-width", 0.0)
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _set_axes_ticks_external(axes, *, show_domain: bool = True) -> None:
    """Apply external tick styling to the PCA axes."""
    axes.x.ticks.show = True
    axes.y.ticks.show = True
    axes.x.domain.show = show_domain
    axes.y.domain.show = show_domain
    axes.x.ticks.near = 5
    axes.x.ticks.far = 0
    axes.y.ticks.near = 5
    axes.y.ticks.far = 0
    axes.x.ticks.labels.offset = 10
    axes.y.ticks.labels.offset = 10
    axes.x.label.offset = 30
    axes.y.label.offset = 30
    axes.label.offset = 20

    # Equivalent to increasing default margin by 10 px on each side.
    axes._xmin_range += 10
    axes._xmax_range -= 10
    axes._ymin_range += 10
    axes._ymax_range -= 10


def _add_axes_box_outline(canvas, axes) -> object:
    """Draw a box outline around the PCA plotting axes and place it behind."""
    style = {
        "stroke": "#262626",
        "stroke-width": 2.0,
        "fill": "none",
    }
    inset = _stroke_width(style) / 2.0
    pad = float(axes.padding)
    left = float(axes._xmin_range) - pad + inset
    right = float(axes._xmax_range) + pad - inset
    top = float(axes._ymin_range) - pad + inset
    bottom = float(axes._ymax_range) + pad - inset
    if not (left < right and top < bottom):
        raise IPyradError("PCA plot outline region collapsed.")

    overlay = canvas.cartesian(
        margin=0,
        padding=0,
        show=False,
        xshow=False,
        yshow=False,
        xmin=0,
        xmax=canvas.width,
        ymin=canvas.height,
        ymax=0,
    )
    overlay.rectangle(left, right, top, bottom, style=style)

    render_targets = canvas._scenegraph._relationships["render"]._targets[canvas]
    if overlay in render_targets and axes in render_targets:
        render_targets.remove(overlay)
        render_targets.insert(render_targets.index(axes), overlay)
    return overlay


def write_pca_svg_plot(
    result: "PCAFamilyResult",
    outfile: Path | str,
    *,
    width: int = 400,
    height: int = 300,
    marker_size: int = 10,
) -> None:
    """Write a default SVG PCA plot using the first two principal components."""
    _require_pca_plot_axes(result)
    toyplot, toyplot_svg = require_toyplot()

    aligned = _align_replicate_coords(result)
    variances = _mean_variances(result)
    sample_to_group = _sample_to_group(result)
    groups = list(result.extracter.imap)
    nreplicates = len(aligned)

    centroid_styles, replicate_styles, legend_items = _build_marker_styles(
        toyplot,
        groups=groups,
        nreplicates=nreplicates,
        size=marker_size,
    )
    centroid_markers = [centroid_styles[sample_to_group[name]] for name in result.samples]
    replicate_markers = [replicate_styles[sample_to_group[name]] for name in result.samples]

    xlab = f"PC1 ({variances[0] * 100:.1f}% explained)"
    ylab = f"PC2 ({variances[1] * 100:.1f}% explained)"

    legend_width = min(140, max(90, width // 4))
    canvas = toyplot.Canvas(width, height)
    axes = canvas.cartesian(
        xlabel=xlab,
        ylabel=ylab,
        bounds=(60, -(60 + legend_width), 60, -60),
        padding=20,
    )
    axes.x.spine.style["stroke-width"] = 1.5
    axes.y.spine.style["stroke-width"] = 1.5
    axes.x.ticks.labels.style["font-size"] = "12px"
    axes.y.ticks.labels.style["font-size"] = "12px"
    axes.x.label.style["font-size"] = "14px"
    axes.y.label.style["font-size"] = "14px"
    axes.x.ticks.locator = toyplot.locator.Extended(only_inside=True)
    axes.y.ticks.locator = toyplot.locator.Extended(only_inside=True)
    _set_axes_ticks_external(axes)
    _add_axes_box_outline(canvas, axes)

    reps = sorted(aligned)
    if nreplicates == 1:
        coords = aligned[reps[0]]
        axes.scatterplot(
            coords[:, 0],
            coords[:, 1],
            marker=centroid_markers,
            title=result.samples,
        )
    else:
        for rep in reps:
            coords = aligned[rep]
            axes.scatterplot(
                coords[:, 0],
                coords[:, 1],
                marker=replicate_markers,
            )
        centroids = np.mean(
            np.stack([aligned[rep] for rep in reps], axis=0),
            axis=0,
        )
        axes.scatterplot(
            centroids[:, 0],
            centroids[:, 1],
            marker=centroid_markers,
            title=result.samples,
        )

    canvas.legend(
        legend_items,
        bounds=(-legend_width - 60, -60, 60, -60),
    )
    toyplot_svg.render(canvas, str(outfile))
