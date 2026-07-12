"""Streamlit app for interactively exploring ipyrad2 assembly outputs.

Run with:

    streamlit run ipyrad2/apps/assembly_browser/app.py
"""

from __future__ import annotations

import argparse
import os
import signal
import threading
from pathlib import Path

import pandas as pd
import streamlit as st

try:
    import plotly.express as px
except ImportError:  # pragma: no cover - UI dependency guard
    px = None

from ipyrad2.apps.assembly_browser.data import AssemblyStore
from ipyrad2.apps.assembly_browser.filters import FilterParams
from ipyrad2.apps.assembly_browser.filters import apply_filters


DEFAULT_PATH = "ip2-pe/ip2-pe_outfiles"


def _parse_app_args() -> argparse.Namespace:
    """Parse arguments passed after `streamlit run app.py --`."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--assembly-dir", default=None)
    args, _unknown = parser.parse_known_args()
    return args


def _shutdown_process(delay_seconds: float = 0.5) -> None:
    """Terminate the Streamlit server after the current rerun renders."""
    pid = os.getpid()

    def terminate() -> None:
        os.kill(pid, signal.SIGTERM)

    timer = threading.Timer(delay_seconds, terminate)
    timer.daemon = True
    timer.start()


@st.cache_resource(show_spinner=False)
def _load_store(path: str) -> AssemblyStore:
    return AssemblyStore.from_path(path)


@st.cache_data(show_spinner=False)
def _metadata(path: str) -> dict:
    return _load_store(path).metadata()


@st.cache_data(show_spinner=True)
def _apply_filters(path: str, params: FilterParams):
    return apply_filters(_load_store(path), params)


def _metric_grid(values: dict[str, object]) -> None:
    columns = st.columns(len(values))
    for col, (label, value) in zip(columns, values.items(), strict=False):
        col.metric(label, value)


def _plot_histogram(df: pd.DataFrame, column: str, title: str, *, nbins: int = 40) -> None:
    if px is None:
        st.info("Install plotly to view charts.")
        return
    if df.empty or column not in df:
        st.info("No data available for this plot.")
        return
    fig = px.histogram(df, x=column, nbins=nbins, title=title)
    fig.update_layout(margin=dict(l=10, r=10, t=45, b=10), bargap=0.04)
    st.plotly_chart(fig, use_container_width=True)


def _plot_bar(df: pd.DataFrame, x: str, y: str, title: str) -> None:
    if px is None:
        st.info("Install plotly to view charts.")
        return
    if df.empty or x not in df or y not in df:
        st.info("No data available for this plot.")
        return
    fig = px.bar(df, x=x, y=y, title=title)
    fig.update_layout(margin=dict(l=10, r=10, t=45, b=10))
    st.plotly_chart(fig, use_container_width=True)


def _format_int(value: int | float) -> str:
    return f"{int(value):,}"


def _baseline_cards(store: AssemblyStore) -> None:
    summary = store.summary
    if not summary:
        return
    _metric_grid(
        {
            "Samples": _format_int(summary.get("samples", 0)),
            "Final loci": _format_int(summary.get("final_loci_written", 0)),
            "SNP sites": _format_int(summary.get("final_snp_sites_written", 0)),
            "Occupancy": f"{summary.get('alignment_matrix_occupancy_fraction', 0):.3f}",
        }
    )


def _sidebar(path_default: str) -> tuple[str, FilterParams]:
    st.sidebar.header("Assembly")
    path = st.sidebar.text_input("Output directory or HDF5", value=path_default)
    metadata = _metadata(path)

    if st.sidebar.button("Shutdown browser", type="secondary"):
        st.sidebar.warning("Shutting down Streamlit...")
        _shutdown_process()

    sample_names = metadata["samples"]
    selected_samples = st.sidebar.multiselect(
        "Samples",
        options=sample_names,
        default=sample_names,
    )

    st.sidebar.header("Filters")
    max_sample_missing = st.sidebar.slider(
        "Max sample missing",
        min_value=0.0,
        max_value=1.0,
        value=1.0,
        step=0.01,
    )
    max_site_missing = st.sidebar.slider(
        "Max site missing",
        min_value=0.0,
        max_value=1.0,
        value=1.0,
        step=0.01,
    )
    min_sample_coverage = st.sidebar.slider(
        "Min sample coverage",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.01,
    )
    min_maf = st.sidebar.slider(
        "Min minor allele frequency",
        min_value=0.0,
        max_value=0.5,
        value=0.0,
        step=0.01,
    )
    min_depth = st.sidebar.number_input(
        "Min genotype depth",
        min_value=0,
        value=0,
        step=1,
        disabled=not metadata["has_sample_dp"],
    )
    min_site_qual = st.sidebar.number_input(
        "Min site quality",
        min_value=0.0,
        value=0.0,
        step=1.0,
        disabled=not metadata["has_site_qual"],
    )
    chunk_size = st.sidebar.number_input(
        "HDF5 chunk size",
        min_value=1_000,
        max_value=500_000,
        value=50_000,
        step=1_000,
    )

    params = FilterParams(
        samples=tuple(selected_samples),
        max_sample_missing=max_sample_missing,
        max_site_missing=max_site_missing,
        min_sample_coverage=min_sample_coverage,
        min_minor_allele_frequency=min_maf,
        min_genotype_depth=int(min_depth),
        min_site_qual=float(min_site_qual),
        chunk_size=int(chunk_size),
    )
    return path, params


def main() -> None:
    st.set_page_config(
        page_title="ipyrad2 Assembly Browser",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    st.title("ipyrad2 Assembly Browser")
    app_args = _parse_app_args()
    path_default = (
        str(Path(app_args.assembly_dir).expanduser().resolve())
        if app_args.assembly_dir
        else str((Path.cwd() / DEFAULT_PATH).resolve())
    )

    try:
        path, params = _sidebar(path_default)
        store = _load_store(path)
        metadata = store.metadata()
        result = _apply_filters(path, params)
    except Exception as exc:  # pragma: no cover - UI error boundary
        st.error(str(exc))
        return

    st.caption(metadata["hdf5"])
    _baseline_cards(store)
    st.divider()

    _metric_grid(
        {
            "Retained samples": f"{result.totals['retained_samples']} / {result.totals['input_samples']}",
            "Retained SNPs": f"{result.totals['retained_snps']:,} / {result.totals['input_snps']:,}",
            "Retained loci": f"{result.totals['retained_loci']:,} / {result.totals['input_loci']:,}",
            "Depth-masked genotypes": f"{result.totals['depth_masked_genotypes']:,}",
        }
    )

    tabs = st.tabs(["Overview", "Filter Impact", "Samples", "Loci", "Sites", "Export"])

    with tabs[0]:
        left, right = st.columns([1, 1])
        with left:
            st.subheader("Assembly Summary")
            summary = pd.DataFrame([store.summary]).T.rename(columns={0: "value"})
            st.dataframe(summary, use_container_width=True)
        with right:
            st.subheader("Current Filter Summary")
            st.dataframe(result.filter_counts, use_container_width=True, hide_index=True)

    with tabs[1]:
        left, right = st.columns([1, 1])
        with left:
            _plot_bar(result.filter_counts, "filter", "sites", "Sites by filter")
        with right:
            retained = result.site_summary[result.site_summary["pass_filter"]]
            _plot_histogram(retained, "maf", "Retained SNP MAF distribution")

    with tabs[2]:
        sample_df = result.sample_summary.sort_values("missing_fraction", ascending=False)
        left, right = st.columns([1, 1])
        with left:
            _plot_bar(sample_df, "sample", "missing_fraction", "Sample missingness")
        with right:
            if result.dropped_samples:
                st.subheader("Dropped Samples")
                st.write(", ".join(result.dropped_samples))
            else:
                st.subheader("Dropped Samples")
                st.write("None")
        st.dataframe(sample_df, use_container_width=True, hide_index=True)

    with tabs[3]:
        locus_df = result.locus_summary
        left, right = st.columns([1, 1])
        with left:
            _plot_histogram(locus_df, "snps", "SNPs per locus")
        with right:
            _plot_histogram(locus_df, "retained_snps", "Retained SNPs per locus")
        st.dataframe(locus_df, use_container_width=True, hide_index=True)

    with tabs[4]:
        site_df = result.site_summary
        left, right = st.columns([1, 1])
        with left:
            _plot_histogram(site_df, "missing_fraction", "Site missingness")
        with right:
            _plot_histogram(site_df, "called_samples", "Called samples per site")
        st.dataframe(site_df, use_container_width=True, hide_index=True)

    with tabs[5]:
        st.download_button(
            "Download sample summary CSV",
            data=result.sample_summary.to_csv(index=False),
            file_name=f"{metadata['prefix']}.sample-summary.csv",
            mime="text/csv",
        )
        st.download_button(
            "Download locus summary CSV",
            data=result.locus_summary.to_csv(index=False),
            file_name=f"{metadata['prefix']}.locus-summary.csv",
            mime="text/csv",
        )
        st.download_button(
            "Download site summary CSV",
            data=result.site_summary.to_csv(index=False),
            file_name=f"{metadata['prefix']}.site-summary.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    main()
