#!/usr/bin/env python

"""Deprecated result-summary helpers for legacy BPP workflows."""

from __future__ import annotations

import itertools
import warnings

import pandas as pd

from ...utils.exceptions import IPyradError


warnings.warn(
    "ipyrad2.analysis.methods.bpp_summary is deprecated and retained only for "
    "legacy notebook workflows.",
    DeprecationWarning,
    stacklevel=2,
)


try:
    import toytree
except ImportError as exc:  # pragma: no cover - optional dependency
    raise IPyradError(
        "You are missing required packages to use legacy BPP summary helpers.\n"
        "First run: conda install toytree -c conda-forge"
    ) from exc


def parse_a00_output(path) -> pd.DataFrame:
    """Parse one legacy BPP A00 summary file into a dataframe."""
    nodes = []
    rows = []
    in_nodes = False
    in_table = False
    with open(path, encoding="utf-8") as infile:
        for raw in infile:
            line = raw.strip().split()
            if "(+1)" in line:
                in_nodes = True
            elif in_nodes and line == []:
                in_nodes = False
            elif in_nodes:
                nodes.append(line[3])

            if "param" in line and "rho1" in line:
                rows.append(line)
                in_table = True
            elif "lnL" in line:
                rows.append(line)
                in_table = False
            elif in_table:
                rows.append(line)

    if len(rows) < 4 or not nodes:
        raise IPyradError(f"Failed to parse BPP A00 output: {path}")

    rows.pop(1)
    rows.pop(-2)
    nodes = pd.DataFrame(nodes)
    nodes.index = [index + 1 for index in nodes.index]
    params = pd.DataFrame(rows[1:], columns=rows[0]).T
    params.columns = params.iloc[0]
    params = params.iloc[1:]

    def _header(value: str) -> str:
        if value == "lnL":
            return value
        parameter, idx = value.split(":")
        return f"{parameter}_{idx}{nodes.loc[int(idx)][0]}"

    params.columns = [_header(value) for value in params.columns]
    return params.astype(float)


def summarize_a01_files(paths) -> pd.DataFrame:
    """Average posterior species-delimitation tables across legacy runs."""
    tables = []
    for path in paths:
        with open(path, encoding="utf-8") as infile:
            dat = infile.read().split("posterior\n")[1]
            table, _rest = dat.split("Order of ancestral nodes:")
            data = [line.strip().split() for line in table.strip().split("\n")]
            df = pd.DataFrame(
                data=data,
                columns=["x", "delim", "prior", "posterior"],
            )
            df = df.drop(columns=["x"])
            df["nspecies"] = [item.count("1") + 1 for item in df["delim"]]
            df["posterior"] = df["posterior"].astype(float)
            tables.append(df)
    result = tables[0].copy()
    for table in tables[1:]:
        result["posterior"] += table["posterior"]
    result["posterior"] /= len(tables)
    return result


def summarize_a10_files(outfiles, mcmcfiles):
    """Return majority-rule trees and posterior tree samples for legacy A10 runs."""
    trees = []
    treelists = []
    for treefile in outfiles:
        with open(treefile, encoding="utf-8") as infile:
            line = None
            for line in infile:
                pass
        if line is None:
            raise IPyradError(f"BPP tree summary file was empty: {treefile}")
        newick = line.split(";")[0] + ";"
        newick = newick.replace(" #", "")
        tree = toytree.tree(newick)
        for node in tree.treenode.traverse():
            node.support = int(round(node.support * 100))
        trees.append(tree)

    for treefile in mcmcfiles:
        with open(treefile, encoding="utf-8") as infile:
            treelists.append(infile.readlines())

    return trees, toytree.mtree(list(itertools.chain(*treelists)))
