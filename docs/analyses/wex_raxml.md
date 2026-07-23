# Seqex concatenation and tree inference

The Seqex tool can filter assembled loci and concatenate them into a
supermatrix for phylogenetic inference. This tutorial demonstrates that
workflow using the empirical dataset from the
[southeastern *Pedicularis* tutorial](../assembly/tutorial-pedic.md).

## Extract and concatenate loci

The command below retains sites represented by at least eight samples, removes
samples with more than 90% missing data, and uses `-C` to concatenate the
filtered loci into one PHYLIP alignment.

```bash
ipyrad2 seqex \
    -d SRP021469/OUT/assembly.hdf5 \
    -o SRP021469/output-seqex \
    -n assembly_min8 \
    -m 8 \
    -r 0.90 \
    -C
```

Without `-w`, Seqex considers every complete locus in the assembly. This is
usually the desired behavior for a de novo assembly, where each assembled RAD
locus is stored on its own scaffold. A scaffold name supplied with `-w` also
selects complete loci. Coordinate regions, such as `scaffold:start-end` or
regions in a BED file, instead clip overlapping loci to the selected interval
before filtering.

For this example, Seqex writes a 13-sample alignment containing 1,965,012
sites:

```literal
wrote 29,339 filtered loci as one concatenated PHYLIP alignment
stats report: SRP021469/output-seqex/assembly_min8.stats.txt
```

The human-readable report records the selected windows, whether coordinate
clipping was applied, filtering settings, output layout, total bases,
non-missing occupancy, per-sample occupancy, and per-locus statistics.

```bash
less SRP021469/output-seqex/assembly_min8.stats.txt
```

The same information is available for programmatic use in
`assembly_min8.stats.json`.

## Infer a concatenation tree

Run RAxML-NG on the concatenated alignment. The `--all` workflow performs a
tree search and bootstrap analysis; `--bs-trees 100` requests 100
non-parametric bootstrap replicates.

```bash
raxml-ng \
    --all \
    --msa SRP021469/output-seqex/assembly_min8.phy \
    --model GTR+G \
    --bs-trees 100 \
    --workers 2
```

## Root and view the tree

The default RAxML-NG result files are written beside the input alignment. The
following commands use [Toytree](https://eaton-lab.org/toytree/) to root the
support tree on the *przewalskii* samples and display it in the terminal.

```bash
toytree root \
    -i SRP021469/output-seqex/assembly_min8.phy.raxml.support \
    -o SRP021469/output-seqex/assembly_min8.phy.raxml.support.rooted \
    -n "~prz" \
    --mad

toytree view \
    -i SRP021469/output-seqex/assembly_min8.phy.raxml.support.rooted \
    --ladderize
```

```literal
                               ┌───32082_przewalskii
┌──────────────────────────────┤
│                              └───33588_przewalskii
│
│                       ┌─────────29154_superba
│                    ┌──┤
│                    │  └─────────30686_cyathophylla
│                ┌───┤
│                │   │            ┌─41954_cyathophylloides
│                │   └────────────┤
└────────────────┤                └─41478_cyathophylloides
                 │
                 │      ┌───────33413_thamno
                 │      │
                 └──────┤┌─────────30556_thamno
                        ││
                        └┤   ┌────35855_rex
                         │┌──┤
                         ││  └────40578_rex
                         └┤
                          │┌────────35236_rex
                          └┤
                           │      ┌──38362_rex
                           └──────┤
                                  └──39618_rex
```

For a publication-quality PDF, draw the rooted tree directly:

```bash
toytree draw \
    -i SRP021469/output-seqex/assembly_min8.phy.raxml.support.rooted \
    -o SRP021469/output-seqex/assembly_min8.phy.raxml.support.rooted.pdf \
    --node-labels support \
    --ladderize
```

Or customize it interactively in Python:

```python
import toytree

newick = "SRP021469/output-seqex/assembly_min8.phy.raxml.support"
tree = toytree.tree(newick).root("~prz").ladderize()
canvas, axes, mark = tree.draw(
    tip_labels_align=True,
    node_labels="support",
    node_labels_style={
        "anchor-shift": "10px",
        "baseline-shift": "10px",
    },
)
toytree.save(canvas, "./tree-drawing.pdf")
```
