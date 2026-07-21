

# window extracter concatenation + tree inference

The window extracter (wex) tool in ipyrad2 makes it easy to generate a concatenated alignment
from all or a filtered subset of assembled loci in a dataset. This tutorial demonstrates
this workflow using an empirical dataset from the [SE denovo tutorial]().

## Example

### run wex
See the [wex]() documentation page for a more detailed description of the options for this tool, and
the [wex recipes]() page for more examples of what this tool can implement.

Here let's filter this denovo assembled dataset to keep all sites that have data for >= 8 samples,
and let's exclude any sample that has missing data for >90% of sites.

```bash
ipyrad2 wex \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-wex \
  -n assembly_min8 \
  -m 8 \
  -r 0.90
```

The logging report shows:
```literal
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | -------------------------------------------------------
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | ----- ipyrad2 wex: extract alignments from windows -----
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | -------------------------------------------------------
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | CMD: ipyrad2 wex -d SRP021469/OUT/assembly.hdf5 -m 8 -o SRP021469/output-wex -n assembly_min8 -r 0.9
2026-07-20 16:59:25 | INFO     | window_extracter.py  | No windows specified; selecting the full length of all scaffolds. Use -w to subset scaffold windows and -P to view scaffold names.
2026-07-20 16:59:25 | INFO     | window_extracter.py  | selected 45448 windows from 45448 scaffolds
2026-07-20 16:59:43 | INFO     | window_extracter.py  | wrote alignment (13, 1965012) to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.phy
2026-07-20 16:59:43 | INFO     | window_extracter.py  | wrote stats/log to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.stats.txt
```

### wex stats

This generated an alignment that is 13 taxa x 1.96M sites. The stats file shows additional information.

```bash
cat SRP021469/output-wex/assembly_min8.stats.txt
```

```literal
CMD: ipyrad2 wex -d SRP021469/OUT/assembly.hdf5 -m 8 -o SRP021469/output-wex -n assembly_min8 -r 0.9

# Extract Summary
infile                    SRP021469/OUT/assembly.hdf5
outfile                   /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.phy
out_format                phy
windows_selected          45,448
selected_windows_preview  locus_1_1:1-69, locus_2_1:1-69, locus_3_1:1-72, locus_4_25:1-71, locus_5_1:1-69, locus_6_1:1-69, locus_7_1:1-69, locus_8_1:1-69, locus_9_1:1-62, locus_10_1:1-69, ... (45448 total)

# Filtering Summary
populations                     all
min_sample_coverage_filter      all=8
max_sample_missing              0.900000
samples_selected_initial        13
samples_dropped_by_max_missing  0
samples_final                   13

# Alignment Summary
nsamples_before_filtering              13
nsites_in_windows_before_filtering     3,100,317
nvariants_in_windows_before_filtering  183,379
nsamples_after_filtering               13
nsites_in_windows_after_filtering      1,965,012
nvariants_in_windows_after_filtering   126,655

# Sample Summary
sample                  population  percent_missing  dropped_by_max_missing
29154_superba           all         26.241           no
30556_thamno            all         6.652            no
30686_cyathophylla      all         14.528           no
32082_przewalskii       all         49.630           no
33413_thamno            all         33.084           no
33588_przewalskii       all         43.547           no
35236_rex               all         4.721            no
35855_rex               all         4.861            no
38362_rex               all         4.950            no
39618_rex               all         13.964           no
40578_rex               all         3.844            no
41478_cyathophylloides  all         6.613            no
41954_cyathophylloides  all         10.684           no
```

Here we did not use the `-w` option to limit the analysis to specific scaffolds. That option is more
appropriate when analyzing reference assembled rather than denovo assembled datasets. So you can see
in the stats file above that it selected all 45,448 scaffolds in the denovo assembly, where each locus
is stored as a distinct scaffold.

Across these scaffolds, wex found 3.1M sites and 183K variant sites before filtering,
which was reduced to 1.96M sites and 126K variant sites after filtering. This filtering
primarily represents the effect of the `-m` filter, which excluded sites that did not have
data for 8 samples or more.

Finally, we can see in the Sample Summary section of the stats file that no samples were dropped
because none exceeded the `-r` maximum missing threshold of 0.90. The percent missing data among
samples varied from 3.8 to 49.6% across the samples in the final filtered alignment.

### run raxml-ng (tree inference)

Here we use `raxml-ng` to infer a phylogenetic tree from the supermatrix alignment generated by `wex`. We specify the `--all` option which will perform a tree search and bootstrap analysis to calculate support values. We specify the input (`--msa`) as the phy file produced in the previous step, indicate the substitution model choice, here using the common default `GTR+G`, and tell it to perform 100 non-parametric bootstrap replicate searches (`--bs-trees`). This will likely take 20 minutes or more to run.

```bash
raxml-ng \
    --all \
    --msa SRP021469/output-wex/assembly_min8.phy \
    --model GTR+G \
    --bs-trees 100 \
    --workers 2
```


### plot tree

When it finishes you an examine the results by plotting a tree. The default result files are saved to the same folder as the input file.
One easy way to plot a tree is using [`toytree`](https://eaton-lab.org/toytree/), which can plot trees either in the terminal, or as high quality vector graphics in formats like PDF. It can also be used to perform operations on a tree, such as re-rooting on an outgroup, as shown below.

Let's start by rooting the tree on an outgroup, which in this case is the taxa labeled 'przewalksii'.

```bash
toytree root \
    -i SRP021469/output-wex/assembly_min8.phy.raxml.support \
    -o SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted \
    -n "~prz" \
    --mad
```

Then we can print a tree visualization to the terminal:

```bash
toytree view \
    -i SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted \
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

Or you can open a python or jupyter session and use toytree interactively in Python to generate
a tree figure with many more styling options. See the toytree docs. [Code below is Python not bash]

```python
import toytree

# path to the NEWICK tree file
NEWICK = "SRP021469/output-wex/assembly_min8.phy.raxml.support"

# parse the newick
tree = toytree.tree(NEWICK)

# root and ladderize the tree
tree = tree.root("~prz").ladderize()

# draw the tree
c, a, m = tree.draw(
    tip_labels_align=True,
    node_labels="support",
    node_labels_style={
        'anchor-shift': '10px',
        'baseline-shift': '10px',
    }
)

# save to PDF
toytree.save(canvas, "./tree-drawing.pdf")
```

![]()
