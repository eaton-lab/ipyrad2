# pca

The `ipyrad2 pca` tool makes it easy to extract a filtered unlinked SNP sites from an assembled dataset,
apply imputation to fill missing cells, fit a dimensionality reduction model (e.g., PCA, t-SNE, UMAP),
and visualize the results.

This tutorial provides a detailed explanation of this tools workflow and the effects of parameter options.
For a quick demonstration of the pca tool you may wish to jump straight to the [Empirical example](#empirical-example)
at the bottom of the page.

## When to Use

Use `ipyrad2 pca` when you want:

- principal component coordinates from filtered SNP data
- quick nonlinear embeddings with t-SNE or UMAP from the same filtered input logic
- per-sample missing-data and imputation summaries alongside the numerical output

`pca` is usually the best first pass. `tsne` and `umap` are useful when you want exploratory embeddings, and
to reduce higher dimensional data to 2D plots, but they are more sensitive to parameter choices and should
be interpreted more cautiously.


## Prerequisites

- the ipyrad2 analysis tools installed: `conda install ipyrad2[analysis] -c conda-forge -c bioconda`
- a SNP-capable HDF5 file produced by `assemble` or `analysis vcf-to-hdf5`
- enough retained samples and SNPs after filtering to support the selected method

## Inputs and Filtering

The shared filtering and sample-selection controls are:

- `-m, --min-sample-coverage`: minimum number of samples with data required at a SNP
- `-r, --max-sample-missing`: drop samples whose missing-data fraction exceeds this threshold, then rerun SNP filtering
- `-a, --min-minor-allele-frequency`: remove low-frequency SNPs after coverage filtering
- `-e, --exclude`: exclude one or more named samples
- `-R, --include-reference`: include `assembly_reference_sequence`
- `-i, --imap`: sample-to-population map for filtering and imputation grouping
- `-g, --minmap`: per-population minimum coverage checks applied on top of `-m` when `imap` is used

By default the command subsamples one SNP per RAD locus. Use `--no-subsample` to keep linked SNPs.

`--seed` affects SNP subsampling, imputation, and method initialization. For PCA replicates, ipyrad2 derives deterministic per-replicate seeds from the base seed.

## Methods

`-M, --method` chooses one of three methods:

- `pca`: principal components analysis
- `tsne`: t-SNE embedding
- `umap`: UMAP embedding

Current method-specific controls are:

- `--replicates`: only valid with `-M pca`
- `--plot`: only valid with `-M pca`
- `--plot-format`: save the Toyplot canvas as `pdf`, `png`, `html`, or `svg`; requires `--plot`
- `--perplexity` and `--max-iter`: only used with `-M tsne`
- `--n-neighbors`: only used with `-M umap`

Important current rules:

- only the PCA method supports running and visualizing multiple replicates
- PCA plotting in the command line currently supports only the first two principal components (TODO: expand plotting options)
<!--
For t-SNE:

- `perplexity` must be greater than zero
- `perplexity` must be smaller than the number of retained samples
- `max_iter` must be at least 250

For UMAP:

- `n_neighbors` must be at least 2
-->

## Imputation and Missing Data

All PCA-family methods in ipyrad2 require a fully imputed genotype matrix. Missing genotypes are not allowed to remain in the numerical matrix passed to PCA, t-SNE, or UMAP.

Supported imputation modes are:

- `sample`: sample missing diploid genotypes from allele frequencies within each imputation group
- `zero-fill`: replace missing genotypes with homozygous reference calls
- TODO: add kmeans clustering method

`sample` is the default and generally the preferred choice. If you provide `imap`, sample-mode imputation uses those groups. If you do not provide `imap`, all retained samples are treated as one imputation group.

`zero-fill` is usually a poor default for exploratory structure analyses because it can pull missing-heavy samples toward the reference state.

Users should inspect `sample_data_summary.tsv` after every run and consider dropping samples that require too much imputation before trusting the structure. Heavy imputation is often a warning that missingness, rather than biology, may be shaping the ordination.

## How ipyrad2 PCA-Family Works

The active implementation is:

1. filter the SNP HDF5 to extract SNPs that pass filtering thresholds
2. optionally subsample a single SNP per locus (i.e., unlinked SNP view)
3. impute the genotype matrix
4. run the chosen numerical method on the imputed matrix
5. write coordinates and summary tables
6. optionally generate a plot

Method-specific details:

- `pca`: uses a direct SVD-based PCA on the centered matrix and writes explained-variance ratios
- `tsne`: uses scikit-learn t-SNE with `init="pca"`
- `umap`: uses `umap-learn` UMAP with `init="spectral"`

PCA replicates are mainly useful when the selected SNP view can change across replicates, especially under
one-SNP-per-locus subsampling and stochastic imputation. PCA itself is otherwise deterministic for a fixed matrix.
Plotting performed on PCA replicates shows low opacity results for all replicates, and a high opacity points for
the centroid of each sample across replicates.

## Outputs

Every run writes:

- `<name>.coords.tsv`
- `<name>.sample_data_summary.tsv`
- `<name>.stats.txt`

PCA runs additionally write:

- `<name>.variance.tsv`

If you add `--plot`, PCA also writes:

- `<name>.plot.<format>` (SVG by default)

### `coords.tsv`

This table contains one row per sample per replicate and can be used to generate a plot.

### `sample_data_summary.tsv`

This table records per-sample missingness and imputation information. Examining it can help you decide whether
some samples should be excluded, and whether to try different parameter options.

For PCA with multiple replicates, the numeric columns in this file are averaged across replicates. The imputation algorithm column stays constant.

### `variance.tsv`

This file is written only for `-M pca`. It contains `replicate`, `axis`, `explained_variance_ratio`. t-SNE and UMAP do not write a variance file.

### `plot.<format>`

This file is written only when you use `--plot` with `-M pca`. Use
`--plot-format` to select `pdf`, `png`, `html`, or `svg`; SVG is the default.

- it is a scatter plot of `PC1` versus `PC2` in the selected format
- the plotting axes are drawn with an external-tick style and a boxed outline
- points are colored by `imap` group, or by the default `all` group if no `imap` is provided
- single-replicate runs show one point per sample
- multi-replicate runs show translucent replicate clouds plus one centroid per sample

This plot is meant as a quick first look, not a replacement for checking the coordinate tables directly.
[Coming soon: API instructions for more detailed plotting.]

### `stats.txt`

This is a human-readable summary file. It includes all of the information necessary to repeat the analysis:

- tool name and input file
- initial, dropped, and retained samples
- `imap`, `minmap`, and reference-inclusion settings
- whether SNPs were subsampled
- random seed
- imputation method
- selected method and replicate count
- exported SNP and linkage-block counts
- number of axes written
- imputation algorithm plus imputed SNP and genotype fractions
- `perplexity`, `max_iter`, or `n_neighbors` when relevant
- shared SNP-extracter filter statistics

## Method Interpretation Notes

### PCA

PCA is usually the most stable and interpretable first pass.

- axis directions are based on variance in the imputed genotype matrix
- `variance.tsv` helps you judge how much structure each principal component explains
- multiple PCA replicates are mainly useful for sensitivity checks, especially when SNP subsampling changes the input matrix

### t-SNE

t-SNE is an exploratory embedding, not a variance-partitioning method.

- distances and cluster spacing are not directly comparable to PCA axes
- the result can change noticeably with `perplexity`, sample count, and imputation
- there is no explained-variance table for t-SNE

### UMAP

UMAP is also an exploratory embedding.

- local and global structure can shift with `n_neighbors`
- the embedding may emphasize visual separation that should not be over-interpreted biologically
- there is no explained-variance table for UMAP

For both t-SNE and UMAP, heavy missing-data imputation deserves extra caution because nonlinear embeddings can amplify preprocessing effects.

## Command Patterns

Basic PCA run:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT
```

Run PCA and write the default SVG plot:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT \
  --plot
```

Save the PCA plot as PNG:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT \
  --plot \
  --plot-format png
```

Run PCA replicates on the default unlinked SNP view:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT \
  --replicates 3 \
  --seed 7
```

Run t-SNE:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o TSNE_OUT \
  -M tsne \
  --perplexity 8 \
  --max-iter 1000
```

Run UMAP:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o UMAP_OUT \
  -M umap \
  --n-neighbors 10
```

Keep linked SNPs and explicitly request zero-fill imputation:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT \
  --no-subsample \
  --impute-method zero-fill
```

Use population-aware filtering and sample-mode imputation groups:

```bash
ipyrad2 analysis pca \
  -d snps.hdf5 \
  -o PCA_OUT \
  -i pops.tsv \
  -g minmap.tsv
```


## Common Failures

- Wrong input type: `analysis pca` expects SNP-capable HDF5, not raw VCF.
- Empty or too-small result after filtering: all methods require at least two retained samples and at least one retained SNP.
- Invalid replicate count: PCA replicates must be at least 1, and t-SNE or UMAP only support one run.
- Invalid t-SNE settings: perplexity must be positive and smaller than the number of retained samples, and `max_iter` must be at least 250.
- Invalid UMAP settings: `n_neighbors` must be at least 2.
- Invalid plotting request: `--plot` currently works only with `-M pca`, and non-default `--plot-format` values require `--plot`.
- Existing output files: use `--force` to overwrite an existing result set.
- Missing dependencies: scikit-learn is required for all methods, `umap-learn` is additionally required for UMAP, and `toyplot` is only required when `--plot` is requested.


## Empirical example

This example demonstrates the `ipyrad2 pca` tool applied to the empirical dataset from the [SE denovo tutorial]().
It uses the following IMAP file that maps each sample in the dataset to a population/taxon name.

```literal
32082_przewalskii       przewalskii
33588_przewalskii       przewalskii
29154_superba           superba
30686_cyathophylla      cyathophylla
41954_cyathophylloides  cyathophylloides
41478_cyathophylloides  cyathophylloides
33413_thamno            thamnophila_subsp_cupuliformis
30556_thamno            thamnophila_subsp_thamnophila
35855_rex               rex_subsp_rex
40578_rex               rex_subsp_rex
35236_rex               rex_subsp_rockii
38362_rex               rex_subsp_lipskyana
39618_rex               rex_subsp_lipskyana
```

### pca with min sample filtering

Here let's filter this denovo assembled dataset to keep all variant sites that have data for >= 8 samples,
and let's exclude any sample that has missing data for >90% of sites.

```bash
ipyrad2 pca \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-pca \
  -n assembly_min8 \
  -m 8 \
  -r 0.90 \
  --imap SRP021469/IMAP.txt \
  --plot \
  --plot-format html
```

??? info "logging report"

    ```
    2026-07-20 20:38:02 | INFO     | cli_analysis.py      | -------------------------------------------------
    2026-07-20 20:38:02 | INFO     | cli_analysis.py      | ---- ipyrad2 pca: numerical PCA-family methods ----
    2026-07-20 20:38:02 | INFO     | cli_analysis.py      | -------------------------------------------------
    2026-07-20 20:38:02 | INFO     | cli_analysis.py      | CMD: ipyrad2 pca -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-pca -n assembly_min8 -m 8 -r 0.90
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | SNP extraction summary
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic samples: 13
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic pre_filter_snps: 204914
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic pre_filter_percent_missing: 31.497 (linked genotype cells missing before site filtering)
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic masked_genotypes_by_min_depth: 0
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_indels_present: 0
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_non_biallelic: 6773
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_mincov: 68550
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_minmap: 68550
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_min_site_qual: 0
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_invariant_after_subsampling: 21183
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic filter_by_minor_allele_frequency: 0
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic post_filter_snps: 121676
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic post_filter_snp_containing_linkage_blocks: 28085
    2026-07-20 20:38:11 | INFO     | snps_extracter.py    | filter statistic post_filter_percent_missing: 16.977 (linked post-filter genotype cells missing before optional subsampling)
    2026-07-20 20:38:11 | INFO     | common.py            | pca SNP imputation: algorithm=sample prepared_matrix_scope=subsampled_unlinked (one SNP per linkage block selected for the downstream method) snp_columns_with_missing=23486/28085 (83.6%) missing_genotype_cells=63965/365105 (17.5%)
    2026-07-20 20:38:11 | INFO     | common.py            | pca prepared SNP summary: prepared_matrix_scope=subsampled_unlinked (one SNP per linkage block selected for the downstream method) samples=13 linked_post_filter_snps=121676 linked_post_filter_snp_containing_linkage_blocks=28085 prepared_snps=28085 prepared_snp_containing_linkage_blocks=28085 subsample=True
    2026-07-20 20:38:11 | INFO     | pca.py               | wrote PCA-family coordinates to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_min8.coords.tsv
    2026-07-20 20:38:11 | INFO     | pca.py               | wrote PCA-family sample data summary to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_min8.sample_data_summary.tsv
    2026-07-20 20:38:11 | INFO     | pca.py               | wrote PCA explained variance to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_min8.variance.tsv
    2026-07-20 20:38:11 | INFO     | pca.py               | wrote PCA-family stats to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_min8.stats.txt
    ```

This generates the plot shown below. Note: you can hover over points to see sample names.


<div class="toyplot" id="tb8301c8706c84869b1a6f3a371a071a4"><svg class="toyplot-canvas-Canvas" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" width="400.0px" height="300.0px" viewBox="0 0 400.0 300.0" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;border-color:#292724;border-style:none;border-width:1.0;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" id="t2dd5c07787cb41a881f1514f470ec095"><g class="toyplot-coordinates-Cartesian" id="tfe9cae058dbf46878480ead7ced14811"><clipPath id="t0fe83a881e8d45d6adae037dfbb17e1e"><rect x="0" y="0" width="400.0" height="300.0"></rect></clipPath><g clip-path="url(#t0fe83a881e8d45d6adae037dfbb17e1e)"><g class="toyplot-mark-Range" style="fill:none;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:2.0" id="tda137bb61f9440ada9fdf6633408db7c"><g class="toyplot-Series"><rect class="toyplot-Datum" x="51.0" y="51.000000000000014" width="198.00000000000003" height="198.0" style="fill:none;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:2.0"></rect></g></g></g></g><g class="toyplot-coordinates-Cartesian" id="t8991a2f7bf604d60a222a2a50c9e40ea"><clipPath id="t18c0ad9375ad4e268ffff65dd8bc90f6"><rect x="50.0" y="50.0" width="200.0" height="200.0"></rect></clipPath><g clip-path="url(#t18c0ad9375ad4e268ffff65dd8bc90f6)"><g class="toyplot-mark-Point" id="t11b0f09d2cce417f93877104e4757772"><g class="toyplot-Series"><g style="fill:rgb(88.7%,32.2%,29%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(86.04524411287719, 144.28922306122197)"><title>29154_superba</title><circle r="5.0"></circle></g><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(80.23100481154518, 210.31429979336778)"><title>30556_thamno</title><rect x="-5.0" y="-5.0" width="10" height="10"></rect></g><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(85.3647200369503, 133.00390046292)"><title>30686_cyathophylla</title><circle r="5.0"></circle></g><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(223.9884393063584, 177.01291892401895)"><title>32082_przewalskii</title><circle r="5.0"></circle></g><g style="fill:rgb(62%,0.4%,25.9%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(85.30831180397706, 197.9937800199105)"><title>33413_thamno</title><circle r="5.0"></circle></g><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(223.94911432466793, 177.1728740875959)"><title>33588_przewalskii</title><circle r="5.0"></circle></g><g style="fill:rgb(98.7%,64.6%,36.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(80.23161189904562, 211.5473973752185)"><title>35236_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(78.04540430226768, 217.4845504501303)"><title>35855_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(76.01156069364161, 223.9884393063584)"><title>38362_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(76.12015000846742, 223.42627817917742)"><title>39618_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(78.02894343074975, 216.83009744288134)"><title>40578_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(77.65834296009835, 76.18554512803183)"><title>41478_cyathophylloides</title><circle r="5.0"></circle></g><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(77.82820791480829, 76.0115606936416)"><title>41954_cyathophylloides</title><circle r="5.0"></circle></g></g></g></g><g class="toyplot-coordinates-Axis" id="t0cd2839afa9e4538bfe3524592aee9be" transform="translate(70.0,230.0)translate(0,20.0)"><line x1="6.011560693641615" y1="0" x2="153.9884393063584" y2="0" style="stroke-width:1.5"></line><g><line x1="10.495929129482146" y1="5.0" x2="10.495929129482146" y2="0" style=""></line><line x1="53.93654096366529" y1="5.0" x2="53.93654096366529" y2="0" style=""></line><line x1="97.37715279784842" y1="5.0" x2="97.37715279784842" y2="0" style=""></line><line x1="140.81776463203158" y1="5.0" x2="140.81776463203158" y2="0" style=""></line></g><g><g transform="translate(10.495929129482146,10.0)"><text x="-8.67" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">-25</text></g><g transform="translate(53.93654096366529,10.0)"><text x="-6.672000000000001" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">25</text></g><g transform="translate(97.37715279784842,10.0)"><text x="-6.672000000000001" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">75</text></g><g transform="translate(140.81776463203158,10.0)"><text x="-10.008" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">125</text></g></g><g transform="translate(80.0,30.0)"><text x="-74.30499999999999" y="11.977" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:14.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">PC1 (32.5% explained)</text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line x1="0" x2="0" y1="-5.0" y2="7.5" style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0"></line><text x="0" y="-10.0" style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle"></text></g></g><g class="toyplot-coordinates-Axis" id="t103e5f9c774340a3a72022fc53e7feee" transform="translate(70.0,230.0)rotate(-90.0)translate(0,-20.0)"><line x1="6.011560693641619" y1="0" x2="153.9884393063584" y2="0" style="stroke-width:1.5"></line><g><line x1="3.8122782439135996" y1="-5.0" x2="3.8122782439135996" y2="0" style=""></line><line x1="54.21070269811773" y1="-5.0" x2="54.21070269811773" y2="0" style=""></line><line x1="104.60912715232185" y1="-5.0" x2="104.60912715232185" y2="0" style=""></line><line x1="155.00755160652596" y1="-5.0" x2="155.00755160652596" y2="0" style=""></line></g><g><g transform="translate(3.8122782439135996,-10.0)"><text x="-8.67" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">-50</text></g><g transform="translate(54.21070269811773,-10.0)"><text x="-3.3360000000000003" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">0</text></g><g transform="translate(104.60912715232185,-10.0)"><text x="-6.672000000000001" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">50</text></g><g transform="translate(155.00755160652596,-10.0)"><text x="-10.008" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">100</text></g></g><g transform="translate(80.0,-30.0)"><text x="-74.30499999999999" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:14.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">PC2 (23.1% explained)</text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line x1="0" x2="0" y1="5.0" y2="-7.5" style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0"></line><text x="0" y="10.0" style="alignment-baseline:hanging;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle"></text></g></g></g><g class="toyplot-coordinates-Table" id="t7180f23ca0d24ab08347410f7085ce7d"><g transform="translate(285.0,70.0)"><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,70.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">cyathophylla</text></g><g transform="translate(285.0,90.0)"><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,90.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">cyathophylloides</text></g><g transform="translate(285.0,110.0)"><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,110.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">przewalskii</text></g><g transform="translate(285.0,130.0)"><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,130.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_lipskyana</text></g><g transform="translate(285.0,150.0)"><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,150.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_rex</text></g><g transform="translate(285.0,170.0)"><g style="fill:rgb(98.7%,64.6%,36.4%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,170.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_rockii</text></g><g transform="translate(285.0,190.0)"><g style="fill:rgb(88.7%,32.2%,29%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,190.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">superba</text></g><g transform="translate(285.0,210.0)"><g style="fill:rgb(62%,0.4%,25.9%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,210.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">thamnophila_subsp_cupuliformis</text></g><g transform="translate(285.0,230.0)"><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><rect x="-5.0" y="-5.0" width="10.0" height="10.0"></rect></g></g><g transform="translate(295.0,230.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">thamnophila_subsp_thamnophila</text></g></g></svg><div class="toyplot-behavior"><script>(function()
{
var modules={};
modules["toyplot/tables"] = (function()
    {
        var tables = [];

        var module = {};

        module.set = function(owner, key, names, columns)
        {
            tables.push({owner: owner, key: key, names: names, columns: columns});
        }

        module.get = function(owner, key)
        {
            for(var i = 0; i != tables.length; ++i)
            {
                var table = tables[i];
                if(table.owner != owner)
                    continue;
                if(table.key != key)
                    continue;
                return {names: table.names, columns: table.columns};
            }
        }

        module.get_csv = function(owner, key)
        {
            var table = module.get(owner, key);
            if(table != undefined)
            {
                var csv = "";
                csv += table.names.join(",") + "\n";
                for(var i = 0; i != table.columns[0].length; ++i)
                {
                  for(var j = 0; j != table.columns.length; ++j)
                  {
                    if(j)
                      csv += ",";
                    csv += table.columns[j][i];
                  }
                  csv += "\n";
                }
                return csv;
            }
        }

        return module;
    })();
modules["toyplot/root/id"] = "tb8301c8706c84869b1a6f3a371a071a4";
modules["toyplot/root"] = (function(root_id)
    {
        return document.querySelector("#" + root_id);
    })(modules["toyplot/root/id"]);
modules["toyplot/canvas/id"] = "t2dd5c07787cb41a881f1514f470ec095";
modules["toyplot/canvas"] = (function(canvas_id)
    {
        return document.querySelector("#" + canvas_id);
    })(modules["toyplot/canvas/id"]);
modules["toyplot/menus/context"] = (function(root, canvas)
    {
        var wrapper = document.createElement("div");
        wrapper.innerHTML = "<ul class='toyplot-context-menu' style='background:#eee; border:1px solid #b8b8b8; border-radius:5px; box-shadow: 0px 0px 8px rgba(0%,0%,0%,0.25); margin:0; padding:3px 0; position:fixed; visibility:hidden;'></ul>"
        var menu = wrapper.firstChild;

        root.appendChild(menu);

        var items = [];

        var ignore_mouseup = null;
        function open_menu(e)
        {
            var show_menu = false;
            for(var index=0; index != items.length; ++index)
            {
                var item = items[index];
                if(item.show(e))
                {
                    item.item.style.display = "block";
                    show_menu = true;
                }
                else
                {
                    item.item.style.display = "none";
                }
            }

            if(show_menu)
            {
                ignore_mouseup = true;
                menu.style.left = (e.clientX + 1) + "px";
                menu.style.top = (e.clientY - 5) + "px";
                menu.style.visibility = "visible";
                e.stopPropagation();
                e.preventDefault();
            }
        }

        function close_menu()
        {
            menu.style.visibility = "hidden";
        }

        function contextmenu(e)
        {
            open_menu(e);
        }

        function mousemove(e)
        {
            ignore_mouseup = false;
        }

        function mouseup(e)
        {
            if(ignore_mouseup)
            {
                ignore_mouseup = false;
                return;
            }
            close_menu();
        }

        function keydown(e)
        {
            if(e.key == "Escape" || e.key == "Esc" || e.keyCode == 27)
            {
                close_menu();
            }
        }

        canvas.addEventListener("contextmenu", contextmenu);
        canvas.addEventListener("mousemove", mousemove);
        document.addEventListener("mouseup", mouseup);
        document.addEventListener("keydown", keydown);

        var module = {};
        module.add_item = function(label, show, activate)
        {
            var wrapper = document.createElement("div");
            wrapper.innerHTML = "<li class='toyplot-context-menu-item' style='background:#eee; color:#333; padding:2px 20px; list-style:none; margin:0; text-align:left;'>" + label + "</li>"
            var item = wrapper.firstChild;

            items.push({item: item, show: show});

            function mouseover()
            {
                this.style.background = "steelblue";
                this.style.color = "white";
            }

            function mouseout()
            {
                this.style.background = "#eee";
                this.style.color = "#333";
            }

            function choose_item(e)
            {
                close_menu();
                activate();

                e.stopPropagation();
                e.preventDefault();
            }

            item.addEventListener("mouseover", mouseover);
            item.addEventListener("mouseout", mouseout);
            item.addEventListener("mouseup", choose_item);
            item.addEventListener("contextmenu", choose_item);

            menu.appendChild(item);
        };
        return module;
    })(modules["toyplot/root"],modules["toyplot/canvas"]);
modules["toyplot/io"] = (function()
    {
        var module = {};
        module.save_file = function(mime_type, charset, data, filename)
        {
            var uri = "data:" + mime_type + ";charset=" + charset + "," + data;
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = filename;

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
        };
        return module;
    })();
modules["toyplot.coordinates.Axis"] = (
        function(canvas)
        {
            function sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return sign(segment.domain.min) * Math.pow(base, mix(log(segment.domain.min, base), log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            var axes = {};

            function display_coordinates(e)
            {
                var current = canvas.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            canvas.addEventListener("click", display_coordinates);

            var module = {};
            module.show_coordinates = function(axis_id, projection)
            {
                axes[axis_id] = projection;
            }

            return module;
        })(modules["toyplot/canvas"]);
(function(tables, context_menu, io, owner_id, key, label, names, columns, filename)
        {
            tables.set(owner_id, key, names, columns);

            var owner = document.querySelector("#" + owner_id);
            function show_item(e)
            {
                return owner.contains(e.target);
            }

            function choose_item()
            {
                io.save_file("text/csv", "utf-8", tables.get_csv(owner_id, key), filename + ".csv");
            }

            context_menu.add_item("Save " + label + " as CSV", show_item, choose_item);
        })(modules["toyplot/tables"],modules["toyplot/menus/context"],modules["toyplot/io"],"t11b0f09d2cce417f93877104e4757772","data","point",["x", "y0"],[[-18.612756877622616, -25.304927010407003, -19.396037829700933, 140.15940282402158, -19.460963518580012, 140.11413990065327, -25.30422825470951, -27.820545940476546, -30.161493181723344, -30.036507240871167, -27.83949234894425, -28.266052260284816, -28.070538261358745], [31.251050585205967, -34.25206528316931, 42.44715713071582, -1.213948288451492, -22.028945307809092, -1.3726389274042072, -35.475414619189934, -41.3656315646688, -47.818103965009406, -47.260386999381055, -40.71635233189867, 98.81633528480441, 98.98894428625083]],"toyplot");
(function(axis, axis_id, projection)
        {
            axis.show_coordinates(axis_id, projection);
        })(modules["toyplot.coordinates.Axis"],"t0cd2839afa9e4538bfe3524592aee9be",[{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 147.07868922425496, "min": -37.080779581956726}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 160.0, "min": 0.0}, "scale": "linear"}]);
(function(axis, axis_id, projection)
        {
            axis.show_coordinates(axis_id, projection);
        })(modules["toyplot.coordinates.Axis"],"t103e5f9c774340a3a72022fc53e7feee",[{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 104.95298062145825, "min": -53.78214030021685}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 160.0, "min": 0.0}, "scale": "linear"}]);
})();</script></div></div>


### pca with minmap filtering

To apply a minmap filter use this MINMAP file that assigns a minimum number of samples that must
be present for each population in the IMAP file.

```literal
przewalskii                     1
superba                         1
cyathophylla                    1
cyathophylloides                1
thamnophila_subsp_cupuliformis  1
thamnophila_subsp_thamnophila   1
rex_subsp_rex                   1
rex_subsp_rockii                1
rex_subsp_lipskyana             1
```

Then run `pca` using an imap/minmap filter and with a different output name.

```bash
ipyrad2 pca \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-pca \
  -n assembly_imap \
  --imap SRP021469/IMAP.txt \
  --minmap SRP021469/MINMAP.txt \
  --plot \
  --plot-format html
```

??? info "logging report"

    ```
    2026-07-20 21:02:38 | INFO     | cli_analysis.py      | -------------------------------------------------
    2026-07-20 21:02:38 | INFO     | cli_analysis.py      | ---- ipyrad2 pca: numerical PCA-family methods ----
    2026-07-20 21:02:38 | INFO     | cli_analysis.py      | -------------------------------------------------
    2026-07-20 21:02:38 | INFO     | cli_analysis.py      | CMD: ipyrad2 pca -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-pca -n assembly_minmap -r 0.90 --plot --imap SRP021469/IMAP.txt --minmap SRP021469/MINMAP.txt
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | SNP extraction summary
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic samples: 13
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic pre_filter_snps: 204914
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic pre_filter_percent_missing: 31.497 (linked genotype cells missing before site filtering)
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic masked_genotypes_by_min_depth: 0
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_indels_present: 0
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_non_biallelic: 6773
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_mincov: 1585
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_minmap: 170096
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_min_site_qual: 0
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_invariant_after_subsampling: 21183
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic filter_by_minor_allele_frequency: 0
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic post_filter_snps: 33050
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic post_filter_snp_containing_linkage_blocks: 7612
    2026-07-20 21:02:47 | INFO     | snps_extracter.py    | filter statistic post_filter_percent_missing: 3.404 (linked post-filter genotype cells missing before optional subsampling)
    2026-07-20 21:02:47 | INFO     | common.py            | pca SNP imputation: algorithm=sample prepared_matrix_scope=subsampled_unlinked (one SNP per linkage block selected for the downstream method) snp_columns_with_missing=2823/7612 (37.1%) missing_genotype_cells=3310/98956 (3.3%)
    2026-07-20 21:02:47 | INFO     | common.py            | pca prepared SNP summary: prepared_matrix_scope=subsampled_unlinked (one SNP per linkage block selected for the downstream method) samples=13 linked_post_filter_snps=33050 linked_post_filter_snp_containing_linkage_blocks=7612 prepared_snps=7612 prepared_snp_containing_linkage_blocks=7612 subsample=True
    2026-07-20 21:02:47 | INFO     | pca.py               | wrote PCA-family coordinates to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_minmap.coords.tsv
    2026-07-20 21:02:47 | INFO     | pca.py               | wrote PCA-family sample data summary to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_minmap.sample_data_summary.tsv
    2026-07-20 21:02:47 | INFO     | pca.py               | wrote PCA explained variance to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_minmap.variance.tsv
    2026-07-20 21:02:47 | INFO     | pca.py               | wrote PCA SVG plot to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_minmap.plot.svg
    2026-07-20 21:02:47 | INFO     | pca.py               | wrote PCA-family stats to /home/deren/Documents/tools/ipyrad2/SRP021469/output-pca/assembly_minmap.stats.txt
    ```

You can see from the logging report and/or stats file that this applied a more stringent filter, reducing the dataset to
7612 unlinked SNPs, with only 3.3% of the remaining sites needing to be imputed.


<div class="toyplot" id="tda8a209305a04a3ab7b5bb9e66698db4"><svg class="toyplot-canvas-Canvas" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" width="400.0px" height="300.0px" viewBox="0 0 400.0 300.0" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;border-color:#292724;border-style:none;border-width:1.0;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" id="t59f06f32b8fd49259d5a2ab793d3f00a"><g class="toyplot-coordinates-Cartesian" id="t69fe84f2f0a94220b0943dcbff959e74"><clipPath id="t27e7f1eb7ee446cca5bf33b05e9d5a06"><rect x="0" y="0" width="400.0" height="300.0"></rect></clipPath><g clip-path="url(#t27e7f1eb7ee446cca5bf33b05e9d5a06)"><g class="toyplot-mark-Range" style="fill:none;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:2.0" id="t29490049e7af4edeacca3706fcaa2dfd"><g class="toyplot-Series"><rect class="toyplot-Datum" x="51.0" y="51.000000000000014" width="198.00000000000003" height="198.0" style="fill:none;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:2.0"></rect></g></g></g></g><g class="toyplot-coordinates-Cartesian" id="t845cd85d6ea342a383170b73f702b6e9"><clipPath id="tbfd01d43b61a49efa09af48750679024"><rect x="50.0" y="50.0" width="200.0" height="200.0"></rect></clipPath><g clip-path="url(#tbfd01d43b61a49efa09af48750679024)"><g class="toyplot-mark-Point" id="t1b517aa90bb34042a4caa7369e9a6ca6"><g class="toyplot-Series"><g style="fill:rgb(88.7%,32.2%,29%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(85.6580803660619, 142.88483723412526)"><title>29154_superba</title><circle r="5.0"></circle></g><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(80.1984614937885, 211.13869542620017)"><title>30556_thamno</title><rect x="-5.0" y="-5.0" width="10" height="10"></rect></g><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(85.56696070320726, 132.15997640183542)"><title>30686_cyathophylla</title><circle r="5.0"></circle></g><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(223.9884393063584, 176.70600093802872)"><title>32082_przewalskii</title><circle r="5.0"></circle></g><g style="fill:rgb(62%,0.4%,25.9%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(85.1081452394738, 198.00506546461213)"><title>33413_thamno</title><circle r="5.0"></circle></g><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(223.84845952200138, 176.73560527668556)"><title>33588_przewalskii</title><circle r="5.0"></circle></g><g style="fill:rgb(98.7%,64.6%,36.4%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(80.51366248216601, 210.05344347611694)"><title>35236_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(78.34927176272431, 216.6978384952394)"><title>35855_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(76.01156069364163, 223.9884393063584)"><title>38362_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(76.07734831401767, 223.15114886427767)"><title>39618_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(78.54265604696892, 215.84526625524032)"><title>40578_rex</title><circle r="5.0"></circle></g><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(77.89629011794347, 76.01156069364163)"><title>41478_cyathophylloides</title><circle r="5.0"></circle></g><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;opacity:1.0;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" class="toyplot-Datum" transform="translate(78.12291011484943, 76.06741507897328)"><title>41954_cyathophylloides</title><circle r="5.0"></circle></g></g></g></g><g class="toyplot-coordinates-Axis" id="t71085bbdcfa0444d8c8bc12140533b6d" transform="translate(70.0,230.0)translate(0,20.0)"><line x1="6.011560693641619" y1="0" x2="153.98843930635843" y2="0" style="stroke-width:1.5"></line><g><line x1="10.708238569908499" y1="5.0" x2="10.708238569908499" y2="0" style=""></line><line x1="53.88903007058378" y1="5.0" x2="53.88903007058378" y2="0" style=""></line><line x1="97.06982157125906" y1="5.0" x2="97.06982157125906" y2="0" style=""></line><line x1="140.25061307193434" y1="5.0" x2="140.25061307193434" y2="0" style=""></line></g><g><g transform="translate(10.708238569908499,10.0)"><text x="-8.67" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">-25</text></g><g transform="translate(53.88903007058378,10.0)"><text x="-6.672000000000001" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">25</text></g><g transform="translate(97.06982157125906,10.0)"><text x="-6.672000000000001" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">75</text></g><g transform="translate(140.25061307193434,10.0)"><text x="-10.008" y="10.265999999999998" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">125</text></g></g><g transform="translate(80.0,30.0)"><text x="-74.30499999999999" y="11.977" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:14.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">PC1 (32.7% explained)</text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line x1="0" x2="0" y1="-5.0" y2="7.5" style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0"></line><text x="0" y="-10.0" style="alignment-baseline:alphabetic;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle"></text></g></g><g class="toyplot-coordinates-Axis" id="t69864239977e4889a4791e7d64b9c4bc" transform="translate(70.0,230.0)rotate(-90.0)translate(0,-20.0)"><line x1="6.011560693641615" y1="0" x2="153.98843930635837" y2="0" style="stroke-width:1.5"></line><g><line x1="4.457736475738312" y1="-5.0" x2="4.457736475738312" y2="0" style=""></line><line x1="54.65805439143564" y1="-5.0" x2="54.65805439143564" y2="0" style=""></line><line x1="104.85837230713295" y1="-5.0" x2="104.85837230713295" y2="0" style=""></line><line x1="155.0586902228303" y1="-5.0" x2="155.0586902228303" y2="0" style=""></line></g><g><g transform="translate(4.457736475738312,-10.0)"><text x="-8.67" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">-50</text></g><g transform="translate(54.65805439143564,-10.0)"><text x="-3.3360000000000003" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">0</text></g><g transform="translate(104.85837230713295,-10.0)"><text x="-6.672000000000001" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">50</text></g><g transform="translate(155.0586902228303,-10.0)"><text x="-10.008" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">100</text></g></g><g transform="translate(80.0,-30.0)"><text x="-74.30499999999999" y="0.0" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:14.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">PC2 (23.2% explained)</text></g><g class="toyplot-coordinates-Axis-coordinates" style="visibility:hidden" transform=""><line x1="0" x2="0" y1="5.0" y2="-7.5" style="stroke:rgb(43.9%,50.2%,56.5%);stroke-opacity:1.0;stroke-width:1.0"></line><text x="0" y="10.0" style="alignment-baseline:hanging;fill:rgb(43.9%,50.2%,56.5%);fill-opacity:1.0;font-size:10px;font-weight:normal;stroke:none;text-anchor:middle"></text></g></g></g><g class="toyplot-coordinates-Table" id="t2b5575dec43f4f0f8fdfb3a169269af1"><g transform="translate(285.0,70.0)"><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,70.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">cyathophylla</text></g><g transform="translate(285.0,90.0)"><g style="fill:rgb(28.3%,63.1%,70.1%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,90.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">cyathophylloides</text></g><g transform="translate(285.0,110.0)"><g style="fill:rgb(63.2%,85.2%,64.4%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,110.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">przewalskii</text></g><g transform="translate(285.0,130.0)"><g style="fill:rgb(93%,97.2%,64%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,130.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_lipskyana</text></g><g transform="translate(285.0,150.0)"><g style="fill:rgb(99.7%,91.3%,60.3%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,150.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_rex</text></g><g transform="translate(285.0,170.0)"><g style="fill:rgb(98.7%,64.6%,36.4%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,170.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">rex_subsp_rockii</text></g><g transform="translate(285.0,190.0)"><g style="fill:rgb(88.7%,32.2%,29%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,190.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">superba</text></g><g transform="translate(285.0,210.0)"><g style="fill:rgb(62%,0.4%,25.9%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><circle r="5.0"></circle></g></g><g transform="translate(295.0,210.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">thamnophila_subsp_cupuliformis</text></g><g transform="translate(285.0,230.0)"><g style="fill:rgb(36.9%,31%,63.5%);fill-opacity:0.75;stroke:rgb(14.9%,14.9%,14.9%);stroke-opacity:1.0;stroke-width:1.5" transform="translate(-5.55, -4.440892098500626e-16)"><rect x="-5.0" y="-5.0" width="10.0" height="10.0"></rect></g></g><g transform="translate(295.0,230.0)"><text x="0" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:normal;stroke:none;vertical-align:baseline;white-space:pre">thamnophila_subsp_thamnophila</text></g></g></svg><div class="toyplot-behavior"><script>(function()
{
var modules={};
modules["toyplot/tables"] = (function()
    {
        var tables = [];

        var module = {};

        module.set = function(owner, key, names, columns)
        {
            tables.push({owner: owner, key: key, names: names, columns: columns});
        }

        module.get = function(owner, key)
        {
            for(var i = 0; i != tables.length; ++i)
            {
                var table = tables[i];
                if(table.owner != owner)
                    continue;
                if(table.key != key)
                    continue;
                return {names: table.names, columns: table.columns};
            }
        }

        module.get_csv = function(owner, key)
        {
            var table = module.get(owner, key);
            if(table != undefined)
            {
                var csv = "";
                csv += table.names.join(",") + "\n";
                for(var i = 0; i != table.columns[0].length; ++i)
                {
                  for(var j = 0; j != table.columns.length; ++j)
                  {
                    if(j)
                      csv += ",";
                    csv += table.columns[j][i];
                  }
                  csv += "\n";
                }
                return csv;
            }
        }

        return module;
    })();
modules["toyplot/root/id"] = "tda8a209305a04a3ab7b5bb9e66698db4";
modules["toyplot/root"] = (function(root_id)
    {
        return document.querySelector("#" + root_id);
    })(modules["toyplot/root/id"]);
modules["toyplot/canvas/id"] = "t59f06f32b8fd49259d5a2ab793d3f00a";
modules["toyplot/canvas"] = (function(canvas_id)
    {
        return document.querySelector("#" + canvas_id);
    })(modules["toyplot/canvas/id"]);
modules["toyplot/menus/context"] = (function(root, canvas)
    {
        var wrapper = document.createElement("div");
        wrapper.innerHTML = "<ul class='toyplot-context-menu' style='background:#eee; border:1px solid #b8b8b8; border-radius:5px; box-shadow: 0px 0px 8px rgba(0%,0%,0%,0.25); margin:0; padding:3px 0; position:fixed; visibility:hidden;'></ul>"
        var menu = wrapper.firstChild;

        root.appendChild(menu);

        var items = [];

        var ignore_mouseup = null;
        function open_menu(e)
        {
            var show_menu = false;
            for(var index=0; index != items.length; ++index)
            {
                var item = items[index];
                if(item.show(e))
                {
                    item.item.style.display = "block";
                    show_menu = true;
                }
                else
                {
                    item.item.style.display = "none";
                }
            }

            if(show_menu)
            {
                ignore_mouseup = true;
                menu.style.left = (e.clientX + 1) + "px";
                menu.style.top = (e.clientY - 5) + "px";
                menu.style.visibility = "visible";
                e.stopPropagation();
                e.preventDefault();
            }
        }

        function close_menu()
        {
            menu.style.visibility = "hidden";
        }

        function contextmenu(e)
        {
            open_menu(e);
        }

        function mousemove(e)
        {
            ignore_mouseup = false;
        }

        function mouseup(e)
        {
            if(ignore_mouseup)
            {
                ignore_mouseup = false;
                return;
            }
            close_menu();
        }

        function keydown(e)
        {
            if(e.key == "Escape" || e.key == "Esc" || e.keyCode == 27)
            {
                close_menu();
            }
        }

        canvas.addEventListener("contextmenu", contextmenu);
        canvas.addEventListener("mousemove", mousemove);
        document.addEventListener("mouseup", mouseup);
        document.addEventListener("keydown", keydown);

        var module = {};
        module.add_item = function(label, show, activate)
        {
            var wrapper = document.createElement("div");
            wrapper.innerHTML = "<li class='toyplot-context-menu-item' style='background:#eee; color:#333; padding:2px 20px; list-style:none; margin:0; text-align:left;'>" + label + "</li>"
            var item = wrapper.firstChild;

            items.push({item: item, show: show});

            function mouseover()
            {
                this.style.background = "steelblue";
                this.style.color = "white";
            }

            function mouseout()
            {
                this.style.background = "#eee";
                this.style.color = "#333";
            }

            function choose_item(e)
            {
                close_menu();
                activate();

                e.stopPropagation();
                e.preventDefault();
            }

            item.addEventListener("mouseover", mouseover);
            item.addEventListener("mouseout", mouseout);
            item.addEventListener("mouseup", choose_item);
            item.addEventListener("contextmenu", choose_item);

            menu.appendChild(item);
        };
        return module;
    })(modules["toyplot/root"],modules["toyplot/canvas"]);
modules["toyplot/io"] = (function()
    {
        var module = {};
        module.save_file = function(mime_type, charset, data, filename)
        {
            var uri = "data:" + mime_type + ";charset=" + charset + "," + data;
            uri = encodeURI(uri);

            var link = document.createElement("a");
            if(typeof link.download != "undefined")
            {
              link.href = uri;
              link.style = "visibility:hidden";
              link.download = filename;

              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            }
            else
            {
              window.open(uri);
            }
        };
        return module;
    })();
modules["toyplot.coordinates.Axis"] = (
        function(canvas)
        {
            function sign(x)
            {
                return x < 0 ? -1 : x > 0 ? 1 : 0;
            }

            function mix(a, b, amount)
            {
                return ((1.0 - amount) * a) + (amount * b);
            }

            function log(x, base)
            {
                return Math.log(Math.abs(x)) / Math.log(base);
            }

            function in_range(a, x, b)
            {
                var left = Math.min(a, b);
                var right = Math.max(a, b);
                return left <= x && x <= right;
            }

            function inside(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(in_range(segment.range.min, range, segment.range.max))
                        return true;
                }
                return false;
            }

            function to_domain(range, projection)
            {
                for(var i = 0; i != projection.length; ++i)
                {
                    var segment = projection[i];
                    if(in_range(segment.range.bounds.min, range, segment.range.bounds.max))
                    {
                        if(segment.scale == "linear")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            return mix(segment.domain.min, segment.domain.max, amount)
                        }
                        else if(segment.scale[0] == "log")
                        {
                            var amount = (range - segment.range.min) / (segment.range.max - segment.range.min);
                            var base = segment.scale[1];
                            return sign(segment.domain.min) * Math.pow(base, mix(log(segment.domain.min, base), log(segment.domain.max, base), amount));
                        }
                    }
                }
            }

            var axes = {};

            function display_coordinates(e)
            {
                var current = canvas.createSVGPoint();
                current.x = e.clientX;
                current.y = e.clientY;

                for(var axis_id in axes)
                {
                    var axis = document.querySelector("#" + axis_id);
                    var coordinates = axis.querySelector(".toyplot-coordinates-Axis-coordinates");
                    if(coordinates)
                    {
                        var projection = axes[axis_id];
                        var local = current.matrixTransform(axis.getScreenCTM().inverse());
                        if(inside(local.x, projection))
                        {
                            var domain = to_domain(local.x, projection);
                            coordinates.style.visibility = "visible";
                            coordinates.setAttribute("transform", "translate(" + local.x + ")");
                            var text = coordinates.querySelector("text");
                            text.textContent = domain.toFixed(2);
                        }
                        else
                        {
                            coordinates.style.visibility= "hidden";
                        }
                    }
                }
            }

            canvas.addEventListener("click", display_coordinates);

            var module = {};
            module.show_coordinates = function(axis_id, projection)
            {
                axes[axis_id] = projection;
            }

            return module;
        })(modules["toyplot/canvas"]);
(function(tables, context_menu, io, owner_id, key, label, names, columns, filename)
        {
            tables.set(owner_id, key, names, columns);

            var owner = document.querySelector("#" + owner_id);
            function show_item(e)
            {
                return owner.contains(e.target);
            }

            function choose_item()
            {
                io.save_file("text/csv", "utf-8", tables.get_csv(owner_id, key), filename + ".csv");
            }

            context_menu.add_item("Save " + label + " as CSV", show_item, choose_item);
        })(modules["toyplot/tables"],modules["toyplot/menus/context"],modules["toyplot/io"],"t1b517aa90bb34042a4caa7369e9a6ca6","data","point",["x", "y0"],[[-19.268467964423504, -25.590282227818864, -19.3739774510354, 140.907334901688, -19.905250093092324, 140.74524919240352, -25.22530398468893, -27.731500193954645, -30.43838789545281, -30.362210944904064, -27.507576225074207, -28.256017727142698, -27.993609386500754], [32.327592455634594, -35.653907488942835, 43.00965710938873, -1.3586122420131528, -22.572685589468314, -1.388098448362023, -34.57298609726385, -41.19086751215904, -48.4523761179036, -47.61842677570342, -40.341697352090705, 98.93401978223604, 98.87838827664929]],"toyplot");
(function(axis, axis_id, projection)
        {
            axis.show_coordinates(axis_id, projection);
        })(modules["toyplot.coordinates.Axis"],"t71085bbdcfa0444d8c8bc12140533b6d",[{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 147.8682548903218, "min": -37.39930788408665}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 160.0, "min": 0.0}, "scale": "linear"}]);
(function(axis, axis_id, projection)
        {
            axis.show_coordinates(axis_id, projection);
        })(modules["toyplot.coordinates.Axis"],"t69864239977e4889a4791e7d64b9c4bc",[{"domain": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 104.9215921156792, "min": -54.439948451346766}, "range": {"bounds": {"max": Infinity, "min": -Infinity}, "max": 160.0, "min": 0.0}, "scale": "linear"}]);
</script></div></div>


Although it contains fewer SNPs overall, an
advantage of this dataset compared to the previous one is that fewer SNPs were imputed, and when they were imputed, population
information was applied so that the imputed allele is sampled from the allele frequency of samples only in its same population.