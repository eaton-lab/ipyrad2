# pca

## Summary

`ipyrad2 analysis pca` is the PCA-family command for SNP-capable HDF5 data. It currently supports three numerical methods:

- `pca`
- `tsne`
- `umap`

All three methods start from the same filtered SNP view, require a fully imputed genotype matrix, and write coordinate tables for downstream inspection. PCA can also write a default SVG plot when requested.

In ipyrad2, this command means:

1. filter SNPs and samples from an SNP-capable HDF5
2. optionally subsample one SNP per RAD locus
3. impute all remaining missing genotypes
4. run PCA, t-SNE, or UMAP on the resulting matrix
5. write coordinates plus method-specific summaries

By default this command writes numerical result tables only. If you add `--plot` with `-M pca`, it also writes a default SVG PCA figure.

## When to Use

Use `analysis pca` when you want:

- principal component coordinates from filtered SNP data
- quick nonlinear embeddings with t-SNE or UMAP from the same filtered input logic
- a light-weight way to compare structure across linked versus unlinked SNP views
- per-sample missing-data and imputation summaries alongside the numerical output

`pca` is usually the best first pass. `tsne` and `umap` are useful when you want exploratory embeddings, but they are more sensitive to parameter choices and should be interpreted more cautiously.

## Prerequisites

- an SNP-capable HDF5 file produced by `assemble` or `analysis vcf-to-hdf5`
- enough retained samples and SNPs after filtering to support the selected method
- scikit-learn installed for all PCA-family methods
- `umap-learn` installed if you want `-M umap`
- `toyplot` installed if you want `--plot`

If `umap-learn` is not installed, only the UMAP method errors. If `toyplot` is not installed, only `--plot` errors. Ordinary PCA, t-SNE, and UMAP runs do not require toyplot unless plotting is requested.

## Inputs and Filtering

`analysis pca` only accepts SNP-capable HDF5 input. It does not run directly from VCF.

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
- `--perplexity` and `--max-iter`: only used with `-M tsne`
- `--n-neighbors`: only used with `-M umap`

Important current rules:

- PCA supports one or more replicates
- PCA plotting currently supports only the first two principal components
- t-SNE supports exactly one run
- UMAP supports exactly one run

For t-SNE:

- `perplexity` must be greater than zero
- `perplexity` must be smaller than the number of retained samples
- `max_iter` must be at least 250

For UMAP:

- `n_neighbors` must be at least 2

## Imputation and Missing Data

All PCA-family methods in ipyrad2 require a fully imputed genotype matrix. Missing genotypes are not allowed to remain in the numerical matrix passed to PCA, t-SNE, or UMAP.

Supported imputation modes are:

- `sample`: sample missing diploid genotypes from allele frequencies within each imputation group
- `zero-fill`: replace missing genotypes with homozygous reference calls
- `none`: currently accepted as a deprecated alias for `zero-fill`

`sample` is the default and generally the preferred choice.

If you provide `imap`, sample-mode imputation uses those groups. If you do not provide `imap`, all retained samples are treated as one imputation group.

`zero-fill` is usually a poor default for exploratory structure analyses because it can pull missing-heavy samples toward the reference state. Since `none` currently behaves the same way, users should not treat it as “leave missing data unchanged.”

Users should inspect `sample_data_summary.tsv` after every run and consider dropping samples that require too much imputation before trusting the structure. Heavy imputation is often a warning that missingness, rather than biology, may be shaping the ordination.

## How ipyrad2 PCA-Family Works

The active implementation is:

1. filter the SNP HDF5 with the shared SNP extracter
2. select a linked or unlinked SNP view
3. impute the genotype matrix
4. run the chosen numerical method on the imputed matrix
5. write coordinates and summary tables

Method-specific details:

- `pca`: uses a direct SVD-based PCA on the centered matrix and writes explained-variance ratios
- `tsne`: uses scikit-learn t-SNE with `init="pca"`
- `umap`: uses `umap-learn` UMAP with `init="spectral"`

PCA replicates are mainly useful when the selected SNP view can change across replicates, especially under one-SNP-per-locus subsampling and stochastic imputation. PCA itself is otherwise deterministic for a fixed matrix.

## Outputs

Every run writes:

- `<name>.coords.tsv`
- `<name>.sample_data_summary.tsv`
- `<name>.stats.txt`

PCA runs additionally write:

- `<name>.variance.tsv`

If you add `--plot`, PCA also writes:

- `<name>.plot.svg`

### `coords.tsv`

This table contains one row per sample per replicate:

- `sample`
- `replicate`
- `method`
- `axis1`
- `axis2`
- additional axes as available

For PCA with multiple replicates, each sample appears once per replicate. For t-SNE and UMAP, the replicate column is still present, but only replicate `0` is written.

### `sample_data_summary.tsv`

This table records per-sample missingness and imputation:

- `sample`
- `missing_fraction`
- `post_imputation_missing_fraction`
- `imputation_algorithm`
- `imputed_genotype_fraction`

For PCA-family methods, `post_imputation_missing_fraction` should be `0` because the numerical matrix is fully imputed.

For PCA with multiple replicates, the numeric columns in this file are averaged across replicates. The imputation algorithm column stays constant.

### `variance.tsv`

This file is written only for `-M pca`. It contains:

- `replicate`
- `axis`
- `explained_variance_ratio`

t-SNE and UMAP do not write a variance file.

### `plot.svg`

This file is written only when you use `--plot` with `-M pca`.

- it is a default SVG scatter plot of `PC1` versus `PC2`
- the plotting axes are drawn with an external-tick style and a boxed outline
- points are colored by `imap` group, or by the default `all` group if no `imap` is provided
- single-replicate runs show one point per sample
- multi-replicate runs show translucent replicate clouds plus one centroid per sample

This plot is meant as a quick first look, not a replacement for checking the coordinate tables directly.

### `stats.txt`

This is a human-readable summary file. It includes:

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
- Invalid plotting request: `--plot` currently works only with `-M pca`.
- Existing output files: use `--force` to overwrite an existing result set.
- Missing dependencies: scikit-learn is required for all methods, `umap-learn` is additionally required for UMAP, and `toyplot` is only required when `--plot` is requested.

## Related Pages

- [Analysis Guide](./index.md)
- [dapc](./dapc.md)
- [popgen](./popgen.md)
- [Writing Outputs](../writing-outputs/index.md)
- [Files and Data Types](../getting-started/files-and-data-types.md)
