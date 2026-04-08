# dapc

## Summary

`ipyrad2 analysis dapc` runs a sklearn-backed DAPC-style clustering workflow on SNP-capable HDF5 data.

In ipyrad2, this workflow means:

1. filter SNPs and samples from an SNP-capable HDF5
2. optionally subsample one SNP per RAD locus
3. impute all remaining missing genotypes
4. run PCA on the genotype matrix
5. cluster the PCA scores with K-means
6. fit linear discriminant analysis to those inferred cluster labels

This is a DAPC-style method, not a claim of full feature parity with adegenet. The current command writes numerical result tables only. It does not make plots.

## When to Use

Use `dapc` when you want:

- a clustering-oriented summary of structure in SNP data
- low-dimensional discriminant coordinates for retained samples
- cluster membership tables and hard assignments
- a quick internal `K` scan using K-means BIC

This command is a good fit when you want to stay inside ipyrad2 and inspect cluster separation without exporting to another package first.

## Prerequisites

- an SNP-capable HDF5 file produced by `assemble` or `analysis vcf-to-hdf5`
- enough retained samples and SNPs after filtering to support PCA and clustering
- scikit-learn installed in the environment

If scikit-learn is missing, the command errors and tells you to install the analysis extras.

## Inputs and Filtering

`dapc` only accepts SNP-capable HDF5 input. It does not run directly from VCF.

The main filtering and sample-selection controls are:

- `-m, --min-sample-coverage`: minimum number of samples with data required at a SNP
- `-r, --max-sample-missing`: drop samples whose missing-data fraction exceeds this threshold, then rerun SNP filtering
- `-a, --min-minor-allele-frequency`: remove low-frequency SNPs after coverage filtering
- `-e, --exclude`: exclude one or more named samples
- `-R, --include-reference`: include `assembly_reference_sequence`
- `-i, --imap`: sample-to-population map for filtering and imputation grouping
- `-g, --minmap`: per-population minimum coverage checks applied on top of `-m` when `imap` is used

Two important points:

- By default `dapc` subsamples one SNP per RAD locus. Use `--no-subsample` to keep linked SNPs.
- `imap` does not provide supervised class labels for DAPC in this implementation. Clusters are inferred internally by K-means. `imap` is only used for filtering and for population-aware sample-mode imputation.

`--seed` affects SNP subsampling, imputation, PCA random state, and K-means initialization.

## Imputation and Missing Data

DAPC in ipyrad2 always runs on a fully imputed genotype matrix. Missing genotypes are not allowed to remain in the matrix that enters PCA.

Supported imputation modes are:

- `sample`: sample missing diploid genotypes from allele frequencies within each imputation group
- `none`: current compatibility alias for zero-fill, which replaces missing genotypes with homozygous reference calls

`sample` is the default and generally the preferred choice.

If you provide `imap`, sample-mode imputation uses those groups. If you do not provide `imap`, all retained samples are treated as one imputation group.

`none` is generally a bad choice for clustering. In the current implementation it does not preserve missing data. It behaves like zero-fill, which can pull missing-heavy samples toward the reference state and distort both PCA structure and final cluster separation.

Users should inspect `sample_data_summary.tsv` after every run and consider dropping samples that require too much imputation before trusting the clustering result. Heavy imputation is a data-quality warning, not just a preprocessing detail.

## Choosing K and PCA Axes

You must provide exactly one of:

- `-k` for a fixed number of clusters
- `--k-range MIN:MAX` for a scan across multiple `K` values

Current `K` rules are:

- `K` must be at least 2
- `K` must be smaller than the number of retained samples
- for `--k-range`, the maximum `K` must also be smaller than the number of retained samples

When `--k-range` is used, ipyrad2:

1. runs PCA once
2. fits K-means for each candidate `K`
3. scores each fit with a lower-is-better BIC approximation
4. selects the `K` with the lowest BIC

`k_scan.tsv` is always written. Fixed-`K` runs produce a one-row table. `K`-range runs produce one row per tested `K`, with a `selected` column marking the chosen solution.

`--n-pcs` controls how many PCA axes are retained before discriminant analysis.

If `--n-pcs` is omitted, the current default is:

- at most 20 PCs
- never more than the available PCs
- never fewer than `Kmax - 1`

If `--n-pcs` is provided, it must be:

- at least `Kmax - 1`
- no greater than the available number of PCs

## How ipyrad2 DAPC Works

The active implementation is:

1. filter the SNP HDF5 with the shared SNP extracter
2. select a linked or unlinked SNP view
3. impute the genotype matrix
4. run PCA on the imputed matrix
5. run K-means on the PCA scores
6. fit LDA on those PCA scores using the K-means labels
7. write discriminant coordinates and cluster membership summaries

This means the discriminant analysis is clustering-first and supervised-second. The LDA labels come from K-means, not from user-supplied populations.

The membership values written by `dapc` are LDA class probabilities for those inferred clusters. They should not be interpreted as admixture proportions from a generative ancestry model.

## Outputs

Every run writes:

- `<name>.coords.tsv`
- `<name>.membership.tsv`
- `<name>.assignments.tsv`
- `<name>.k_scan.tsv`
- `<name>.sample_data_summary.tsv`
- `<name>.stats.txt`

### `coords.tsv`

This table contains one row per retained sample and one column per discriminant axis:

- `sample`
- `axis1`
- `axis2`
- additional axes as available

For small `K`, there may be only one or a few discriminant axes.

### `membership.tsv`

This table contains soft cluster membership values:

- `sample`
- `cluster1`
- `cluster2`
- additional cluster columns up to the selected `K`

These values come from `lda.predict_proba(...)` on the PCA scores.

### `assignments.tsv`

This table contains hard assignments derived from the membership matrix:

- `sample`
- `assigned_cluster`
- `assignment_score`

`assignment_score` is the largest membership value for that sample.

### `k_scan.tsv`

This table records the tested `K` values and the BIC score used for model selection:

- `k`
- `bic`
- `selected`

Lower BIC is better in the current implementation.

### `sample_data_summary.tsv`

This table is especially important for DAPC because missing data are always imputed before PCA:

- `sample`
- `missing_fraction`
- `post_imputation_missing_fraction`
- `imputation_algorithm`
- `imputed_genotype_fraction`

For DAPC, `post_imputation_missing_fraction` should be `0` because the numerical matrix is fully imputed.

Use this table to identify samples with unusually heavy imputation before interpreting cluster membership or plotting discriminant axes downstream.

### `stats.txt`

This is a human-readable summary file. It includes:

- tool name and input file
- initial, dropped, and retained samples
- `imap`, `minmap`, and reference-inclusion settings
- whether SNPs were subsampled
- random seed
- imputation method
- selected `K`
- requested `K` range or `NA`
- retained number of PCs
- counts of filtered and exported SNPs or linkage blocks
- shared SNP-extracter filter statistics

## Command Patterns

Basic fixed-`K` run with the default imputation strategy:

```bash
ipyrad2 analysis dapc \
  -d snps.hdf5 \
  -o DAPC_OUT \
  -k 2
```

Scan a range of `K` values and let BIC choose:

```bash
ipyrad2 analysis dapc \
  -d snps.hdf5 \
  -o DAPC_OUT \
  --k-range 2:5
```

Set an explicit number of retained PCA axes:

```bash
ipyrad2 analysis dapc \
  -d snps.hdf5 \
  -o DAPC_OUT \
  -k 3 \
  --n-pcs 10
```

Keep linked SNPs instead of using one SNP per RAD locus:

```bash
ipyrad2 analysis dapc \
  -d snps.hdf5 \
  -o DAPC_OUT \
  -k 3 \
  --no-subsample
```

Apply population-aware filtering and sample-mode imputation groups:

```bash
ipyrad2 analysis dapc \
  -d snps.hdf5 \
  -o DAPC_OUT \
  --k-range 2:4 \
  -i pops.tsv \
  -g minmap.tsv
```

## Interpretation Notes

- The discriminant axes are optimized to separate inferred clusters. They should not be interpreted the same way as PCA axes.
- `K`, `--n-pcs`, SNP linkage, and missing-data handling can all change the apparent clustering strongly.
- Keeping linked SNPs can overweight dense loci. The default one-SNP-per-locus subsampling is there to reduce that problem.
- Membership values are useful for summarizing assignment uncertainty, but they are not admixture coefficients in the same sense as model-based ancestry tools.
- Samples with high `imputed_genotype_fraction` deserve extra scrutiny. If many genotypes had to be imputed, cluster placement may reflect missing-data structure as much as biology.

## Common Failures

- Wrong input type: `dapc` expects SNP-capable HDF5, not raw VCF.
- Empty or too-small result after filtering: DAPC needs at least two retained samples and at least one retained SNP after filtering.
- Invalid `K`: `K` must be at least 2 and smaller than the number of retained samples.
- Invalid `--n-pcs`: it cannot be smaller than `Kmax - 1` or larger than the available PCs.
- Existing output files: use `--force` to overwrite an existing result set.
- Missing dependency: scikit-learn must be installed for this command.

## Related Pages

- [Analysis Guide](./index.md)
- [popgen](./popgen.md)
- [Writing Outputs](../writing-outputs/index.md)
- [Files and Data Types](../getting-started/files-and-data-types.md)
