# popgen

## Summary

`ipyrad2 analysis popgen` computes built-in population-genetic summary statistics from assembled sequence HDF5 or SNP-capable HDF5 inputs.

It supports two backends:

- Sequence-backed HDF5: `pi`, `dxy`, `fst`, `tajima_d`, `theta_w`, `heterozygosity`, `fis`, `sfs`
- SNP-backed HDF5: `fst`, `heterozygosity`, `fis`, `sfs`

Windowed analyses are currently sequence-only.

## When to Use

Use `popgen` when you want population-level summaries directly inside ipyrad2 rather than exporting to another package first.

This command is most useful when you want:

- quick within-population diversity summaries
- pairwise divergence or differentiation summaries between populations
- a folded, genome-wide site-frequency spectrum
- optional windowed summaries on sequence-backed assemblies

## Prerequisites

- An input HDF5 that contains either sequence data, SNP data, or both
- A population map file if you want multi-population summaries
- An understanding of how much missing data remains in your assembled dataset

If you do not provide `--imap`, all retained samples are treated as one population named `all`.

## Inputs and Backends

`popgen` chooses its backend from the datasets stored in the HDF5:

- Sequence backend: used when `phy` and `phymap` are present
- SNP backend: used when sequence data are absent but `genos` and `snpsmap` are present

The sequence backend supports the full current statistic panel. The SNP backend supports only the SNP-backed subset because several estimators in this phase require aligned sequence data.

`--stats all` means:

- On sequence-backed input: run `pi`, `dxy`, `fst`, `tajima_d`, `theta_w`, `heterozygosity`, `fis`, and `sfs`
- On SNP-backed input: run `fst`, `heterozygosity`, `fis`, and `sfs`

Windowing is only available on sequence-backed HDF5. If you request genomic or locus windows on SNP-only HDF5, the command errors.

## Population Grouping and Filtering

Population membership is controlled with `--imap`. Each retained sample must belong to exactly one population for the per-sample and per-population summary tables.

`--minmap` controls the per-population minimum number of called samples required for site-based summaries. This matters most for pairwise and within-population summaries when missingness differs across populations.

The main dataset filters are:

- `--min-sample-coverage`: remove sites that do not meet the required sample coverage
- `--max-sample-missing`: drop samples with too much missing data before the analysis
- `--min-minor-allele-frequency`: remove low-frequency SNPs in SNP-backed runs
- `--exclude`: remove named samples before population summaries are calculated
- `-R, --include-reference`: include the assembly reference sequence if present
- `--subsample-unlinked`: for SNP-backed runs, subsample one SNP per linkage block
- `--seed`: control the random seed for unlinked SNP subsampling

Missing data are not imputed in `popgen`.

Current missing-data behavior is:

- Sequence backend: `N` and gap `-` are treated as missing
- SNP backend: missing genotype calls are excluded from allele-count summaries

As a general rule, stronger filtering is usually safer than imputation for these estimators. Most of these statistics are designed around observed allele counts, and imputation can distort rare alleles, heterozygosity, and population differences.

## Command Patterns

Sequence-backed full panel for grouped samples:

```bash
ipyrad2 analysis popgen \
  -d assembly.hdf5 \
  -i populations.tsv \
  --stats all \
  -o POPGEN_OUT
```

SNP-backed heterozygosity, `Fis`, `FST`, and SFS with one unlinked SNP per locus:

```bash
ipyrad2 analysis popgen \
  -d snps.hdf5 \
  -i populations.tsv \
  --stats all \
  --subsample-unlinked \
  --seed 7 \
  -o POPGEN_OUT
```

Sequence-backed genomic windows:

```bash
ipyrad2 analysis popgen \
  -d assembly.hdf5 \
  -i populations.tsv \
  --stats pi,fst,fis \
  --window-size 50000 \
  --step-size 10000 \
  -o POPGEN_WINDOWS
```

Sequence-backed anonymous RAD-style locus windows:

```bash
ipyrad2 analysis popgen \
  -d assembly.hdf5 \
  -i populations.tsv \
  --stats pi,dxy,fst \
  --loci-per-window 50 \
  --locus-step 10 \
  -o POPGEN_LOCUS_WINDOWS
```

## Outputs

`popgen` writes a manifest plus one or more TSV tables, depending on which statistics were requested.

Current outputs are:

- `<name>.manifest.txt`
- `<name>.sample_stats.tsv`
- `<name>.population_stats.tsv` when population-level scalar statistics are requested
- `<name>.pairwise_stats.tsv` when pairwise statistics are requested
- `<name>.sfs.tsv` when `sfs` is requested
- `<name>.window_population_stats.tsv` when windowed population summaries are requested
- `<name>.window_pairwise_stats.tsv` when windowed pairwise summaries are requested

Two important current behaviors:

- `popgen` does not write a separate `stats.txt`
- `popgen` does not currently write `sample_data_summary.tsv` as a separate file

Instead, the manifest embeds a `Sample Data Summary` section showing:

- `sample`
- `missing_fraction`
- `post_imputation_missing_fraction`
- `imputation_algorithm`
- `imputed_genotype_fraction`

For `popgen`, that embedded table documents that the data were not imputed.

### `sample_stats.tsv`

This is the per-sample summary table. It records:

- the population assignment for each retained sample
- total sites examined for that sample
- called and missing site counts
- called and missing fractions
- homozygous and heterozygous called-site counts
- observed heterozygosity per sample

This table is often the first place to check for sample-level missingness problems before interpreting population summaries.

### `population_stats.tsv`

This table contains within-population summaries. The exact columns depend on the requested statistics, but can include:

- `pi`
- `theta_w`
- `tajima_d`
- `observed_heterozygosity`
- `expected_heterozygosity`
- `fis`

It also includes the site counts actually used for each class of estimator, which matters because different statistics can retain different site sets under missing-data filters.

### `pairwise_stats.tsv`

This table contains pairwise population comparisons. Depending on the backend and requested statistics it can include:

- `dxy`
- `fst`

Each row also reports `sites_used`, which is critical for interpretation when populations have uneven missing data.

### `sfs.tsv`

This table reports a folded, genome-wide site-frequency spectrum by population using:

- `population`
- `minor_allele_count`
- `site_count`

It is currently genome-wide only. Windowed SFS is not written in this phase.

## Windowed Analyses

Windowing is currently supported only for sequence-backed HDF5.

There are two windowing modes:

- Genomic windows: `--window-size` with optional `--step-size`
- Consecutive-locus windows: `--loci-per-window` with optional `--locus-step`

Genomic windows use scaffold coordinates. Locus windows slide across consecutive `phymap` loci and are useful for RAD-style data where fixed genomic spans may be less informative.

Windowed runs still write the genome-wide outputs. The window tables are additional outputs, not replacements.

Current window metadata columns include:

- `window_id`
- `window_mode`
- `scaffold`
- `start`
- `end`
- `first_locus`
- `last_locus`
- `nloci`
- `sites_total`

`window_mode` is either:

- `genomic`
- `locus`

## Statistic Guide

### `pi`

`pi` is within-population nucleotide diversity, or mean pairwise difference per site within a population.

In `ipyrad2`, `pi` is calculated from per-site allele counts using unbiased gene diversity and then averaged across retained sites with enough called chromosomes. For the sequence backend, ambiguous IUPAC heterozygote codes contribute split allele counts. For the SNP backend, `pi` is not currently reported.

Missing data can affect `pi` by changing which sites remain analyzable and by reducing the number of called chromosomes at retained sites. In practice, stronger filtering is usually preferable to imputation for `pi`, because imputation can flatten real allele-frequency variation.

### `dxy`

`dxy` is the mean sequence divergence between two populations.

`ipyrad2` calculates sitewise `dxy` as `1 - sum(p_i q_i)`, where `p_i` and `q_i` are allele frequencies in the two populations, and then averages across retained sites shared by both populations.

Missing data can matter a lot for `dxy` if one population loses more sites than the other, because the effective comparison set changes. Users should generally filter more strictly rather than impute when using `dxy`.

### `fst`

`fst` is a pairwise measure of population differentiation. In this implementation it is Hudson-style `FST`, computed by summing the sitewise numerator and denominator and then taking the ratio of sums.

`ipyrad2` uses the sitewise components:

- numerator: `dxy - (pi1 + pi2) / 2`
- denominator: `dxy`

and reports `FST = sum(numerator) / sum(denominator)` across retained sites.

Reference: Hudson, Slatkin, and Maddison 1992.

Missing data can strongly affect `FST`, especially when coverage differs across populations, because the retained site set may shift in ways that look like structure. Strong missing-data filtering is usually safer than imputation here.

### `theta_w`

`theta_w` is Watterson's theta, a diversity estimator based on the number of segregating sites.

`ipyrad2` calculates per-site Watterson's theta from the number of segregating sites divided by the harmonic-number constant `a1`, then normalizes by the number of sites used. In the current sequence backend, this estimator is computed from fully called sites within each population.

Reference: Watterson 1975.

Missing data strongly affects `theta_w` because segregating-site counts are very sensitive to incomplete site retention. Users should generally filter more strictly and avoid imputation for this statistic.

### `tajima_d`

Tajima's D compares two diversity summaries: diversity from pairwise differences and diversity from the number of segregating sites.

`ipyrad2` calculates Tajima's D from total pairwise diversity, the number of segregating sites, and sample size using the standard Tajima formula. As with `theta_w`, the current implementation uses fully called sites within each population for this estimator.

Reference: Tajima 1989.

Missing data strongly affects Tajima's D because both the segregating-site count and the usable-site set change. Strong filtering is usually preferable to imputation.

### `heterozygosity`

`heterozygosity` reports two related summaries:

- observed heterozygosity (`Ho`)
- expected heterozygosity (`He`)

Observed heterozygosity is the fraction of called genotypes that are heterozygous. Expected heterozygosity is calculated from allele counts as unbiased gene diversity.

This statistic is available on both sequence-backed and SNP-backed data. On sequence-backed data, heterozygous IUPAC ambiguity codes are treated as heterozygous calls.

Missing data has a moderate effect because both `Ho` and `He` depend on the called genotype set. In most cases, filtering is preferred to imputation because imputation can inflate or suppress heterozygote counts.

### `fis`

`fis` is the within-population inbreeding coefficient.

`ipyrad2` calculates it as:

```text
Fis = 1 - Ho / He
```

where `Ho` is observed heterozygosity and `He` is expected heterozygosity. This same formula is also written to the manifest when `fis` is requested.

Because `fis` is a ratio of heterozygosity summaries, it can be sensitive to missing data and especially unstable when expected heterozygosity is very low. Strong filtering is generally better than imputation.

### `sfs`

`sfs` is a folded, genome-wide site-frequency spectrum. In the current implementation it is reported as counts of biallelic sites by minor-allele-count bin within each population.

This means the output is not an unfolded ancestral-state spectrum. It is a folded spectrum based on minor allele count, and only genome-wide SFS is currently written.

Missing data matters a great deal for the SFS, especially in the rare-allele bins, because missing calls change which sites remain biallelic and how minor-allele counts are tallied. Users should usually apply fairly strong missing-data filtering and avoid imputation for SFS calculations.

## References

- Watterson, G. A. 1975. On the number of segregating sites in genetical models without recombination.
- Tajima, F. 1989. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.
- Hudson, R. R., Slatkin, M., and Maddison, W. P. 1992. Estimation of levels of gene flow from DNA sequence data.

## Example Outputs

### Example manifest

```text
# paste real example manifest text here
```

### Example `sample_stats.tsv`

```text
# paste real sample_stats.tsv text here
```

### Example `population_stats.tsv`

```text
# paste real population_stats.tsv text here
```

### Example `pairwise_stats.tsv`

```text
# paste real pairwise_stats.tsv text here
```

### Example `sfs.tsv`

```text
# paste real sfs.tsv text here
```

### Example `window_population_stats.tsv`

```text
# paste real window_population_stats.tsv text here
```

### Example `window_pairwise_stats.tsv`

```text
# paste real window_pairwise_stats.tsv text here
```

## Related Pages

- [Analysis Guide](./index.md)
- [assemble](../assembly/assemble.md)
- [Writing Outputs](../writing-outputs/index.md)
- [Files and Data Types](../getting-started/files-and-data-types.md)
