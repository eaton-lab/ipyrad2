# snpex

`ipyrad2 snpex` can be used to export subsampled and filtered SNP datasets in a variety of
common formats, including PLINK, PHYLIP, NEXUS, and FASTA, as well as several more complex supported
formats, including TreeMix and EEMS that group individual data into population frequencies.
The subsampling, filtering, and optional imputation options allow for easily generating multiple
datasets with different thresholds for dealing with missing data.

## Quick guide

!!! tip "Use snpex to filter/subsample/impute and write SNP data"

```bash
ipyrad2 snpex \
  --data HDF5 \
  --out SNPS/ \
  --min-sample-coverage 20 \
  --max-sample-missing 0.5 \
  --out-format plink fasta
```

Note: If you are using [Analysis](../analyses) tools inside `ipyrad2`, sich as its `pca`, `popgen`, `baba`, and related tools,
you do not need to use `snpex` as its filtering options are implemented within those commands directly.

## Input and output

- `-d, --data`: input SNP HDF5
- `-o, --out`: ...
- `-O, --out-format`: write one or more additional formats selected from `plink`, `phylip`, `nexus`, `fasta`, `treemix`, and `eems`
- `-i, --imap`: sample-to-population map for subsetting and population-aware filtering
- `-g, --minmap`: per-population minimum coverage checks applied on top of `-m` when `imap` is used

## Filtering

If an `--imap` is entered the dataset will be subset to include only the samples listed in the file.

- `-m, --min-sample-coverage`: minimum number of samples with data required at a SNP
- `-r, --max-sample-missing`: drop samples whose missing-data fraction exceeds this threshold, then rerun SNP filtering
- `-a, --min-minor-allele-frequency`: drop low-frequency SNPs after coverage filtering
- `-e, --exclude`: exclude one or more named samples
- `-R, --include-reference`: include `assembly_reference_sequence`
- `--no-subsample`: By default `snpex` writes one SNP per RAD locus. Use `--no-subsample` to keep linked SNPs.
- `--seed` Use this when you want reproducible one-SNP-per-locus subsampling.

## Imputation

By default `snpex` preserves missing genotypes. If you pass `--impute-method`, imputation is applied
once before any outputs are written, and every written format uses that same imputed SNP view.

Supported imputation modes are:

- `sample`: sample missing diploid genotypes from within-group allele frequencies
- `zero-fill`: replace missing genotypes with homozygous reference calls


## Core Outputs

Every run writes these baseline files:

- `NAME.genos.npy`: diploid genotype matrix coded as `0`, `1`, `2`, or `255` for missing unless imputed
- `NAME.snps.npy`: SNP character matrix using nucleotide and IUPAC SNP codes
- `NAME.snpsmap.tsv`: SNP metadata with `loc`, `loc_idx`, `loc_pos`, `scaff`, and `pos`
- `NAME.samples.txt`: retained sample order
- `NAME.sample_data_summary.tsv`: per-sample missingness before and after optional imputation
- `NAME.stats.txt`: export summary, filter statistics, imputation summary, and written formats

These baseline outputs are useful even if you do not request any additional external format.

## Optional Export Formats

Use `-O/--out-format` followed by one or more format names to add external-tool files:

- `plink`: `NAME.bed`, `NAME.bim`, `NAME.fam`
- `phylip`: `NAME.phy`
- `nexus`: `NAME.nex`
- `fasta`: `NAME.fa`
- `treemix`: `NAME.treemix.gz`
- `eems`: `NAME.eems`

For example, `-O plink fasta` writes both formats from the same selected SNP matrix. If `-O` is omitted, `snpex` writes only the core outputs listed above.

## Format Notes

### PLINK

PLINK export writes BED/BIM/FAM from the selected SNP view.

### PHYLIP, NEXUS, and FASTA

These formats write SNP characters as nucleotide/IUPAC alignments:

- homozygous reference sites write the reference base
- homozygous alternate sites write the alternate base
- heterozygotes write the corresponding IUPAC ambiguity code
- missing sites write `N` unless imputed beforehand

These are SNP alignments, not locus alignments. If you want clipped sequence windows or whole loci instead of one-character-per-SNP matrices, use `seqex`.

### TreeMix

TreeMix export writes gzipped counts in the usual `ancestral,derived` text format.

- if you supplied `imap`, columns are populations
- if you did not supply `imap`, each retained sample is treated as its own column

### EEMS

EEMS export writes only the genetic dissimilarity matrix, named `NAME.eems`.

- rows and columns follow the sample order in `NAME.samples.txt`
- if global SNP imputation was requested, the EEMS matrix is built from that imputed genotype matrix
- otherwise `snpex` follows the usual EEMS SNP convention of filling missing genotypes with per-site means before forming the pairwise dissimilarity matrix

This command does not write the spatial `.coord` or habitat `.outer` files. Those remain external inputs you prepare separately.

## Common Patterns

### filter for missing data and write unlinked SNPs matrix

```bash
ipyrad2 snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -m 10
```

### filter for missing data and write linked SNPs to PLINK format

```bash
ipyrad2 snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -m 10 \
  --no-subsample \
  -O plink
```

### filter, impute, and write to FASTA

```bash
ipyrad2 snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -m 10 \
  --impute-method sample \
  -O fasta
```

### filter and impute population-aware and write to TREEMIX

```bash
ipyrad2 snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -i pops.tsv \
  -g minmap.tsv \
  -O treemix
```

## Common Failures

- Wrong input type: `snpex` expects SNP-capable HDF5, not raw VCF. Convert VCF first with `analysis vcf-to-hdf5`.
- Empty result after filtering: coverage, missingness, MAF, `imap`, and `minmap` can combine to remove all SNPs or all samples.
- Reference-related errors: PLINK export and global imputation both require the HDF5 `reference` dataset.
- IMAP mismatch: sample names in `imap` must match names in the HDF5 exactly.
- Existing output files: use `--force` to overwrite an existing export set.

## Related Pages

- [Writing Outputs](./index.md)
- [seqex](./seqex.md)
- [Analysis Guide](../analyses/index.md)
