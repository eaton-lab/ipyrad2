# snpex

`ipyrad2 snpex` exports filtered SNP datasets. It always writes the core filtered matrices
used inside `ipyrad2`, and it can also write external-tool formats such as PLINK, PHYLIP,
NEXUS, FASTA, TreeMix, and EEMS from the same selected SNP view.

This is the SNP export command to use when you need explicit control over missing-data handling, linked versus unlinked SNPs, population-aware filtering with `imap` and `minmap`, or files for software outside the built-in `analysis` commands.

## When to Use

Use `snpex` when you need SNP files for software outside ipyrad2, such as:

- PLINK-style clustering, association, or preprocessing workflows
- SNP alignments for phylogenetic tools that accept PHYLIP, NEXUS, or FASTA
- TreeMix population-count input
- EEMS genetic dissimilarity input

If you are staying inside `ipyrad2 analysis`, the built-in methods already apply similar filtering and organize their own inputs internally, so you usually do not need `snpex`.

## Prerequisites

- an SNP-capable HDF5 file produced by `assemble` or `analysis vcf-to-hdf5`
- a filtering strategy for samples, populations, and linked versus unlinked SNPs
- a reference-aware HDF5 if you want PLINK export or global SNP imputation

## Inputs and Filtering

`snpex` starts from one SNP-capable HDF5 and applies the same SNP filtering logic used by the built-in numerical methods.

- `-d, --data`: input SNP HDF5
- `-m, --min-sample-coverage`: minimum number of samples with data required at a SNP
- `-r, --max-sample-missing`: drop samples whose missing-data fraction exceeds this threshold, then rerun SNP filtering
- `-a, --min-minor-allele-frequency`: drop low-frequency SNPs after coverage filtering
- `-e, --exclude`: exclude one or more named samples
- `-R, --include-reference`: include `assembly_reference_sequence`
- `-i, --imap`: sample-to-population map for subsetting and population-aware filtering
- `-g, --minmap`: per-population minimum coverage checks applied on top of `-m` when `imap` is used

By default `snpex` writes one SNP per RAD locus. Use `--no-subsample` to keep linked SNPs. Use `--seed` when you want reproducible one-SNP-per-locus subsampling.

## Imputation

By default `snpex` preserves missing genotypes. If you pass `--impute-method`, imputation is applied once before any outputs are written, and every written format uses that same imputed SNP view.

Supported imputation modes are:

- `sample`: sample missing diploid genotypes from within-group allele frequencies
- `zero-fill`: replace missing genotypes with homozygous reference calls

Because `snpex` always writes both diploid genotypes and SNP character matrices, global SNP imputation requires the HDF5 `reference` dataset so the imputed character matrix can be rebuilt consistently.

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

Use one or more per-format flags to add external-tool files:

- `--plink`: `NAME.bed`, `NAME.bim`, `NAME.fam`
- `--phylip`: `NAME.phy`
- `--nexus`: `NAME.nex`
- `--fasta`: `NAME.fa`
- `--treemix`: `NAME.treemix.gz`
- `--eems`: `NAME.eems`

`snpex` does not force you to choose only one export target. A single run can write several of these formats from the same selected SNP matrix.

## Format Notes

### PLINK

PLINK export writes BED/BIM/FAM from the selected SNP view. It requires the HDF5 `reference` dataset because the BIM file must name the reference and alternate alleles.

### PHYLIP, NEXUS, and FASTA

These formats write SNP characters as nucleotide/IUPAC alignments:

- homozygous reference sites write the reference base
- homozygous alternate sites write the alternate base
- heterozygotes write the corresponding IUPAC ambiguity code
- missing sites write `N` unless imputed beforehand

These are SNP alignments, not locus alignments. If you want sequence windows or whole loci instead of one-character-per-SNP matrices, use `wex` or `lex`.

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

### Baseline unlinked SNP export

```bash
ipyrad2 analysis snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/
```

### Keep linked SNPs and write PLINK plus phylogenetic alignments

```bash
ipyrad2 analysis snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  --no-subsample \
  --plink \
  --phylip \
  --nexus
```

### Apply one imputation method to every written output

```bash
ipyrad2 analysis snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  --plink \
  --fasta \
  --impute-method sample
```

### Export a population-filtered TreeMix dataset

```bash
ipyrad2 analysis snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -i pops.tsv \
  -g minmap.tsv \
  --treemix
```

### Export an EEMS genetic matrix from the same filtered SNP set

```bash
ipyrad2 analysis snpex \
  -d assembly.hdf5 \
  -o SNP_OUT/ \
  -i pops.tsv \
  -g minmap.tsv \
  --eems
```

## Common Failures

- Wrong input type: `snpex` expects SNP-capable HDF5, not raw VCF. Convert VCF first with `analysis vcf-to-hdf5`.
- Empty result after filtering: coverage, missingness, MAF, `imap`, and `minmap` can combine to remove all SNPs or all samples.
- Reference-related errors: PLINK export and global imputation both require the HDF5 `reference` dataset.
- IMAP mismatch: sample names in `imap` must match names in the HDF5 exactly.
- Existing output files: use `--force` to overwrite an existing export set.

## Related Pages

- [Writing Outputs](./index.md)
- [wex](./wex.md)
- [lex](./lex.md)
- [Analysis Guide](../analyses/index.md)
