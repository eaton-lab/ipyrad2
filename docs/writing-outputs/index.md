# Writing Outputs

The export tools in `ipyrad2` are designed to allow you to assemble a dataset just once, and to
use the generated HDF5 database file to export many curated datasets for downstream analyses.
The export tools, `seqex` and `snpex`, allow you to filter the dataset by scaffold,
sample or population names, and critically, by patterns of missing data, and to write the
resulting data to output files formatted for various external tools.

This design is especially useful for RAD-seq data because missing data usually needs to be handled explicitly rather than hidden. Different downstream methods tolerate different patterns of missing
samples, missing sites, linked SNPs, and genotype imputation. The export layer makes those choices
easy to control for each target analysis instead of forcing one global preprocessing decision on
every output.

Similar data filtering patterns are also implemented in the `ipyrad2` [Analysis](Analysis) tools,
but these cover only a subset of downstream analyses for RAD-seq data. Using the export tools, you can create formatted datasets for most common downstream sequence- and SNP-based software.

## Why This Layer Exists

RAD-seq datasets often need different filtering choices depending on the software you plan to run next. In practice that means deciding, for each export:

- which samples to keep or exclude
- whether to include `assembly_reference_sequence`
- how much missing data is acceptable per sample
- how much site or SNP coverage is required
- whether population-aware minimum coverage rules should apply
- whether linked SNPs should be subsampled to one SNP per RAD locus
- whether missing genotypes should stay missing or be imputed for a specific file format

## The Export Tools

- [`seqex`](./seqex.md): writes multi-locus, concatenated, or split alignments from complete loci or exact coordinate intersections after applying locus-wise and site-wise filters.
- [`snpex`](./snpex.md): writes filtered SNP matrices after applying site-wise filters, with options for subsampling unlinked SNPs and exporting to various formats including PLINK, GENO, PHYLIP, TreeMix, and EEMS.


## Command examples

```bash
# write concatenated alignment of chromosome 1 in a reference assembly
ipyrad2 seqex \
    -d DATA.hdf5 \
    -o ALIGNMENTS/ \
    -n Chr01 \
    -w Chr01 \
    --min-sample-coverage 10 \
    --concatenate

# select 1000 random denovo loci and concatenate into an alignment
ipyrad2 seqex \
    -d DATA.hdf5 \
    -o ALIGNMENTS/ \
    -n loci_10K \
    -N 1000 \
    --min-sample-coverage 10 \
    --concatenate

# select 1000 unlinked SNPs random denovo loci and concatenate into an alignment
ipyrad2 snpex \
    -d DATA.hdf5 \
    -o ALIGNMENTS/ \
    -n snps_10K \
    -N 1000 \
    --min-sample-coverage 10

```

<!--
## Comparison

| Command | Input | Output | Main filtering controls | Missing-data handling | Imputation availability | Typical downstream use |
| --- | --- | --- | --- | --- | --- | --- |
| `seqex` | HDF5 | one or many locus alignments | windows, number of loci, minimum locus length, `-m`, `-r`, `-e`, `-R`, `imap`, `minmap` | explicit locus, sample, and site filtering before export | none | multilocus sequence analyses and locus-based phylogenetic workflows |
| `snpex` | HDF5 | filtered SNP matrices or supported formatted  | `-m`, `-r`, `-a`, `-e`, `-R`, `imap`, `minmap`, linked vs unlinked SNPs | explicit sample and SNP filtering before export | optional, and applied to every written output in that run | SNP-based external tools such as PLINK, TreeMix, EEMS, or SNP-alignment workflows |
-->


## Missing Data, Linkage, and Imputation

The key idea is not simply to export data, but to export the right version of the data for the next tool.

For sequence exports, `seqex` lets you control which loci, samples, and sites survive filtering before alignments are written. That makes missing data explicit at the alignment stage instead of leaving it to downstream programs to interpret blindly.

For SNP exports, `snpex` adds the choices that matter most for RAD data:

- one-SNP-per-locus subsampling by default to reduce linkage effects
- an explicit option to keep linked SNPs when that is appropriate
- sample and SNP filtering before matrix export
- optional output formats for several external tools from the same SNP selection
- optional global SNP imputation, when the downstream tool expects fully filled genotypes or SNP characters

That means you can create one SNP export for a linkage-sensitive method, another for a program that wants linked SNPs retained, another that needs phylogenetic SNP alignments, and another that needs imputed genotype files, all from the same structured dataset.

## Writing Outputs vs Analysis Guide

Use Writing Outputs when the next step happens outside ipyrad2 and you need files in a format another program will read directly.

If you plan to stay inside ipyrad2, these export tools are often unnecessary. The built-in methods in the [Analysis Guide](../analyses/index.md) support similar filtering ideas, but they handle dataset preparation and internal file organization for you. In other words:

- use Writing Outputs for external software workflows
- use the built-in `analysis` methods when you want ipyrad2 to manage filtering and method inputs internally

## Where to Go Next

- Read [`seqex`](./seqex.md) or [`snpex`](./snpex.md) for the command that matches the kind of export you need.
- Return to [Assemble](../assembly/assemble.md) if you are not yet at the HDF5 stage.
- Read the [Analysis Guide](../analyses/index.md) if you want to stay inside ipyrad2 instead of writing external-format files.
- See [Recipes](../recipes/index.md) for worked examples once that section is filled in. For now, Recipes is a visible TODO section where export and downstream workflow examples will live.
