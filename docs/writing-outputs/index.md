# Writing Outputs

Writing Outputs is the export layer for workflows that need files for software outside ipyrad2. The commands in this section are `ipyrad2 analysis wex`, `ipyrad2 analysis lex`, and `ipyrad2 analysis snpex`. Together they let you turn assembled HDF5 data into sequence alignments, locus alignments, or filtered SNP matrices that match the needs of a particular downstream tool.

This design is especially useful for RAD-seq data because missing data usually needs to be handled explicitly rather than hidden. Different downstream methods tolerate different patterns of missing samples, missing sites, linked SNPs, and genotype imputation. The export layer makes those choices easy to control for each target analysis instead of forcing one global preprocessing decision on every output.

## Why This Layer Exists

RAD-seq datasets often need different filtering choices depending on the software you plan to run next. In practice that means deciding, for each export:

- which samples to keep or exclude
- whether to include `assembly_reference_sequence`
- how much missing data is acceptable per sample
- how much site or SNP coverage is required
- whether population-aware minimum coverage rules should apply
- whether linked SNPs should be subsampled to one SNP per RAD locus
- whether missing genotypes should stay missing or be imputed for a specific file format

That is the main point of `wex`, `lex`, and `snpex`: they let you tailor the exported dataset to the assumptions of a particular external program.

## The Three Export Tools

- [`wex`](./wex.md): writes one alignment from selected genomic windows in an assembly HDF5 file.
- [`lex`](./lex.md): writes sampled whole-locus alignments from an assembly HDF5 file.
- [`snpex`](./snpex.md): writes filtered SNP matrices from an SNP-capable HDF5 file, with optional PLINK, phylogenetic, TreeMix, and EEMS exports.

By combining these three tools, their filtering controls, and their output-format options, ipyrad2 can stage data for most common downstream sequence- and SNP-based software.

## Comparison

| Command | Input data | Output unit | Main filtering controls | Missing-data handling | Imputation availability | Typical downstream use |
| --- | --- | --- | --- | --- | --- | --- |
| `analysis wex` | assembly HDF5 | one alignment from selected windows | windows, `-m`, `-r`, `-e`, `-R`, `imap`, `minmap` | explicit sample and site filtering before export | none | sequence-based tools that need one alignment over chosen regions |
| `analysis lex` | assembly HDF5 | many delimited locus alignments | windows, number of loci, minimum locus length, `-m`, `-r`, `-e`, `-R`, `imap`, `minmap` | explicit sample and site filtering before locus export | none | multilocus sequence analyses and locus-based phylogenetic workflows |
| `analysis snpex` | SNP-capable HDF5 | filtered SNP matrices plus optional PLINK, phylogenetic, TreeMix, and EEMS files | `-m`, `-r`, `-a`, `-e`, `-R`, `imap`, `minmap`, linked vs unlinked SNPs | explicit sample and SNP filtering before export | optional, and applied to every written output in that run | SNP-based external tools such as PLINK, TreeMix, EEMS, or SNP-alignment workflows |

## Missing Data, Linkage, and Imputation

The key idea is not simply to export data, but to export the right version of the data for the next tool.

For sequence exports, `wex` and `lex` let you control which samples and sites survive filtering before alignments are written. That makes missing data explicit at the alignment stage instead of leaving it to downstream programs to interpret blindly.

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

- Read [`wex`](./wex.md), [`lex`](./lex.md), or [`snpex`](./snpex.md) for the command that matches the kind of export you need.
- Return to [Assemble](../assembly/assemble.md) if you are not yet at the HDF5 stage.
- Read the [Analysis Guide](../analyses/index.md) if you want to stay inside ipyrad2 instead of writing external-format files.
- See [Recipes](../recipes/index.md) for worked examples once that section is filled in. For now, Recipes is a visible TODO section where export and downstream workflow examples will live.
