# What Is ipyrad2?

ipyrad2 is a command-line toolkit for assembling and analyzing RAD-seq and related reduced-representation genomic
datasets (e.g., RAD, ddRAD, GBS, 2bRAD, 3RAD).

In practical terms, it helps you move from reads or mapped BAMs to assembled loci, variant calls, exported sequence
or SNP datasets, and downstream population-genetic analyses, all within one workflow.

## Major new features in ipyrad2

- Significant speed, memory, and disk-usage improvements relative to ipyrad1
- An atomized workflow makes it easier to handle very large datasets
- More logging makes the internal steps more transparent
- Automated kmer-based detection of restriction cutsite motifs
- Denovo assembly implements graph-based paralog splitting
- Variant calling now uses a more standardized bcftools pipeline
- Faster and easier options for filtering and writing outputs
- Supports assembly of RAD and WGS samples together in one dataset
- and more...


## The subcommand structure

ipyrad2 organizes the workflow into explicit command stages. Reads can be demultiplexed, trimmed, used to build a _denovo_ pseudoreference when needed, mapped to a reference, assembled into loci and variants, and then exported or analyzed directly. The top-level subcommands reflect that structure: `demux`, `trim`, `denovo`, `map`, `assemble`, and several analysis tools.

One of the main differences from the older ipyrad interface is that ipyrad2 centers the command line around named subcommands with clearer inputs, outputs, logs, and stats. The older ipyrad format, which uses a single command and a params file to organize the workflow, can still be used to implement the new pipeline, if preferred, through a companion tool made available as ``ipyrad2-classic``.

The assembled outputs of ipyrad2 assembly includes a VCF as well as a HDF5 database file. The latter can serve as a substrate for downstream analyses implemented in ipyrad2, or to generate filtered and formatted output files for external tools by using `seqex` and `snpex`.


## Support for RAD & WGS data

Increasingly, whole genome sequence (WGS) data is available for many organisms, and it is desirable to be able to combine this data with RAD-seq type samples. This is supported in ipyrad2, while retaining the benefits of _reduced genome representation_. Specifically, RAD-seq based samples are used to delimit locus regions (beds) along a reference genome, and for WGS samples, reads mapping to these beds are used to call variants and assemble loci.

This is one of the most important additions in ipyrad2.


## Where to go next

- Read [The ipyrad Ethos](./ethos.md) for the design goals and project philosophy.
- Read [Installation](./installation.md) if you need to set up the environment.
- Read [Files and Data Types](./files-and-data-types.md) for the current file model, HDF5 outputs, mixed RAD/WGS workflows, and supported workflow entry points.
- Read [Assembly](../assembly/index.md) for the end-to-end workflow.
- Read [trim](../assembly/trim.md), [map](../assembly/map.md), or the [Analysis Guide](../analyses/index.md) when you are ready for step-specific details.
