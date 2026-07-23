# denovo

`ipyrad2 denovo` is an optional step used to build a pseudoreference genome FASTA from trimmed sample FASTQs. It is used when you do not have a suitable external reference genome, or your samples are highly divergent and may be biased by using a single reference.

In the normal assembly workflow, `denovo` sits after [`trim`](./trim.md) and before [`map`](./map.md). Its main output is a pseudoreference FASTA that you can map reads against. It is not the final assembled dataset, and it does not replace the later [`assemble`](./assemble.md) step.

![ipyrad2 assembly workflow from input reads to assembled outputs](../images/Fig1-assembly.png){ width="100%" }

## When to Use

Use `denovo` when: no suitable external reference genome exists; the available reference is too divergent to use confidently for read mapping; or you want to build a shared pseudoreference from the data themselves before mapping and assembly.

Do not use `denovo` when you already have a trusted reference FASTA that is appropriate for your samples. In that case, skip directly from [`trim`](./trim.md) to [`map`](./map.md).


## Overview

Reads are first clustered within samples at a high threshold (``--similarity-within``) to dereplicate and group reads representing alleles at the same locus into a consensus sequence. The consensus sequences are then clustered across samples at a lower threshold (``--similarity-across``) to group homologous loci, putatively including orthologs and paralogs. Using a similarity graph among samples, we apply a graph-splitting algorithm to separate duplicated components. The final set of graph components includes sets of sequences that include at most one sequence per sample, all of which are within the across-sample similarity threshold of each other. An alignment is performed on each set to generate a consensus that will serve as a locus in the pseudoreference. (See further details below).

## Prerequisites

Trimmed sample-level FASTQ or FASTQ.gz files, usually from [`trim`](./trim.md). All inputs in one run must be consistently single-end or consistently paired-end. (Note: you do not need to use all of your samples during this step, e.g., we recommend generally using 10-20 samples. You can still combine SE and PE data in the map stage using the pseudoreference assembled in this stage.)

## Command Patterns

The smallest useful run is:

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o output-denovo
```

That tells `ipyrad2` to parse sample names from the FASTQ filenames, run the denovo clustering workflow, and write a curated pseudoreference output set to `output-denovo/`.

### Core Inputs

- `-d, --fastqs`: one or more FASTQ paths or glob patterns
- `-o, --out`: output directory, default `./output-denovo`
- `-f, --force`: overwrite denovo outputs created by this command
- `-b, --allow-reverse-complement`: cluster both strands rather than plus strand only

### Clustering and Consensus

- `-s, --within-similarity`: within-sample clustering threshold, default `0.95`
- `-S, --across-similarity`: across-sample clustering threshold, default `0.85`
- `-m, --min-derep-size`: minimum duplicate count retained during dereplication, default `5`
- `-i, --min-length`: minimum retained sequence length after merge or join, default `35`
- `-g, --min-merge-overlap`: minimum overlap required to merge paired reads, default `20`
- `-e, --max-merge-diffs`: maximum mismatches allowed in the merged region, default `4`
- `--no-alignment`: skip MAFFT in the final locus step and use the longest stripped sequence per locus

### Sample Naming

- `-dx, --delim-str`: delimiter used when parsing sample names
- `-di, --delim-idx`: index of the retained delimiter-split token

### Runtime

- `-c, --cores`: maximum total cores to use, default `6`
- `-t, --threads`: threads per `vsearch` job, default `3`
- `--imap`: optional IMAP file used to choose one representative sample per group for denovo
- `--use-all-samples`: disable automatic downselection and use every parsed sample
- `--keep-intermediates`: retain the internal denovo working directory instead of cleaning it on success


## Workflow

1. parse sample names from input trimmed FASTQ (single or paired-end)
2a. (if paired-end): merge overlapping pairs and join unmerged pairs while recording spacer position
2b. dereplicate and cluster within samples to collapse into unique loci
3. write sample consensus fastas
4. cluster across-samples using vsearch all-by-all search among spacer-stripped consensus fastas
5. build graph of connected components
6. reconcile duplicated components that are below within-sample clustering threshold
7. split graph clusters into components using a maximum spanning forests algorithm
8. align sequences within each component using mafft and store the consensus as a pseudoreference locus
9. write pseudoreference loci to FASTA and write stats summary of graph splitting


## Performance

Two key parameter settings most significantly impact the runtime of the `denovo` step: `--imap` and `--min-dereplication-size`. Assembling a pseudoreference does not require every read from every sample. These options allow subselecting which samples and reads will be used to construct it.

### Dereplication

During this step, reads are dereplicated within each sample to retain only one copy of
each unique sequence, and reads that occur fewer than (`--min_dereplication_size`)
are discarded. This is an important heuristic. We find that a min dereplication size
of 5 typically performs well. For very low depth datasets you may want to reduce this.

### Sample selection

The goal of psuedoreference construction is to represent every ortholog in the dataset
by detecting it from at least one sample during this step. Analyzing and comparing
multiple samples increases the chances you will find more orthologs, and detect paralogs
that can be separated. However, there are usually diminishing returns from adding more and
more samples for this goal. Therefore, we recommend generally using ~10 of the most phylogenetically diverse, and high depth, samples in your dataset when running `denovo`.
Adding more samples can improve the result, but will increase runtimes non-linearly.

The current default auto-subsampling behavior is:

- keep all parsed samples when there are 10 or fewer
- when there are more than 10 samples, rank them by total input FASTQ size
- use the top 50% largest samples as the eligible pool
- if that pool has 10 samples, keep them all; if it has more than 10, randomly select 10
- if the eligible pool has fewer than 10 samples, fill the remaining slots from the next-largest samples until 10 are selected

You can override this in two ways:

- `--use-all-samples`: disable downselection and use every parsed input sample
- `--imap`: provide a two-column IMAP file and `denovo` will select one representative sample per group, choosing the largest available sample within each group

The IMAP path is the recommended way to steer denovo toward phylogenetically diverse representatives while still keeping the pseudoreference construction set relatively small. IMAP selection may still exceed 10 representatives, in which case `denovo` emits an advisory warning rather than truncating the set.


## Outputs and Intermediates

The primary result is the pseudoreference FASTA used for downstream mapping, named
`denovo_reference.fa`.

The `denovo` step also writes five curated outputs that can be optionally inspected.

- `denovo.stats.txt`: human-readable run summary with inputs, parameters, binaries, runtime settings, QC summaries, and output paths
- `denovo.loci.mapping.tsv`: mapping from final loci to the contributing sample-level consensus records
- `denovo.loci.stats.tsv`: per-locus summary table for the final denovo loci
- `denovo.sample_graph_summary.tsv`: per-sample burden of split, duplicated, and reconciled graph components
- `denovo.audit/`: compact audit files for duplicated, reconciled, or split components

`denovo.audit/` is intended for empirical diagnosis of difficult graph cases without requiring full intermediates to be preserved. In particular, `denovo.loci.mapping.tsv` and `denovo.loci.stats.tsv` expose reconciliation-related metadata such as the reconciliation mode and final output form, while the audit directory records per-component membership and summary information for duplicated components that were evaluated during graph refinement.

## Runtime and Performance Notes

- `--cores` controls total concurrency across the run.
- `--threads` applies to the `vsearch` stage.
- In the current implementation, the final MAFFT stage chooses its own worker scheduling internally from `--cores` and the number of loci.
- `--no-alignment` is usually much faster on large datasets because it skips MAFFT entirely in the last stage.

## Common Failures and Interpretation Notes

### Mixed single-end and paired-end inputs

One run must be consistently SE or consistently PE. Mixed layouts are rejected.

### Invalid clustering or merge parameters

`denovo` validates the similarity thresholds and length or overlap settings before launching work. For example:

- similarity thresholds must be greater than `0` and less than or equal to `1`
- dereplication size and minimum lengths must be at least `1`
- maximum merge differences must be at least `0`

### Existing outputs already exist

If curated denovo outputs already exist in the output directory, `denovo` stops unless you provide `--force`.

### Denovo output is not the final project output

`denovo_reference.fa` is a pseudoreference for mapping. It is not the final assembled HDF5, VCF, or loci dataset. The normal next step is [`map`](./map.md), followed by [`assemble`](./assemble.md).


## Examples

### Basic denovo pseudoreference build

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o output-denovo
```

### Use custom similarity settings and more concurrency

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o OUT -s 0.95 -S 0.85 -c 12 -t 3
```

### Select the graph refinement algorithm explicitly

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o OUT
```

### Use delimiter-based sample naming

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o OUT -dx _R -di 1
```

### Keep intermediate files for inspection

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o OUT --keep-intermediates
```

### Skip final alignment and write longest-sequence representatives

```bash
ipyrad2 denovo -d TRIMMED/*.fastq.gz -o OUT --no-alignment
```

### Empirical example

See the [Empiricial SE denovo tutorial](tutorial-ped.md) for an example run on real data.

## Related Pages

- [Quick Guide](./index.md)
- [trim](./trim.md)
- [map](./map.md)
- [assemble](./assemble.md)
