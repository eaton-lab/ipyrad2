# map

## Summary

`ipyrad2 map` aligns sample FASTQ files to a reference or denovo pseudoreference and writes coordinate-sorted, indexed BAM files. It uses `bwa-mem2` for alignment and `samtools` for filtering, sorting, indexing, duplicate removal, and stats reporting.

In the normal assembly workflow, `map` comes after [`trim`](./trim.md) and before [`Assemble`](./assemble.md). Its main job is to convert trimmed reads into final BAMs that are ready for locus assembly.

## When to Use

Use `map` when you have sample-level FASTQ files and a reference sequence to align against.

That reference can be:

- an external reference genome
- a denovo pseudoreference from [`denovo`](./denovo.md)

If you already have mapped BAM files, you do not need this step. Start directly at [`Assemble`](./assemble.md).

## Prerequisites

- Sample-level FASTQ files, usually from [`trim`](./trim.md)
- A reference or pseudoreference FASTA supplied with `-r/--reference`
- An activated ipyrad2 environment with executable `bwa-mem2` and `samtools`
- Read access to the FASTQ files and reference FASTA
- Write access to the reference directory if the reference still needs to be indexed

`map` supports plain FASTQ and `.gz`-compressed FASTQ input. `.bz2` is not supported.

## Inputs and Sample Grouping

Use `-d/--fastqs` with one or more FASTQ paths or shell-expanded globs. Inputs can be single-end or paired-end, but all samples in one run must be consistent. Mixed single-end and paired-end inputs are rejected.

Sample names are parsed from FASTQ filenames. If the default parsing is not right for your filenames, use:

- `-dx, --delim-str`: delimiter substring used to split the filename
- `-di, --delim-idx`: which side of the delimiter to keep, default `1`

After parsing, `map` strips a terminal `.trimmed` from the sample key when present. This keeps names from the normal `trim -> map` workflow canonical internally while still allowing later pipeline entry from externally named FASTQs.

You can also provide `-i/--imap` to subset, rename, or merge samples before mapping. The `imap` file is a whitespace-delimited two-column table:

```text
sample  group
```

This table can do three things:

- keep only a subset of parsed samples
- rename samples before BAM writing
- merge multiple parsed samples into one mapping target by concatenating their FASTQs

If some `imap` names do not match the canonical parsed sample names, ipyrad2 warns and skips them. If none match, the run stops.

For a worked example of explicit delimiter-based pairing and sample naming, see [Using -dx and -di to pair and name samples](../recipes/sample-name-parsing.md).

## Command Patterns

The smallest useful run is:

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o MAPPED/
```

That tells ipyrad2 to parse sample names from the FASTQ filenames, index the reference if needed, map reads with `bwa-mem2`, filter and sort alignments with `samtools`, and write final BAMs into `MAPPED/`.

### Core Inputs

- `-d, --fastqs`: one or more FASTQ paths or shell-expanded globs
- `-r, --reference`: reference FASTA or denovo pseudoreference FASTA
- `-o, --out`: output directory for BAMs and map stats, default `./MAPPED`

If the reference is not already indexed for `bwa-mem2`, `map` indexes it automatically before launching mapping jobs.

### Duplicate Removal

- `-m, --mark-dups-by-coords`: remove PCR duplicates by coordinates; intended for WGS-style data
- `-u, --mark-dups-by-umis`: remove PCR duplicates using UMI tags written by `ipyrad2 trim -U`

These modes are mutually exclusive, and duplicate removal is allowed only for paired-end data.

Coordinate-based duplicate removal is not the normal RAD setting. The mapper warns about that explicitly because RAD data naturally share start coordinates within loci. Use `-m` only when the run really contains WGS-style data, not ordinary RAD samples.

UMI-based duplicate removal is appropriate only when your reads were prepared with [`trim`](./trim.md) using `-U` so the i5/index2 value was stored in the read name as a UMI tag.

### Sample Naming and Grouping

- `-i, --imap`: subset, rename, or merge parsed samples using a sample-to-group table
- `-dx, --delim-str`: delimiter used to parse sample names from FASTQ filenames
- `-di, --delim-idx`: index of the retained token when splitting names

Use these when FASTQ filenames do not follow the default parser assumptions or when technical replicates should be merged before mapping.

### Performance and Overwrite

- `-c, --cores`: maximum total cores to use, default `6`
- `-t, --threads`: threads per mapping job, default `3`
- `-f, --force`: overwrite existing `.bam` and `.bam.csi` outputs for matching sample names

ipyrad2 runs up to `cores // threads` mapping jobs in parallel. `--threads` cannot exceed `--cores`.

Without `--force`, samples that already have BAM outputs are skipped instead of being remapped.

### Logging

- `-l, --log-level`: logging verbosity, default `INFO`

At normal verbosity, ipyrad2 reports parsed sample names, duplicate-removal warnings, reference indexing, mapping progress, and where the map stats file was written.

## Outputs and Stats

For each sample, `map` writes:

- `SAMPLE.trimmed.sorted.bam`
- `SAMPLE.trimmed.sorted.bam.csi`

During the run it also uses `OUTDIR/tmpdir/` for temporary sort, fixmate, and stats files, but those temporary files are cleaned up when mapping finishes.

For the run as a whole, it writes:

- `ipyrad_map_stats_N.txt`

The stats report is numbered so repeated runs in the same output directory do not overwrite older summaries. It reports final BAM retention and alignment-quality summaries rather than full aligner logs.

For single-end data, the report includes fields such as:

- reads processed
- reads filtered before BAM writing
- reads retained in the final BAM
- proportions retained
- counts of reads below reporting thresholds for MAPQ, soft clipping, and NM

For paired-end data, it also reports:

- pairs evaluated in the final BAM
- pairs with both mates retained
- singleton counts
- duplicate-removal totals when duplicate marking was enabled
- pair-level reporting summaries for MAPQ, soft clipping, NM, and insert size

The MAPQ, soft-clipping, NM, and TLEN thresholds in this report are reporting thresholds only. They are not additional mapper filters applied during BAM generation.

## Common Failures and Interpretation Notes

### No FASTQ inputs match

If the shell does not expand your `-d` pattern the way you expect, ipyrad2 may see no usable inputs. Check the glob itself first and make sure the files are visible from your current working directory.

### Sample names are parsed incorrectly

If mates are not grouped together correctly, or several files collapse into one sample unexpectedly, adjust `--delim-str` and `--delim-idx`.

### Unsupported compression

Only plain FASTQ and `.gz` FASTQ are supported. `.bz2` input raises an error before mapping begins.

### Mixed SE and PE inputs

If some parsed samples are paired-end and others are single-end, the run stops. That includes cases where an `imap` merge would mix single-end and paired-end inputs into one output sample.

### Reference indexing fails

`map` auto-indexes the reference with `bwa-mem2` when needed. If the reference path does not exist or its directory is not writable, indexing fails before mapping starts.

### Mapper dependencies are missing

Both `bwa-mem2` and `samtools` must exist in the active ipyrad2 environment and be executable.

### `threads` exceeds `cores`

`map` validates this before running. Increase `--cores`, reduce `--threads`, or both.

### Duplicate-removal mode is invalid

The two duplicate-removal modes cannot be selected together. Duplicate removal also cannot be used for single-end data.

### `imap` names do not match parsed samples

If some names in the `imap` file do not match canonical parsed sample names, ipyrad2 warns and skips them. If no names match at all, the run stops and the `imap` file or parsing arguments need to be fixed.

### Disk-space failures during sorting

The mapping pipeline creates temporary BAMs during sorting and duplicate removal. If the output filesystem fills up, `samtools sort` can fail with a disk-space error. ipyrad2 now surfaces those failures more cleanly, but the fix is still to free space or choose a different output location.

## Examples

### Basic reference-based mapping

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o MAPPED/
```

### Map against a denovo pseudoreference

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r DENOVO/pseudoref.fa -o MAPPED/
```

### Subset, rename, or merge samples with `imap`

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o MAPPED/ -i IMAP.tsv
```

### Coordinate-based duplicate removal for WGS-style data

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o MAPPED/ -m
```

### UMI-based duplicate removal after `trim -U`

```bash
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o MAPPED/ -u
```

## Related Pages

- [Quick Guide](./index.md)
- [trim](./trim.md)
- [Denovo](./denovo.md)
- [Assemble](./assemble.md)
