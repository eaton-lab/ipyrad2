# wex

## Summary

`ipyrad2 analysis wex` extracts one concatenated alignment from selected genomic windows in an assembly HDF5 file. If you do not specify any windows, `wex` concatenates sequence data across the full length of all scaffolds in genomic order. If you do specify windows, it concatenates only those selected scaffold regions into one alignment.

That makes `wex` the right export when you want one supermatrix for an external sequence-based tool, such as a phylogenetic tree search in `raxml-ng`.

## When to Use

Use `wex` when your next step needs one alignment rather than many separate loci.

Typical use cases include:

- one genome-wide concatenated alignment for a tree search
- one concatenated alignment per chromosome or scaffold
- one concatenated alignment from a custom set of regions or BED intervals
- one filtered alignment for an external program that expects FASTA, PHYLIP, or NEXUS input

If you want many locus files instead of one concatenated alignment, use [`lex`](./lex.md) instead.

## Prerequisites

- A sequence-capable assembly HDF5 file from `ipyrad2 assemble`
- Scaffold names that you can inspect with `--print-scaffold-table`
- A clear idea of whether you want:
  - all scaffolds
  - one scaffold at a time
  - explicit scaffold ranges
  - BED-selected windows

`wex` is for sequence HDF5 input, not SNP-only HDF5 input.

## Inputs and Window Selection

The real command is:

```bash
ipyrad2 analysis wex -d assembly.hdf5 [OPTIONS]
```

The main inputs are:

- `-d, --data`: assembly HDF5 file
- `-n, --name`: output prefix, default `alignment`
- `-o, --out`: output directory, default `output-wex`
- `-O, --out-format`: `phy`, `nex`, or `fa`, default `phy`

Window selection is controlled by `-w/--windows`.

`-w` accepts:

- scaffold names such as `Chr1`
- regexes such as `Chr[1-3]`
- region strings such as `Chr1:1-500000`
- a BED file path

If `-w` is omitted, `wex` selects the full length of all scaffolds. For region strings, coordinates are 1-based inclusive. BED input uses standard 0-based half-open coordinates and is converted internally.

Before selecting windows, you can print the scaffold table:

```bash
ipyrad2 analysis wex -d assembly.hdf5 --print-scaffold-table
```

That table contains:

- `scaffold_name`
- `scaffold_length`

Use it to decide which scaffold names to pass to `-w`. The table can be piped through short consumers without emitting a broken-pipe error:

```bash
ipyrad2 analysis wex -d assembly.hdf5 -P | head
```

## Filtering and Sample Selection

`wex` exposes the sequence filters that usually matter most for downstream alignments:

- `-m, --min-sample-coverage`: minimum number of samples with data required at a site, default `4`
- `-r, --max-sample-missing`: maximum missing-site frequency allowed before dropping a sample, default `1.0`
- `-e, --exclude`: exclude one or more samples by name
- `-R, --include-reference`: include `assembly_reference_sequence`
- `-i, --imap`: sample-to-population mapping file
- `-g, --minmap`: population-to-minimum-coverage file that overrides `-m`

These filters are applied before the alignment is written. That means the final concatenation can differ substantially depending on how you filter missing data, populations, and excluded samples.

If you use `-i`, then `-R` only works if `assembly_reference_sequence` is also assigned to an IMAP group.

## Output Control

- `-P, --print-scaffold-table`: print scaffold table to stdout and exit
- `-x, --stdout`: write the alignment to stdout instead of a file
- `-f, --force`: overwrite existing outputs with the same name
- `-l, --log-level`: logging verbosity

## Outputs and Stats

One `wex` run writes one alignment plus one stats file.

Depending on `-O`, the main alignment file is:

- `PREFIX.phy`
- `PREFIX.nex`
- `PREFIX.fa`

The stats file is always:

- `PREFIX.stats.txt`

The stats file starts with the `CMD:` invocation used for the run and records:

- samples before filtering
- sites in the selected windows before filtering
- variant sites in the selected windows before filtering
- samples after filtering
- sites in the selected windows after filtering
- variant sites after filtering
- input file
- output file
- selected windows
- IMAP / minmap context

If `--stdout` is used, the alignment is written to stdout, but the stats file is still written unless output is suppressed elsewhere.

## Common Failures and Interpretation Notes

### Wrong database type

`wex` expects an assembly HDF5 with sequence data. If the input is not a sequence-capable HDF5, window extraction will fail before alignment writing.

### Scaffold names do not match

If a scaffold or regex does not match anything in the scaffold table, `wex` raises an error. Use `--print-scaffold-table` first.

### Malformed regions

Region strings must look like `scaffold:start-end`. Start must be positive, and end must be greater than or equal to start.

### Regex cannot be combined with `:`

`wex` does not allow regex matching inside `scaffold:start-end` syntax. If you want several region strings, list them separately:

```bash
-w Chr1:1-1000 Chr2:1-1000
```

### Overlapping windows

Selected windows cannot overlap. This applies both to command-line windows and BED input.

### Selected windows contain no sequence data

If the selected windows do not overlap any assembled sequence blocks, `wex` stops with a zero-data error and you need to choose larger or different windows.

### Selected windows lose all data after filtering

If coverage filtering removes all sites, `wex` stops with a coverage-related zero-data error. This often means `-m`, `-g`, or sample filtering is too strict for the selected regions.

### No samples pass `max_sample_missing`

If all samples are dropped by `-r`, no alignment is written.

### Existing outputs already exist

If `PREFIX.phy`/`.nex`/`.fa` or `PREFIX.stats.txt` already exists, `wex` stops unless `--force` is set.

## Examples

### Genome-wide concatenation

The simplest whole-assembly concatenation omits `-w` entirely:

```bash
ipyrad2 analysis wex \
  -d OUT/assembly.hdf5 \
  -o WEX/ \
  -n genome \
  -O phy
```

This writes one genome-wide concatenated PHYLIP alignment:

```text
WEX/genome.phy
WEX/genome.stats.txt
```

### One concatenation per chromosome

First inspect the scaffold table:

```bash
ipyrad2 analysis wex -d OUT/assembly.hdf5 --print-scaffold-table
```

Then write one concatenated alignment for a single scaffold:

```bash
ipyrad2 analysis wex \
  -d OUT/assembly.hdf5 \
  -o WEX/by_chrom/ \
  -n Chr1 \
  -w Chr1 \
  -O phy
```

This produces:

```text
WEX/by_chrom/Chr1.phy
WEX/by_chrom/Chr1.stats.txt
```

To do this for every scaffold in the table:

```bash
ipyrad2 analysis wex -d OUT/assembly.hdf5 --print-scaffold-table > scaffold_table.tsv

awk -F '\t' 'NR>1 {print $2}' scaffold_table.tsv | while read chrom; do
  ipyrad2 analysis wex \
    -d OUT/assembly.hdf5 \
    -o WEX/by_chrom/ \
    -n "$chrom" \
    -w "$chrom" \
    -O phy
done
```

The `awk` command uses column 2 because the printed scaffold table includes a leading DataFrame index column, followed by `scaffold_name` and `scaffold_length`.

### Concatenate a custom region

```bash
ipyrad2 analysis wex \
  -d OUT/assembly.hdf5 \
  -o WEX/ \
  -n chr1_head \
  -w Chr1:1-500000 \
  -O fa
```

### Concatenate windows from a BED file

```bash
ipyrad2 analysis wex \
  -d OUT/assembly.hdf5 \
  -o WEX/ \
  -n selected_regions \
  -w windows.bed \
  -O nex
```

### Genome-wide tree inference in RAxML-NG

Once you have a concatenated PHYLIP alignment, you can pass it directly to `raxml-ng`:

```bash
raxml-ng \
  --search \
  --msa WEX/genome.phy \
  --model GTR+G \
  --prefix TREES/genome \
  --redo
```

This is a natural use of `wex`: create one filtered concatenation in ipyrad2, then use that exported alignment in an external tree-search program.

### Per-chromosome tree inference in RAxML-NG

If you wrote one `wex` alignment per chromosome, you can search one tree per file:

```bash
for msa in WEX/by_chrom/*.phy; do
  chrom=$(basename "$msa" .phy)
  raxml-ng \
    --search \
    --msa "$msa" \
    --model GTR+G \
    --prefix "TREES/$chrom" \
    --redo
done
```

This is often useful when you want one tree for each chromosome-scale concatenation rather than one tree from the whole assembly.

## Related Pages

- [Writing Outputs](./index.md)
- [lex](./lex.md)
- [snpex](./snpex.md)
- [Recipes](../recipes/index.md)
