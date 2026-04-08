# demux

`ipyrad2 demux` splits pooled reads into per-sample FASTQ files by barcode or index. If your sequencing provider already delivered sample-specific FASTQ files, this step is unnecessary and you can proceed to [trim](./trim.md).

## Overview
Demultiplexing is an optional step in ipyrad2 used to separate reads among samples that were pooled on a sequencing lane. We support demultiplexing using inline (internal) barcodes of fixed or variable length, including single-inline designs (barcode on one read) and dual-inline designs (barcodes on both reads), and we also support demultiplexing using external Illumina i7 index reads. Barcode mismatches can be tolerated up to an “off-by-n” limit, but only while the off-by-n barcode set remains unique (i.e., increasing n is disallowed once distinct barcodes become identical under the off-by-n expansion, typically at n ≥ 2).

Restriction cutsite motifs at the 5′ or 3′ end of single- or paired-end reads are auto-detected using fast k-mer analysis. Inline barcodes are parsed relative to the detected motif, and paired reads are evaluated jointly to improve sorting accuracy. Users can override the auto-detected motifs and will receive a warning if user-specified motifs disagree with the dominant patterns in the data. Reads that do not perfectly match one of the enumerated keys defined by the off-by-n map are assigned to an unassigned category. Demultiplexed reads are written as gzipped per-sample FASTQ files (R1 or R1/R2), and a summary statistics report counts reads assigned per sample and key.

## Prerequisites

- Pooled raw FASTQ files.
- A barcode or index table.
- A clear understanding of whether your run uses inline barcodes or `i7` indices.

## Inputs

### FASTQ inputs

Use `-d/--fastqs` with one or more raw FASTQ paths or shell-expanded globs.

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv
```

Multiple raw input files can be listed in one run. That is useful when a sample sheet should be applied across multiple lanes or technical replicate runs.

### Barcode table

The barcode table is whitespace-delimited, has no header, and should use:

```output
sample  barcode1  barcode2
```

`barcode2` is optional. If it is present, ipyrad2 treats the table as combinatorial inline barcodes and requires paired-end reads. If barcodes are only present on read1 the second barcode column can be absent.

For `--i7` demux, only `barcode1` is used. If extra barcode columns are present, ipyrad2 ignores them in that mode.

If duplicate sample names appear in the barcode file, ipyrad2 treats them as technical replicates. You can either merge them with `-M/--merge-technical-replicates` or let ipyrad2 write them as separate samples with `-technical-replicate-N` appended to the names.

## Command Patterns

The simplest inline-barcode run is:

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv -o DEMUX/
```

### Core inputs

- `-d, --fastqs`: pooled FASTQ inputs
- `-b, --barcodes`: barcode/index table
- `-o, --out`: output directory, default `./DEMUX`
- `-f, --force`: overwrite demux outputs from the current run

### Demultiplexing mode

- `--i7`: demultiplex by i7 index instead of inline barcodes
- `-m, --max-mismatch`: allow barcode/index mismatches, from `0` to `2`
- `-M, --merge-technical-replicates`: merge replicated sample names into one output sample

Use mismatches conservatively. Allowing mismatches can rescue reads from low-quality barcode positions, but it also increases the chance of ambiguous assignments when barcode sequences are too similar.

### Cutsite motifs

These options matter for inline demux:

- `-e1, --cutsite-1`: explicit R1 cutsite motif or motifs
- `-e2, --cutsite-2`: explicit R2 cutsite motif or motifs
- `-E, --disable-infer-cutsite-motifs`: skip motif inference
- `--max-reads-kmer`: reads sampled for motif inference, default `100000`

If you do not supply cutsite motifs, ipyrad2 tries to infer them from barcoded reads. That works well for standard libraries, but explicit motifs are better when:

- you already know the enzyme remnants
- you have a multi-enzyme design such as 3RAD
- inference is unstable on a small pilot dataset

### Performance and testing

- `-c, --cores`: maximum parallel workers, default `4`
- `-k, --chunksize`: reads per write batch, default `10000000`
- `-x, --max_reads`: stop after N reads per file for a test run
- `--pigz`: use `pigz` for final FASTQ compression

`--pigz` leads to big speed improvements, but larger disk usage temporarily before the FASTQs are ultimately compressed.

### Logging

- `-l, --log-level`: logging verbosity

At normal verbosity, ipyrad2 reports the chosen demux mode, inferred or manual cutsite motifs, technical replicate handling, and where the demux stats file was written.

## Outputs

`demux` writes one FASTQ per sample mate:

- `SAMPLE_R1.fastq.gz`
- `SAMPLE_R2.fastq.gz` for paired-end data

When `--pigz` is used, temporary plain `.fastq` files may appear during compression, but the intended final outputs are still the per-sample demultiplexed FASTQs.

It also writes a numbered demux report:

- `ipyrad_demux_stats_N.txt`

That report includes:

- raw file statistics
- sample demux statistics
- barcode detection statistics
- restriction motif inference details for inline demux
- barcode-boundary ambiguity and collision details when relevant

If technical replicates are merged, the stats report includes both replicate-level information and merged-sample totals.

## Common Failures

### Barcode table cannot be parsed

The barcode file must be whitespace-delimited and sample names cannot contain spaces. If parsing fails, check the delimiter and remove spaces from sample names.

### Combinatorial barcodes with single-end data

If the barcode table contains both `barcode1` and `barcode2`, ipyrad2 expects paired-end reads. Using a combinatorial barcode table with single-end inputs is an error.

### Existing outputs block the run

If the destination directory already contains output files with the expected demux names, ipyrad2 stops unless `--force` is set.

### Ambiguous barcodes when mismatches are allowed

`--max-mismatch` can make distinct samples overlap in barcode space. When that happens, ipyrad2 warns that ambiguous barcodes will be assigned arbitrarily to the first matching sample. The safer fix is usually to lower the mismatch setting.

### Motif inference looks wrong

If inferred cutsite motifs do not match the library design, inline demux can fail or produce boundary collisions. Provide `--cutsite-1` and `--cutsite-2` explicitly when the enzyme remnants are known.

### Zero-read samples

ipyrad2 warns when a listed sample receives zero reads. That usually means the sample truly failed, the barcode table is wrong, or the demux mode does not match the run structure.

## Examples

### Basic inline-barcode demux

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv -o DEMUX/ -c 10
```

### Demultiplex by i7 index

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv --i7 -o DEMUX/
```

### Merge technical replicates across runs

```bash
ipyrad2 demux -d RUN1/*.fastq.gz RUN2/*.fastq.gz -b BARCODES.tsv -M -o DEMUX/
```

### Set cutsite motifs manually

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv -e1 TGCAG -e2 CGG -o DEMUX/
```

### Test on a limited number of reads

```bash
ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.tsv -x 500000 -o DEMUX_TEST/
```

## Related Pages

- [Quick Guide](./index.md)
- [trim](./trim.md)
- [Files and Data Types](../getting-started/files-and-data-types.md)
