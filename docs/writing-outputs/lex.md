# lex

## Summary

`ipyrad2 lex` extracts complete, delimited loci from an assembly HDF5
database. It can write each selected locus separately or join accepted loci
into one site-wise concatenated alignment. Outputs use PHYLIP, NEXUS, or BPP
sequence formats.

This tool is especially useful for de novo RAD assemblies, where the HDF5
contains many independent loci. It also works with reference-based assemblies
when complete assembled loci, rather than clipped genomic windows, are the
desired export unit.

## When to Use

Use `lex` for:

- separate alignments of randomly sampled or region-selected loci
- a multi-locus BPP sequence file
- one concatenated alignment made from complete loci
- reproducible random locus selection
- population-aware sequence filtering

Use [`wex`](./wex.md) to export the literal sequence inside genomic windows.
Use [`snpex`](./snpex.md) for SNP matrices and SNP-specific formats.

## Prerequisites

- An assembly HDF5 produced by `ipyrad2 assemble`, version 2.0 or newer.
- Sequence and locus-map datasets in the HDF5. SNP-only databases are not
  valid `lex` inputs.
- Sample, scaffold, and population names matching those stored in the HDF5.

```bash
ipyrad2 lex -d assembly.hdf5 [OPTIONS]
```

## How Loci Are Selected

`lex` treats the delimited loci in the HDF5 locus map as indivisible
candidates. Windows determine which loci are eligible but do not clip their
sequences. A locus spanning `Chr1:90-160` selected by
`-w Chr1:100-120` is exported in full.

Selection proceeds as follows:

1. Find complete loci overlapping the selected scaffolds or windows.
2. Remove loci shorter than `--min-length`.
3. Randomize candidate order; `--random-seed` makes the order reproducible.
4. Apply site-coverage and sample-missingness filters to each locus.
5. Reject a filtered locus if it is now shorter than `--min-length`.
6. Continue until `--max-loci` loci pass or candidates are exhausted.

Thus, `--max-loci` is the maximum number of accepted loci, not the number
examined. If fewer pass, `lex` writes those that passed and logs a warning.

## Quick Examples

Write up to 100 loci as separate PHYLIP files:

```bash
ipyrad2 lex -d assembly.hdf5 -o output-lex -N 100 -L 150 -O phy
```

Write NEXUS files for complete loci overlapping two scaffolds:

```bash
ipyrad2 lex \
  -d assembly.hdf5 \
  -o output-lex \
  -n chromosome_loci \
  -w Chr01 Chr02 \
  -N 50 \
  -O nex
```

Select loci overlapping a region or a BED file:

```bash
ipyrad2 lex -d assembly.hdf5 -w Chr01:1-500000 -N 25
ipyrad2 lex -d assembly.hdf5 -w target_regions.bed -N 25
```

Write one reproducible concatenated alignment:

```bash
ipyrad2 lex \
  -d assembly.hdf5 \
  -o output-lex \
  -n concat \
  -N 100 \
  -s 123 \
  -C \
  -O phy
```

Apply population-aware coverage filters:

```bash
ipyrad2 lex \
  -d assembly.hdf5 \
  -o output-lex \
  -i imap.tsv \
  -g minmap.tsv \
  -N 100
```

Write a standard multi-locus BPP file:

```bash
ipyrad2 lex -d assembly.hdf5 -o output-lex -n bpp_data -N 1000 -O bpp
```

## Parameter Reference

### Core Inputs

| Option | Default | Effect |
| --- | --- | --- |
| `-d Path`, `--data Path` | Required | Assembly HDF5 containing sequences and a locus map. |
| `-n str`, `--name str` | `alignment` | Prefix used for alignment and stats files. |
| `-o Path`, `--out Path` | `output-lex` | Directory for output files. |
| `-O str`, `--out-format str` | `phy` | Format: `phy`, `nex`, or `bpp`. |

For example, `-n oak` produces `oak.stats.txt` and alignment names
beginning with `oak`.

### Locus Sampling

| Option | Default | Effect |
| --- | --- | --- |
| `-w [str ...]`, `--windows [str ...]` | All scaffolds | Restrict eligible loci to those overlapping selected scaffolds or regions. |
| `-N int`, `--max-loci int` | `100` | Maximum number of accepted loci to write. |
| `-s int`, `--random-seed int` | Random | Non-negative seed for reproducible candidate randomization. |
| `-L int`, `--min-length int` | `150` | Minimum locus length before filtering and retained length afterward. |

#### Window Syntax

`--windows` accepts:

- an exact scaffold name, such as `-w Chr01`
- a full-match regular expression, such as `-w 'Chr0[1-3]'`
- a 1-based inclusive region, such as `-w Chr01:1-500000`
- one BED file using standard 0-based, half-open coordinates

Extra BED columns are ignored. If `-w` is omitted, every scaffold is
eligible. A region must resolve to one scaffold, so its scaffold component
cannot be a regex matching several names. Selected windows cannot overlap.

Inspect scaffold names and lengths with:

```bash
ipyrad2 lex -d assembly.hdf5 --print-scaffold-table
```

### Filtering and Samples

Filters are applied independently to each locus. Site coverage is evaluated
first, then sample missingness is measured across retained sites.

| Option | Default | Effect |
| --- | --- | --- |
| `-m int`, `--min-sample-coverage int` | `4` | Retain a site when at least this many selected samples have data. Without IMAP, all samples form one group. |
| `-r float`, `--max-sample-missing float` | `1.0` | Drop a sample from a locus when its missing-site fraction exceeds this value. Values are constrained internally to `0.0` through `1.0`. |
| `-e [str ...]`, `--exclude [str ...]` | None | Exclude exact sample names. This takes precedence over IMAP and `--include-reference`. |
| `-R`, `--include-reference` | Off | Include `assembly_reference_sequence`, which is otherwise excluded unless IMAP contains it. |
| `-i Path`, `--imap Path` | None | Subset samples and assign populations using exact names or glob patterns. |
| `-g Path`, `--minmap Path` | None | Require an integer number of samples with data in every IMAP population at each retained site. |

Without IMAP, `-m` applies to an implicit group named `all`. With IMAP,
only mapped samples are retained. If IMAP is provided without minmap, each
population defaults to one required sample. Supplying `-g` replaces this
default and the global `-m` threshold with per-population requirements.

An IMAP file maps samples or globs to populations:

```text
sample_a    pop1
sample_b    pop1
sample_c*   pop2
```

A minmap file maps the same populations to integer minimums:

```text
pop1    2
pop2    1
```

The group sets must match exactly. `lex` also accepts the legacy combined
format with minimums on a final comment line:

```text
sample_a    pop1
sample_b    pop1
sample_c*   pop2
# pop1:2 pop2:1
```

When `-R` and `-i` are combined, `assembly_reference_sequence` must
appear in an IMAP group. Explicit exclusion still takes precedence.

### Output Control and Logging

| Option | Default | Effect |
| --- | --- | --- |
| `-C`, `--concatenate` | Off | Append accepted loci end to end and write one locus in the selected format. |
| `-P`, `--print-scaffold-table` | Off | Print scaffold names and lengths as TSV, then exit without extracting loci. |
| `-x`, `--stdout` | Off | Write sequences to stdout instead of alignment files. The stats file is still written. |
| `-f`, `--force` | Off | Overwrite existing alignment and stats files. |
| `-l str`, `--log-level str` | `INFO` | Logging level: `TRACE`, `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `-h`, `--help` | Off | Print command help and exit. |

Without concatenation, `--stdout` can emit several consecutive alignment
records. With concatenation, stdout contains one alignment:

```bash
ipyrad2 lex -d assembly.hdf5 -N 25 -s 123 -C -O phy -x | head
```

Logs are written separately from sequence stdout, so the sequence stream can
be piped to another program.

## Output Formats and Names

| Format and mode | Output with `-n alignment` |
| --- | --- |
| Separate PHYLIP | `alignment.Chr01_1-150.phy`, ... |
| Separate NEXUS | `alignment.Chr01_1-150.nex`, ... |
| Multi-locus BPP | `alignment.phy` with multiple locus records |
| Concatenated PHYLIP | `alignment.phy` |
| Concatenated NEXUS | `alignment.nex` |
| Concatenated BPP | `alignment.bpp` with one locus record |
| Stats | `alignment.stats.txt` |

Coordinates remain colon-delimited in stats, but individual filenames use
`_start-end` so they are portable across filesystems.

### Concatenated Output

`-C/--concatenate` preserves PHYLIP, NEXUS, or BPP syntax while turning all
accepted loci into one locus. Loci are appended in selection order.

Sample filtering still occurs per locus. The result includes the union of
samples that passed at least one accepted locus, in original HDF5 order. If a
retained sample did not pass a locus, an `N` block equal to that filtered
locus length is inserted. Samples that fail every accepted locus are omitted.

BPP accepts `--concatenate`, but `lex` warns that BPP is intended for
multi-locus analyses and one concatenated locus is likely unsuitable.

## Stdout Examples

### Sequence Output

Given two four-site loci and three samples:

```bash
ipyrad2 lex \
  -d assembly.hdf5 \
  -o LEX \
  -n demo \
  -w chr1 \
  -N 2 \
  -s 17 \
  -L 4 \
  -m 1 \
  -r 0.5 \
  -C \
  -O phy \
  -x
```

The PHYLIP alignment written to stdout is:

```text
3 8
s1      AAAACCCC
s2      AAAANNNN
s3      NNNNCCCC
```

Here `s2` failed the second locus and `s3` failed the first, so each has
one four-character `N` block. Informational logs go to stderr and are not
part of this stream.

### Scaffold Table

```bash
ipyrad2 lex -d assembly.hdf5 -P
```

```text
	scaffold_name	scaffold_length
0	chr1	8
1	chr2	12500
```

The leading field is the table row index. Short consumers are supported:

```bash
ipyrad2 lex -d assembly.hdf5 -P | head
```

## Stats File

Every extraction writes one `NAME.stats.txt`, including runs using stdout.
The `Summary` section records settings and run-level counts. The
`Accepted loci` table records coordinates, destinations, and before/after
filtering counts for every written locus.

For the stdout example, `LEX/demo.stats.txt` contains data like:

```text
Summary
-------
tool: lex
name: demo
infile: assembly.hdf5
out_format: phy
concatenate: True
random_seed: 17
outfile: STDOUT
concatenated_sites: 8
nloci_requested: 2
nloci_written: 2
min_length_requested: 4
eligible_loci_before_filtering: 2
loci_rejected_after_filtering: 0
windows: ['chr1:5-8']
imap: {'all': ['s1', 's2', 's3']}
min_sample_coverage_filter: {'all': 1}
max_sample_missing_filter: 0.5

Accepted loci
-------------
locus_index  locus_name  scaffold  start  end  outfile  nsamples_before_filtering  nsites_in_windows_before_filtering  nvariants_in_windows_before_filtering  nsamples_after_filtering  nsites_in_windows_after_filtering  nvariants_in_windows_after_filtering
1            chr1:1-4    chr1      1      4    STDOUT   3                          4                                   0                                      2                         4                                  0
2            chr1:5-8    chr1      5      8    STDOUT   3                          4                                   0                                      2                         4                                  0
```

For separate PHYLIP or NEXUS files, summary `outfile` is `multiple`, and
each row names its file. Multi-locus BPP and concatenated runs name their
shared file. `concatenated_sites` is zero without `--concatenate`.

## Common Failures

- **Wrong database:** `lex` requires a sequence-capable assembly HDF5,
  version 2.0 or newer.
- **No matching scaffold:** use `-P` to inspect names. Regions require valid
  positive coordinates in `scaffold:start-end` form.
- **Overlapping windows:** command-line and BED windows cannot overlap.
- **No loci meet the length requirement:** select other windows, reduce
  `-L`, or relax `-m` or `-g`.
- **No sites or samples remain:** relax the coverage or missingness filters.
- **IMAP/minmap mismatch:** samples must match HDF5 names or globs, cannot
  appear in several groups, and minmap groups must equal IMAP groups.
- **Reference not assigned:** with `-R -i`, add
  `assembly_reference_sequence` to one IMAP group.
- **Invalid numeric arguments:** use a non-negative `-s` and positive values
  for `-N` and `-L`.
- **Output exists:** choose another `-n` or `-o`, or use `-f`.

## Related Pages

- [Writing Outputs](./index.md)
- [wex](./wex.md)
- [snpex](./snpex.md)
- [Files and Data Types](../getting-started/files-and-data-types.md)
