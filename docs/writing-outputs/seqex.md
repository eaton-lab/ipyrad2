# seqex

<code>ipyrad2 seqex</code> exports sequence alignments from an assembly
HDF5 file. It applies locus and site coverage filters independently to
every selected locus and can write accepted loci as one multi-locus file,
one concatenated matrix, or separate files for each locus.

## Window selection and clipping

The form of <code>-w/--windows</code> determines whether loci remain complete
or are clipped automatically:

- Omitting <code>-w</code> selects complete loci from every scaffold.
- Scaffold names and regular expressions select complete overlapping loci.
- <code>scaffold:start-end</code> regions clip loci to the exact 1-based,
  inclusive interval boundaries.
- BED files clip loci to each 0-based, half-open BED interval after coordinate
  conversion.

When one coordinate interval crosses several loci, each locus intersection is
an independent unit. Two disjoint intervals within the same source locus also
produce two independent clipped fragments. These units can be filtered,
sampled with <code>-N</code>, concatenated, or written separately like complete
loci.

## Filtering order

For every complete locus or clipped fragment selected by
<code>-w/--windows</code>, seqex:

1. Counts a sample as present when it has at least one called base in the
   raw locus.
2. Rejects the locus unless sample presence satisfies <code>-m</code>, or
   each population satisfies <code>--minmap</code>.
3. Removes sites that fail those same coverage thresholds.
4. Removes samples whose missing fraction exceeds <code>-r</code>, or the
   effective group override from <code>--rmap</code>.
5. Applies <code>-L/--min-length</code> to the retained site count.
6. Randomly retains at most <code>-N/--max-loci</code> accepted loci, when
   requested.

Both <code>-L</code> and <code>-N</code> are disabled by default. A random
seed supplied with <code>-s</code> makes sampling reproducible. Randomly
selected loci are written in their original genomic order.

## Group coverage and missingness

IMAP is a sample-selection and grouping input, not merely a label annotation.
Only mapped samples are selected (subject to <code>--exclude</code>). A sample
must belong to one group; overlapping sample names or glob patterns across
groups are rejected.

Coverage policy is resolved independently of missingness:

| Inputs | Effective locus/site coverage rule |
|---|---|
| No IMAP | At least <code>-m</code> samples in total. |
| IMAP, without a coverage map | At least **one sample per group**, instead of <code>-m</code>. A warning reports this implicit rule. |
| IMAP and <code>--minmap</code> | The specified minimum for every group; overrides <code>-m</code>. All IMAP groups must be listed. |

A supported pops-format IMAP may supply the coverage minima itself; an explicit
<code>--minmap</code> takes precedence. To use IMAP only for group-specific
missingness, still provide an explicit minmap that states your intended coverage
policy. In particular, adding IMAP plus rmap does **not** preserve a scalar
<code>-m 10</code> rule automatically. Ten ingroup calls and zero required outgroup
calls is a different rule from ten calls across all samples.

<code>--rmap</code> (alias <code>--max-sample-missing-map</code>) is an optional
two-column, whitespace-separated map from IMAP group to maximum missing fraction:

~~~text
# group   maximum missing fraction
outgroup  1.0
~~~

It requires <code>--imap</code>. Listed groups override <code>-r</code>; groups
omitted from rmap inherit <code>-r</code>. Thus <code>-r 0.9</code> with the map
above allows at most 90% missingness in other groups, while outgroups have no
missingness-based exclusion. A threshold of 1.0 can retain an all-missing row;
check the reported called-base count and retention reason before inference.
No missing bases are imputed. N is missing, gaps are normalized to N, and IUPAC
ambiguity codes are counted as observed. Threshold comparison is inclusive.

Rmap accepts blank lines and <code>#</code> comments. Duplicate group entries,
unknown groups, malformed rows, and nonfinite or out-of-range fractions are
errors; values must lie in [0, 1]. The Python API also accepts
<code>rmap={"outgroup": 1.0}</code> with an IMAP dictionary or file. Without rmap,
existing scalar <code>-r</code> filtering is unchanged.

The same per-sample threshold is applied within every locus and, with
<code>-C</code>, across the concatenated matrix. It does not change coverage,
sample selection, read depth, or the missingness denominator. Split and
multi-locus outputs apply only the per-locus stage; their reported overall
sample missingness is descriptive and is not an additional global filter.

### Filter-order note

Coverage is evaluated **before** missingness-based sample removal. As in the
existing scalar-r workflow, subsequent row removal can leave sites below their
original coverage requirement. Rmap does not silently change that ordering or
iteratively filter the final matrix. Workflows requiring simultaneous final
coverage and missingness guarantees must independently recheck them. Statistics
record this coverage stage explicitly. A final-coverage enforcement change is
separate from group-specific missingness.

## Output layouts

With neither <code>-C</code> nor <code>-X</code>, one file contains
independent locus records:

- PHYLIP writes consecutive alignment matrices.
- NEXUS writes one named DATA block per locus.
- FASTA writes unique identifiers in the form <code>sample|locus</code>.

Use <code>-C/--concatenate</code> to append loci into one matrix. Samples
omitted from an individual locus receive an N block, and <code>-r</code> is
then applied once more across the complete matrix, using group-specific
<code>--rmap</code> overrides when provided. Concatenated boundaries
are recorded in the stats file.

Use <code>-X/--split</code> to write one file per accepted locus.
<code>-C</code> and <code>-X</code> are mutually exclusive, and split
output cannot be written to stdout.

With an IMAP, <code>--append-population</code> changes output names to
<code>population^sample</code>. Its short form is <code>-a</code>. This
ordering is compatible with BPP, which interprets the value after the caret
as the individual identifier.

## Parallel filtering

Use <code>-c/--cores</code> to filter batches of loci in parallel. The
default is one core. Worker processes read and filter bounded HDF5 batches;
the parent process restores genomic order, applies seeded locus sampling,
and writes output. Consequently, the same <code>-s</code> seed selects the
same loci regardless of the core count.

## Statistics

The stats report includes the total sites and sequence characters written,
the number of non-missing bases, and full-matrix non-missing occupancy.
Full-matrix occupancy treats a sample omitted from a retained locus as
missing. A per-sample table reports population, final output status, loci
written, loci removed by <code>-r</code>, non-missing bases, and occupancy.
Each sample row additionally reports its missing fraction, effective maximum
missing fraction, threshold source (<code>-r</code> or <code>rmap:GROUP</code>),
and final retention/removal reason. Matrix missingness uses all written locus
columns as the denominator, including N padding from loci where the sample was
removed. Reasons distinguish a global missingness failure, absence from all
retained loci after per-locus missingness filtering, absence of data in written
loci, and a retained row with no observed bases. The existing
<code>loci_dropped_by_r</code> counter includes removals by either the scalar
threshold or rmap and refers to the selected, written loci.

The summary records the requested <code>-m</code>, resolved IMAP and effective
minmap, and whether coverage came from scalar-m, provided minima, or the
implicit IMAP default. It also records default-r, group overrides, applicable
missingness stages, and the missing-character definition. These settings are
reported in both text and JSON, including the resolved map contents so results
remain interpretable if the original map file later changes.

The summary records the clipping mode and the original values passed to
<code>-w/--windows</code>, without expanding regexes or BED files. Each
written-locus row records the source locus, resolved selected window, output
coordinates, and whether the sequence was clipped. The same structured fields
are available in the JSON report.

Completion messages describe whether loci were written as independent
records, concatenated into one alignment, or split among separate files.
The stats-report path is logged separately. Samples removed by
<code>-r</code> are also reported.

## Examples

Write all accepted loci to one multi-locus PHYLIP file:

~~~bash
ipyrad2 seqex -d assembly.hdf5 -o output-seqex -m 4
~~~

Select a reproducible sample of 500 loci at least 150 retained sites long:

~~~bash
ipyrad2 seqex -d assembly.hdf5 -o output-seqex -N 500 -s 123 -L 150 -c 4 -O nex
~~~

Write a concatenated FASTA matrix with population-appended names:

~~~bash
ipyrad2 seqex -d assembly.hdf5 -i imap.tsv -g minmap.tsv -a -C -O fa
~~~

Write one PHYLIP file per locus overlapping two scaffolds:

~~~bash
ipyrad2 seqex -d assembly.hdf5 -w Chr01 Chr02 -X
~~~

Write one concatenated alignment clipped to exact BED intervals:

~~~bash
ipyrad2 seqex -d assembly.hdf5 -w windows.bed -C -O phy
~~~

Retain sparse outgroups while requiring ten ingroup samples per locus/site:

~~~text
# imap.tsv (sample or glob, then group)
focal-*    ingroup
external-* outgroup
~~~

~~~text
# minmap.tsv
ingroup  10
outgroup 0
~~~

~~~text
# missingness.tsv (ingroup inherits -r)
outgroup 1.0
~~~

~~~bash
ipyrad2 seqex -d assembly.hdf5 -o rooted-export -C \
  --imap imap.tsv --minmap minmap.tsv -r 0.9 --rmap missingness.tsv
~~~

This example requires ten **ingroup** calls, not ten calls in total. Group names
are user-defined filtering strata; they do not authenticate species, population,
or outgroup identity. Selecting thresholds is separate from assessing whether
sparse outgroups share enough informative sequence to support a stable root.

## See also

- [Writing Outputs](./index.md)
