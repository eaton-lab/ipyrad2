# seqex

<code>ipyrad2 seqex</code> is an experimental exporter for complete,
delimited loci in an assembly HDF5 file. It applies the same locus and site
coverage requirements independently to every locus and can write accepted
loci as one multi-locus file, one concatenated matrix, or separate files.

The established [lex](./lex.md) and [wex](./wex.md) commands remain
available while the behavior and performance of seqex are evaluated.

## Filtering order

For every complete locus selected by <code>-w/--windows</code>, seqex:

1. Counts a sample as present when it has at least one called base in the
   raw locus.
2. Rejects the locus unless sample presence satisfies <code>-m</code>, or
   each population satisfies <code>--minmap</code>.
3. Removes sites that fail those same coverage thresholds.
4. Removes samples whose missing fraction exceeds <code>-r</code>.
5. Applies <code>-L/--min-length</code> to the retained site count.
6. Randomly retains at most <code>-N/--max-loci</code> accepted loci, when
   requested.

Both <code>-L</code> and <code>-N</code> are disabled by default. A random
seed supplied with <code>-s</code> makes sampling reproducible. Randomly
selected loci are written in their original genomic order.

Regions and BED intervals select whole overlapping loci. They do not clip
sequence at interval boundaries; use wex when exact clipping is required.

## Output layouts

With neither <code>-C</code> nor <code>-X</code>, one file contains
independent locus records:

- PHYLIP writes consecutive alignment matrices.
- NEXUS writes one named DATA block per locus.
- FASTA writes unique identifiers in the form <code>sample|locus</code>.

Use <code>-C/--concatenate</code> to append loci into one matrix. Samples
omitted from an individual locus receive an N block, and <code>-r</code> is
then applied once more across the complete matrix. Concatenated boundaries
are recorded in the stats file.

Use <code>-X/--split</code> to write one file per accepted locus.
<code>-C</code> and <code>-X</code> are mutually exclusive, and split
output cannot be written to stdout.

With an IMAP, <code>--append-population</code> changes output names to
<code>sample^population</code>. Its short form is <code>-a</code>.

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

## See also

- [Writing Outputs](./index.md)
- [lex](./lex.md)
- [wex](./wex.md)
