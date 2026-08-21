
# seqex recipes

In the examples below I either show a generic name for an input file (ASSEMBLY.hdf5) when demonstrating options
where I do not show example results, or I use `SRP021469/OUT/assembly.hdf5` as the input file when using an
example database file from the [empirical SE denovo assembly tutorial]().


### windows (reference)
The `-w/--windows` argument can be used to specify one or more specific scaffolds/chromosomes by name
from which you want to extract filtered loci to concatenate.


```bash
ipyrad2 seqex \
  -d ASSEMBLY.hdf \
  -o SEQEX \
  -w CHROM_01 \
  -C
```

You can use regular expression to select multiple scaffolds
```bash
ipyrad2 seqex \
  -d ASSEMBLY.hdf \
  -o SEQEX \
  -w CHROM_0[1-5] \
  -C
```

Or to select multiple scaffolds you can list them manually

```bash
ipyrad2 seqex \
  -d ASSEMBLY.hdf \
  -o SEQEX \
  -w CHROM_01 CHROM_02 CHROM_03 CHROM_04 CHROM_05 \
  -C
```

Or you can enter a BED formatted file

```bash
ipyrad2 seqex \
  -d ASSEMBLY.hdf \
  -o SEQEX \
  -w WINDOWS.bed \
  -C
```

### min-samples-locus

Selected a subset of loci is less straight forward when working with a denovo assembly. Here the locus numbers are not very meaningful,
except that subcomponent numbering indicates likely relationships of paralogs.

Let's say you wanted to sample 1000 random loci that meet your filtering requirements. You could do...

```bash
ipyrad2 seqex \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-seqex \
  -N 1000 \
  -C
```

??? note "Stats file"
    ...

### imap/minmap
```bash
ipyrad2 seqex \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-seqex \
  -n assembly_min8 \
  --imap IMAP.tsv \
  --minmap MINMAP.tsv \
  -C
```

??? note "Stats file"
    ...
