# BPP MSC-model

A multi-locus RAD-seq dataset can be used to estimate parameters of a multi-species coalescent (MSC) model.
A widely used tool for this is BPP [link]. We provide a simple analysis tool that makes it easy to
implement BPP.

## Setup

The ipyrad2 tool only stages a BPP run by generating the necessary files. To run it you can then
call BPP on the generated files as demonstrated below. For this, you will need to install
BPP v.4.4 or greater.

## Requirements

- A completed ipyrad2 assembly HDF5 file
- An IMAP file to optionally subsample and assign individuals to populations
- An understanding of the BPP run settings

## Generated files

- {name}.ctl
- {name}.phy
- {name}.imap
- {name}.stats

## Command

```bash
ipyrad2 bpp \
  -d assembly.hdf5 \
  -o output-bpp/ \
  -n msc-fit \
  --tree species.nwk \
  --imap IMAP.txt \
  --minmap MINMAP.txt \
  --maxloci 1000 \
  --min-length 100 \
  --threads 50 \
  --seed 123 \
  --tauprior 3 0.05
```

The logging report shows:

```

```

## run

```bash
bpp --cfile output-bpp/msc-fit.ctl
```