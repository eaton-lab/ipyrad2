
# ipyrad
Integrative assembly and analysis of RAD-seq datasets

## Usage
```bash
# help message for subcommand options
ipyrad2 -h
```

The core command-line arguments:
```bash
# demultiplex reads by barcode/index to sample files
ipyrad2 demux -d DATA/*.fastq.gz -b barcodes.tsv -o FASTQs/ -c 10 -t 2

# quality, adapter, and restriction overhang trimming
ipyrad2 trim -d FASTQs/*fastq.gz -o TRIMMED/ -c 10 -t 2

# map reads to reference genome to get filtered/sorted/marked bams
ipyrad2 map -d TRIMMED/*.fastq.gz -r REF.fa -o BAMs/ -c 10 -t 2

# assemble loci and call variants
ipyrad2 assemble -d BAMs/*.bam -r REF.fa -o OUT/ -m 4 -d 5 -c 10 -t 2
```

## Installation

```bash
# clone the development repo
git clone ...

# install dependencies
cd ipyrad2/
conda env create -f environment.yml -n ipyrad2
conda activate

# install local dev copy of ipyrad2
pip install -e . --no-deps
```
