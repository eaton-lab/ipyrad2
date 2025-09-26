
# ipyrad
Assemble RAD-seq datasets for evolutionary analysis.


```bash
ipyrad -h
```

The core command-line arguments:
```bash
# demultiplex reads by barcode/index to sample files
ipyrad demux -d ... -b ... -o ... -c 10 -t 2

# trim reads demultiplex reads by barcode/index to sample files
ipyrad trim -d ... -o ... -c 10 -t 2

# map reads to a reference genome to get filtered sorted bams
ipyrad map -d ... -r ... -o ... -c 10 -t 2

# assemble loci
ipyrad assemble -b ... -r ... -o ... -m 4 -d 5 -c 10 -t 2
```

