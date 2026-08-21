# BPP MSC-model

BPP is a powerful tool for analyzing multi-locus sequence data under a variety of evolutionary
models. It can be used to estimate parameters of a species tree or network; infer a species tree
topology; and test species delimitation hypotheses.

Although the short length of individual RAD-seq loci limits the amount of information that each
contains -- usually a few informative SNPs per locus -- they can serve as suitable input to a
multi-locus analysis in BPP. Many hundreds or thousands of loci together can capture the
distribution of coalescent variation in one or multiple populations.

`ipyrad2` makes it easy to run BPP. There are two distinct ways to do this. First, you can
simply generate a properly formatted sequence file for BPP, and then go run BPP externally.
Or, alternatively, you can use our built-in `bpp` tool to not only generate the formatted
sequence file, but also generate the control-file and optionally start the BPP run. Both
modes are demonstrated in this tutorial.

## Requirements

- A completed ipyrad2 assembly HDF5 file
- An IMAP file to optionally subsample and assign individuals to populations
- An understanding of the BPP run settings
- Optionally download BPP v.4.8.7 or newer. Or `ipyrad2` can auto-download it for you.


## (seqex) write seqs only

Here we use `seqex` to export assembled loci in a format appropriate for BPP, which requires
a specific sequence format. We supply arguments supported by seqex to keep only samples that
are listed in an IMAP file, and only loci that data for some number of samples in each population,
as specified in a MINMAP file. This is a way to filter the dataset to minimize missing data and
keep loci that are informative across your dataset. The `-O phy` and `--append-population` args
are key here to format the data for BPP. We also use `-N 200` and `-s 12345` to specify to
randomly sample at most 200 loci that pass filtering, using a repeatable random seed. This is
because a smaller dataset will run faster in BPP for this tutorial.

```bash
ipyrad2 seqex \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-bpp/ \
  -i SRP021469/IMAP2.txt \
  -g SRP021469/MINMAP2.txt \
  -n assembly_min8 \
  -m 8 \
  -N 200 \
  -s 12345 \
  -O phy \
  --append-population
```

This will print a log and generate a stats file, both shown below.

??? "ipyrad2 seqex run for bpp log"

    ```parsed-literal
    2026-07-23 16:33:05 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 16:33:05 | INFO     | cli_analysis.py      | ---- ipyrad2 seqex: extract filtered delimited loci ----
    2026-07-23 16:33:05 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 16:33:05 | INFO     | cli_analysis.py      | CMD: ipyrad2 seqex -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-bpp/ -i SRP021469/IMAP2.txt -g SRP021469/MINMAP2.txt -n assembly_min8 -m 8 -N 200 -s 12345 -O phy --append-population
    2026-07-23 16:33:05 | INFO     | sequence_common.py   | imap file doesn't include minmap info, parsing standard imap file format.
    2026-07-23 16:33:05 | INFO     | sequence_windows.py  | No windows specified; selecting the full length of all scaffolds. Use -w to subset scaffold windows and -P to view scaffold names.
    2026-07-23 16:33:05 | INFO     | sequence_windows.py  | selected 45448 windows from 45448 scaffolds
    2026-07-23 16:33:15 | INFO     | seqex.py             | wrote 200 filtered loci as independent records in one PHYLIP file to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/assembly_min8.phy
    2026-07-23 16:33:15 | INFO     | seqex.py             | wrote stats report to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/assembly_min8.stats.txt
    ```


??? "ipyrad2 seqex stats file"

    ```parsed-literal
    # Seqex Summary
    command: ipyrad2 seqex -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-bpp/ -i SRP021469/IMAP2.txt -g SRP021469/MINMAP2.txt -n assembly_min8 -m 8 -N 200 -s 12345 -O phy --append-population
    data: SRP021469/OUT/assembly.hdf5
    output_layout: multi-locus
    out_format: phy
    cores: 1
    max_loci: 200
    random_seed: 12345
    min_length: none
    clipping_mode: automatic
    coordinate_clipping_applied: false
    windows_selected: 45448
    selected_windows: none
    candidate_loci: 45448
    rejected_raw_length: 0
    rejected_locus_coverage: 34026
    rejected_site_coverage: 0
    rejected_sample_missing: 0
    rejected_filtered_length: 0
    accepted_before_sampling: 11422
    written_loci: 200

    # Output Summary
    total_sites_written: 13251
    total_bases_written: 145761
    full_matrix_bases: 145761
    non_missing_bases: 139379
    non_missing_occupancy: 0.956216
    max_samples: 11
    mean_samples: 11.000000

    # Sample Occupancy
    sample                  population           written_final  loci_written  loci_dropped_by_r  matrix_bases  non_missing_bases  non_missing_occupancy
    ----------------------  -------------------  -------------  ------------  -----------------  ------------  -----------------  ---------------------
    29154_superba           superba              yes                     200                  0         13251              13251               1.000000
    30686_cyathophylla      cyathophylla         yes                     200                  0         13251              13251               1.000000
    32082_przewalskii       przewalskii          yes                     200                  0         13251              10962               0.827258
    33588_przewalskii       przewalskii          yes                     200                  0         13251              11901               0.898121
    35236_rex               rex_subsp_rockii     yes                     200                  0         13251              13251               1.000000
    35855_rex               rex_subsp_rex        yes                     200                  0         13251              12885               0.972379
    38362_rex               rex_subsp_lipskyana  yes                     200                  0         13251              13180               0.994642
    39618_rex               rex_subsp_lipskyana  yes                     200                  0         13251              11719               0.884386
    40578_rex               rex_subsp_rex        yes                     200                  0         13251              13177               0.994416
    41478_cyathophylloides  cyathophylloides     yes                     200                  0         13251              13182               0.994793
    41954_cyathophylloides  cyathophylloides     yes                     200                  0         13251              12620               0.952381

    # Written Loci
    locus_index  locus                source_locus         selected_window      clipped  raw_samples  raw_sites  filtered_samples  filtered_sites  concat_start  concat_end
    -----------  -------------------  -------------------  -------------------  -------  -----------  ---------  ----------------  --------------  ------------  ----------
              1  locus_623_1:1-69     locus_623_1:1-69     locus_623_1:1-69     no                11         69                11              69
              2  locus_968_1:1-69     locus_968_1:1-69     locus_968_1:1-69     no                11         69                11              69
              3  locus_1028_1:1-69    locus_1028_1:1-69    locus_1028_1:1-69    no                11         69                11              69
              4  locus_1226_1:1-69    locus_1226_1:1-69    locus_1226_1:1-69    no                11         69                11              62
              5  locus_1467_1:1-69    locus_1467_1:1-69    locus_1467_1:1-69    no                11         69                11              69
              ...
    ```

Let's peek at the output file. You can see it is a multi-locus phylip file with `population^sample` naming
convention. This sequence can be used as input to BPP.

```bash
head -n 100 SRP021469/output-bpp/assembly_min8.phy
```

```parsed-literal
11 69
superba^29154_superba                        GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCAAAATCGATCACCA
cyathophylla^30686_cyathophylla              GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCAAAATCGATCACCA
przewalskii^32082_przewalskii                GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAGGTCGCCTCAAAATCGATCACCA
przewalskii^33588_przewalskii                GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAGGTCGCCTCAAAATCGATCACCA
rex_subsp_rockii^35236_rex                   GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCGAAATCGATCACCA
rex_subsp_rex^35855_rex                      GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCGAAATCGATCACCA
rex_subsp_lipskyana^38362_rex                GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCGAAATCRATCACCA
rex_subsp_lipskyana^39618_rex                GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCGAAATCGATCACCA
rex_subsp_rex^40578_rex                      GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCGAAATCGATCACCA
cyathophylloides^41478_cyathophylloides      GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCAAAATCGATCACCA
cyathophylloides^41954_cyathophylloides      GGAACTCAATAATCTCTTGGGGATGTGGGTTGTTCTCTTTATCACAAGTCGCCTCAAAATCGATCACCA

11 69
superba^29154_superba                        CAATAAATTACTTTCTAAATGAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
cyathophylla^30686_cyathophylla              CAATAAATTACTTTCTAAATGAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
przewalskii^32082_przewalskii                CAATAAATTACTTTCTCACTGAGATGGCAGGGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
przewalskii^33588_przewalskii                CAATAAATTACTTTCTCACTGAGATGGCAGGGTGAAGTACGAGTTTGAGATATTTGATTCAATTGACAC
rex_subsp_rockii^35236_rex                   CAATAAATTACTTTCTAAATAAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
rex_subsp_rex^35855_rex                      CAATAAATTACTTTCTAAATAAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
rex_subsp_lipskyana^38362_rex                CAATAAATTACTTTCTAAATAAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
rex_subsp_lipskyana^39618_rex                CAATAAATTACTTTCTAAATAAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
rex_subsp_rex^40578_rex                      CAATAAATTACTTTCTAAATAAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
cyathophylloides^41478_cyathophylloides      CAATAAATTACTTTCTAAATGAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC
cyathophylloides^41954_cyathophylloides      CAATAAATTACTTTCTAAATGAGATGGCAGAGTGAAGTATGAGTTTGAGATATTTGATTCAATTGACAC

11 69
superba^29154_superba                        CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTTCCGCGGGCG
cyathophylla^30686_cyathophylla              CTGATTTGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCACTCATTGCCGCGGGCG
przewalskii^32082_przewalskii                CTGATTCGGTTAAAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTGCCGCGGGAG
przewalskii^33588_przewalskii                CTGATTCGGTTAAAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTGCCGCGGGAG
rex_subsp_rockii^35236_rex                   CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTGCCGCGGGCG
rex_subsp_rex^35855_rex                      CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTGCCGCGGGCG
rex_subsp_lipskyana^38362_rex                CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCACTCATTGCCGCGGGCG
rex_subsp_lipskyana^39618_rex                CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCACTCATTGCCGCGGGCG
rex_subsp_rex^40578_rex                      CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTGCCGCGGGCG
cyathophylloides^41478_cyathophylloides      CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTTCCGCGGGCG
cyathophylloides^41954_cyathophylloides      CTGATTCGGTTAGAGGAAAAAAGTGCACTGCTTATCCGACTATTAAACCTGCGCTCATTTCCGCGGGCG

11 62
superba^29154_superba                        CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
cyathophylla^30686_cyathophylla              CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
przewalskii^32082_przewalskii                CACATCAGCCCATACTTCCTGATGATTGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACN
przewalskii^33588_przewalskii                CACATCAGCCCATACTTCCTGATGATTGCGTGTGGATTTCCACAGACCCGGCTGCAAAAACG
rex_subsp_rockii^35236_rex                   CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
rex_subsp_rex^35855_rex                      CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
rex_subsp_lipskyana^38362_rex                CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
rex_subsp_lipskyana^39618_rex                CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
rex_subsp_rex^40578_rex                      CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
cyathophylloides^41478_cyathophylloides      CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
cyathophylloides^41954_cyathophylloides      CACATCAGGCCATACTTCCTGATGATAGCGTGTGGGTTTCCACAGACCCGGCTGCAAAAACG
```

You can now use this file as an input to BPP. You will also have to generate a BPP control (.ctl)
file that specifies the path to this file, the number of loci that are in it, and many specifications
of how you want the program to run.

To expedite that process, we provide a wrapper tool called `ipyrad2 bpp`, described below.


## (bpp) write seqs/ctl/map & run

The `ipyrad2 bpp` accepts similar argument as `seqex` for filtering and formatting the sequence data,
but it also includes many additional argments for specifying the BPP run. It will create a CTL file
and fill it with the path to the written sequence file, and you can use the command argument flags to
further specify the parameters of the BPP run. By default, this tool will start the BPP run using the
resulting files, unless you specify `--write-only`. If BPP is not installed it will download the latest
version into your TMPDIR.


!!! warning "important"

    Note: BPP is a very complex program. You should read the BPP documentation, and you should
    always inspect the BPP CTL file to ensure that your analysis is configured properly.


Here let's run A01 to infer a species tree while keeping species assignments fixed. I increased the
values of some run parameters to run the BPP MCMC chain longer than the default.
This runs for about 20 minutes on my laptop.

```bash
ipyrad2 bpp \
    -d SRP021469/OUT/assembly.hdf5 \
    -o SRP021469/output-bpp/ \
    -n species-tree \
    --tree SRP021469/TREE.nwk \
    --imap SRP021469/IMAP2.txt \
    --minmap SRP021469/MINMAP2.txt \
    --max-loci 1000 \
    --min-length 50 \
    --analysis A01 \
    --species-model-prior 1 \
    --substitution-model HKY \
    --locus-seed 102 \
    --mcmc-seed 202 \
    --tau-prior invgamma 3 0.03 \
    --samplefreq 5 \
    --nsample 50_000 \
    --burnin 10000 \
    --threads 6
```

??? "ipyrad2 bpp log"

    ```literal
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | ---- ipyrad2 bpp: single-run BPP staging and execution ----
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | CMD: ipyrad2 bpp -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-bpp/ -n msc-infer --imap SRP021469/IMAP2.txt --minmap SRP021469/MINMAP2.txt --min-length 50 --threads 6 --analysis A01 --locus-seed 123 --mcmc-seed 123 --tau-prior invgamma 3 0.05 -N 100 --tree SRP021469/TREE.nwk --samplefreq 5 --nsample 50_000 --burnin 10000
    2026-07-23 15:56:33 | INFO     | sequence_common.py   | imap file doesn't include minmap info, parsing standard imap file format.
    2026-07-23 15:56:33 | INFO     | sequence_windows.py  | No windows specified; selecting the full length of all scaffolds. Use -w to subset scaffold windows and -P to view scaffold names.
    2026-07-23 15:56:33 | INFO     | sequence_windows.py  | selected 45448 windows from 45448 scaffolds
    [####################] 100% | Filtering seqex loci - total jobs: 44686
    2026-07-23 15:56:36 | INFO     | seqex.py             | wrote 100 filtered loci as independent records in one PHYLIP file to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/msc-infer.phy
    2026-07-23 15:56:36 | INFO     | seqex.py             | wrote stats report to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/msc-infer.stats.txt
    2026-07-23 15:56:37 | INFO     | bpp.py               | bpp v4.8.7 (/tmp/bpp-4.8.7-linux-x86_64/bin/bpp)
    2026-07-23 16:36:24 | INFO     | bpp.py               | wrote BPP outputs to /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp
    ```


This writes several output files. The first few are generated by seqex running internally,
which writes `{name}.stats.txt` and `{name}.stats.json`, and `{name}.phy`. Then ipyrad2 also writes
the map and configuration file, `{name}.imapfile.txt`, and `{name}.ctl.txt`.
Finally, if BPP starts running, it can generate many result files, which vary depending
on the settings and algorithm run.

In this example, BPP generates a result summary file labeled `{name}.txt`. It includes
information about the data, MCMC chain, convergence diagnostics, and summarized species-tree
results. Per-locus gene trees are disabled by default because they can create many large files;
add `--print-gene-trees` when those outputs are needed.

```bash
tail -n 100 SRP021469/output-bpp/species-tree.txt
```

```parsed-literal
Current finetune:  46.52083 0.01089 0.00217 0.00583 0.00456 0.10146 2.95025 0.19043 0.11375 0.18844 0.82254 0.68950
New finetune:      94.22963 0.01119 0.00219 0.00518 0.00486 0.12955 2.65636 0.16875 0.11217 0.21042 0.72622 0.75644

=> 'finetune = 1 Gage:94.229629 Gspr:0.011194 th1:0.002191 th2:0.005182 tau:0.004862 mix:0.129548 alfa:2.656356 mubr:0.168752 nubr:0.112169 mu_i:0.210419 nu_i:0.726223 brte:0.756444'

   5%  0.51 0.30 0.29 0.31 1.00 0.28 0.27 0.27 0.32 0.27 0.33 0.30 0.31 0.0104 0.0000   0.0370  0.0618   8342.76319  -12504.03254  2:35
  10%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0090 0.0000   0.0323  0.0647   8812.76285  -12505.10366  4:25
  15%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0084 0.0000   0.0308  0.0644   8521.82421  -12505.03592  6:09
  20%  0.51 0.30 0.29 0.29 1.00 0.28 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0092 0.0000   0.0311  0.0636   8446.80192  -12504.98805  8:04
  25%  0.51 0.29 0.28 0.29 1.00 0.28 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0094 0.0000   0.0309  0.0623   8659.26399  -12505.28271  9:56
  30%  0.51 0.29 0.28 0.28 1.00 0.27 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0096 0.0000   0.0300  0.0602   8571.22293  -12505.08117  11:58
  35%  0.51 0.29 0.27 0.28 1.00 0.27 0.27 0.28 0.32 0.28 0.33 0.30 0.31 0.0097 0.0000   0.0292  0.0597   8514.24785  -12505.25006  13:55
  40%  0.51 0.29 0.28 0.28 1.00 0.27 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0098 0.0000   0.0305  0.0606   7593.46964  -12505.28388  15:49
  45%  0.51 0.29 0.28 0.29 1.00 0.28 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0101 0.0000   0.0299  0.0636   8160.77893  -12505.30766  17:39
  50%  0.51 0.30 0.29 0.29 1.00 0.28 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0100 0.0000   0.0298  0.0647   8755.95614  -12505.43479  19:38
  55%  0.51 0.30 0.29 0.29 1.00 0.28 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0108 0.0000   0.0290  0.0648   8401.52223  -12505.50856  21:40
  60%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0105 0.0000   0.0287  0.0657   8093.53117  -12505.37230  23:35
  65%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0106 0.0000   0.0287  0.0657   8182.35936  -12505.33935  25:37
  70%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0104 0.0000   0.0278  0.0686   8961.59914  -12505.68131  27:36
  75%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0103 0.0000   0.0288  0.0678   8256.11679  -12505.64574  29:39
  80%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0105 0.0000   0.0286  0.0681   8129.77656  -12505.47247  31:33
  85%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0102 0.0000   0.0281  0.0682   8310.53737  -12505.61942  33:37
  90%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0100 0.0000   0.0282  0.0677   8491.17808  -12505.53850  35:38
  95%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0100 0.0000   0.0285  0.0674   8243.79313  -12505.49787  37:42
 100%  0.51 0.30 0.29 0.29 1.00 0.29 0.27 0.28 0.32 0.28 0.32 0.30 0.31 0.0098 0.0000   0.0291  0.0672   8340.54215  -12505.51091  39:47

39:47 spent in MCMC

Species in order:
   1. ___cyathophylla
   2. ___cyathophylloides
   3. ___przewalskii
   4. ___rex_subsp_lipskyana
   5. ___rex_subsp_rex
   6. ___rex_subsp_rockii
   7. ___superba

(A) Best trees in the sample (9 distinct trees in all)
    32616  0.65231  0.65231 ((((___cyathophylla, ___superba), ___cyathophylloides), ((___rex_subsp_lipskyana, ___rex_subsp_rockii), ___rex_subsp_rex)), ___przewalskii);
    12711  0.25421  0.90652 ((((___cyathophylla, ___superba), ___cyathophylloides), (___rex_subsp_lipskyana, (___rex_subsp_rex, ___rex_subsp_rockii))), ___przewalskii);
     4014  0.08028  0.98680 ((((___cyathophylla, ___superba), ___cyathophylloides), ((___rex_subsp_lipskyana, ___rex_subsp_rex), ___rex_subsp_rockii)), ___przewalskii);
      298  0.00596  0.99276 (((___cyathophylla, ___superba), (___cyathophylloides, ((___rex_subsp_lipskyana, ___rex_subsp_rockii), ___rex_subsp_rex))), ___przewalskii);
      121  0.00242  0.99518 ((((___cyathophylla, ___superba), ((___rex_subsp_lipskyana, ___rex_subsp_rockii), ___rex_subsp_rex)), ___cyathophylloides), ___przewalskii);
      107  0.00214  0.99732 (((___cyathophylla, ___superba), (___cyathophylloides, (___rex_subsp_lipskyana, (___rex_subsp_rex, ___rex_subsp_rockii)))), ___przewalskii);
      106  0.00212  0.99944 ((((___cyathophylla, ___superba), (___rex_subsp_lipskyana, (___rex_subsp_rex, ___rex_subsp_rockii))), ___cyathophylloides), ___przewalskii);
       27  0.00054  0.99998 (((___cyathophylla, ___superba), (___cyathophylloides, ((___rex_subsp_lipskyana, ___rex_subsp_rex), ___rex_subsp_rockii))), ___przewalskii);
        1  0.00002  1.00000 ((((___cyathophylla, ___cyathophylloides), ___superba), ((___rex_subsp_lipskyana, ___rex_subsp_rockii), ___rex_subsp_rex)), ___przewalskii);

(B) Best splits in the sample of trees (10 splits in all)
 50001 1.000000  0001110
 50001 1.000000  1101111
 50000 0.999980  1000001
 49342 0.986820  1100001
 33036 0.660707  0001010
 12924 0.258475  0000110
  4041 0.080818  0001100
   432 0.008640  0101110
   227 0.004540  1001111
     1 0.000020  1100000

(C) Majority-rule consensus tree
((((___cyathophylla, ___superba) #0.999980, ___cyathophylloides) #0.986820, ((___rex_subsp_lipskyana, ___rex_subsp_rockii) #0.660707, ___rex_subsp_rex) #1.000000) #1.000000, ___przewalskii);

(D) Best tree (or trees from the mastertree file) with support values
((((___cyathophylla, ___superba) #0.999980, ___cyathophylloides) #0.986820, ((___rex_subsp_lipskyana, ___rex_subsp_rockii) #0.660707, ___rex_subsp_rex) #1.000000) #1.000000, ___przewalskii);   [P = 0.652307]
```

Pulling out the Majority-rule consensus tree newick from this file, we can visual the resulting inferred population/species tree using toytree:

```bash
toytree view -i "((((___cyathophylla, ___superba) #0.999980, ___cyathophylloides) #0.986820, ((___rex_subsp_lipskyana, ___rex_subsp_rockii) #0.660707, ___rex_subsp_rex) #1.000000) #1.000000, ___przewalskii);"
```
```parsed-literal
                 ┌──────___cyathophylla
            ┌────┤
      ┌─────┤    └──────___superba
      │     │
      │     └─────___cyathophylloides
┌─────┤
│     │          ┌──────___rex_subsp_lipskyana
│     │     ┌────┤
│     └─────┤    └──────___rex_subsp_rockii
│           │
│           └─────___rex_subsp_rex
│
└──────___przewalskii
```

## Parameters

`ipyrad2 bpp` first uses the seqex engine to select and format loci, then writes
BPP input files and runs BPP unless `--write-only` is set. For the complete BPP
syntax, see the [BPP manual](https://bpp.github.io/bpp-manual/bpp-4-manual/).

!!! warning "Inspect the generated control file"

    Priors and chain lengths are starting values. Inspect `{name}.ctl.txt` and
    compare independent chains using different `--mcmc-seed` values.

### Input data and locus selection

| Argument | Required/default | Effect |
| --- | --- | --- |
| `-d, --data PATH` | required | Reads loci from an ipyrad2 assembly HDF5 file. |
| `-t, --tree NEWICK or PATH` | required | Guide tree whose tips exactly match IMAP population names. |
| `-i, --imap PATH` | required | Assigns samples or glob patterns to populations. |
| `-g, --minmap PATH` | required | Sets minimum sample coverage for each population. |
| `-N, --max-loci INT` | required | Caps the number of loci written after filtering. |
| `-L, --min-length INT` | required | Rejects loci shorter than this length. |
| `-o, --out PATH` | `output-bpp/` | Sets the result directory. |
| `-n, --name STR` | `bpp` | Sets the common filename and BPP job prefix. |

### Analysis model

`--analysis` is case-insensitive, so `A01`, `a01`, `MSC-I`, and `msc-i` are
accepted. The canonical modes are:

| Mode | BPP controls | Purpose |
| --- | --- | --- |
| `A00` | fixed delimitation and fixed tree | Estimate population sizes and divergence times on a fixed tree. |
| `A01` | fixed delimitation, inferred tree | Infer a species tree with fixed species assignments. |
| `A10` | inferred delimitation, fixed tree | Infer species boundaries on a guide tree. |
| `A11` | inferred delimitation and tree | Jointly infer species boundaries and the species tree. |
| `MSC-I` | fixed A00 + introgression | Estimate episodic introgression from extended Newick `&phi=` annotations. |
| `MSC-M` | fixed A00 + migration | Estimate directional migration on a fixed tree. |

| Argument | Default | Effect |
| --- | --- | --- |
| `--analysis MODE` | `A00` | Selects one of the six modes above. |
| `--migration SOURCE,TARGET` | none | Adds one MSC-M migration band; repeat for multiple bands. |
| `--species-model-prior INT` | `1` when relevant | Uses prior codes 0–1 for A01/A10 and 0–3 for A11. |
| `--delimitation-proposal TOKENS` | `0 2` | Uses A10/A11 proposal `0 EPSILON` or `1 ALPHA M`. |
| `--constraint-file PATH` | none | Adds A01/A11 `define`, `constraint`, and `outgroup` statements. Population names are validated and staged in BPP's namespace. |

MSC-I and MSC-M are fixed-tree A00 analyses. `--migration` is required only
for MSC-M. Constraints are supported only for A01 and A11.

### Evolutionary models and priors

| Argument | Default | Effect |
| --- | --- | --- |
| `--theta-prior TOKENS` | `invgamma 3 0.002` | Supports `invgamma A B [estimate or integrate]`, `gamma A B`, and `beta P Q LOWER UPPER`. |
| `--tau-prior DIST A B` | `invgamma 3 0.03` | Sets the root divergence-time prior using inverse-Gamma or Gamma parameters. |
| `--phi-prior A B` | `1 1` | Sets the MSC-I Beta prior on introgression probabilities. |
| `--w-prior A B` | `2 200` | Sets the MSC-M Gamma shape-rate prior on migration rates. |
| `--substitution-model MODEL` | `JC69` | Selects JC69, K80, F81, HKY, T92, TN93, F84, or GTR. |
| `--alpha-prior A B NCAT` | none | Enables discrete-Gamma among-site rate variation; omitting it keeps site rates equal. |
| `--theta-model MODEL` | `linked-none` | Selects `linked-none`, `linked-all`, `linked-inner`, `linked-msci`, or `linked-mscm`. |
| `--clean-data` | off | Removes sites containing ambiguity characters. |

`linked-msci` and `linked-mscm` require the corresponding MSC-I or MSC-M mode.
Theta and tau prior parameters must be positive. Appending `integrate` to an
inverse-Gamma theta prior writes BPP's `int` option.

### Rate and clock models

| Argument | Default | Supported syntax |
| --- | --- | --- |
| `--locus-rate TOKENS` | `0` | `0` fixes equal rates; `1 A B C iid or dir` estimates rates among loci. |
| `--clock TOKENS` | `1` | `1` is strict; modes `2`/`3` accept `A B C iid or dir G or LN`; mode `4 A` is also supported. |

When both settings contain `iid` or `dir`, those values must match. Fixed
external locus-rate files and tip dating are not exposed because their required
inputs cannot be mapped safely through locus selection. Clock mode 3 is rejected
for MSC-I.

### MCMC and execution

| Argument | Default | Effect |
| --- | --- | --- |
| `--burnin INT` | `10000` | Automatic-finetuning burn-in; values must be greater than 200. |
| `--samplefreq INT` | `10` | Iterations between retained samples. |
| `--nsample INT` | `100000` | Number of retained samples. |
| `--threads N [START STEP]` | automatic | Sets BPP's one- or three-value thread specification; the count is also used by seqex. |
| `--locus-seed NONE or INT` | none | Controls only random locus subsampling. |
| `--mcmc-seed NONE or INT` | none | Controls only the BPP MCMC stream; use different values for independent chains. |
| `--checkpoint START [INTERVAL]` | none | Writes BPP checkpoint files. Resume remains a direct BPP operation. |
| `--bpp-binary PATH` | auto | Uses an explicit executable; execution requires BPP 4.8.7 or newer. |

### Additional output options and logging

| Argument | Default | Effect |
| --- | --- | --- |
| `--print-locus-rates` | off | Includes sampled locus rates in MCMC output. |
| `--print-gene-trees` | off | Writes per-locus gene trees; output can be large. |
| `--print-substitution-parameters` | off | Includes sampled substitution parameters. |
| `--write-only` | off | Writes inputs without resolving, downloading, or running BPP. |
| `-f, --force` | off | Allows same-name files to be overwritten. |
| `-l, --log-level LEVEL` | `INFO` | Sets ipyrad2 logging verbosity. |

The default BPP print setting is `1 0 0 0 0`, so optional gene-tree and rate
outputs are disabled unless requested. Heredity scalars, trait analyses, custom
per-locus models, prior-only runs, checkpoint resume, and Bayes-factor workflows
are not exposed by this wrapper.

## Examples

The following commands use generic input filenames to show parameter combinations
appropriate for each supported model. Replace the paths with your own assembly,
tree, IMAP, and MINMAP files. `-t` is shorthand for `--tree`.

### A00: fixed tree and assignments

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t species.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 \
    -L 100 \
    --analysis A00 \
    --theta-prior invgamma 3 0.002 \
    --tau-prior invgamma 3 0.03 \
    --locus-seed 101 \
    --mcmc-seed 201
```

### A01: species-tree inference

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t species.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 \
    -L 100 \
    --analysis A01 \
    --species-model-prior 1 \
    --substitution-model HKY \
    --locus-seed 102 \
    --mcmc-seed 202
```

### A10: species delimitation on a guide tree

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t species.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 -L 100 \
    --analysis A10 \
    --species-model-prior 1 \
    --delimitation-proposal 1 2 1 \
    --locus-seed 103 \
    --mcmc-seed 203
```

### A11: joint delimitation and tree inference

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t species.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 -L 100 \
    --analysis A11 \
    --species-model-prior 3 \
    --delimitation-proposal 1 2 1 \
    --constraint-file constraints.txt \
    --locus-seed 104 \
    --mcmc-seed 204
```

### MSC-I: episodic introgression

The guide tree must be extended Newick containing BPP `&phi=` annotations.

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t msci.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 \
    -L 100 \
    --analysis MSC-I \
    --phi-prior 1 1 \
    --theta-model linked-msci \
    --locus-seed 105 \
    --mcmc-seed 205
```

### MSC-M: directional migration

Migration pairs may reference tips or appropriately labeled ancestral clades.

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t species.nwk \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 \
    -L 100 \
    --analysis MSC-M \
    --migration A,B --migration B,A \
    --w-prior 2 200 \
    --theta-model linked-mscm \
    --locus-seed 106 \
    --mcmc-seed 206
```

### MSC-M: migration from the A/B ancestor into C

To model one unidirectional migration band from the common ancestor of A and B
into C, label that ancestral branch in the guide tree:

```text
((A,B)AB,C)R;
```

The internal label `AB` is required because the migration argument refers to
that ancestor. With `IMAP.txt` and `MINMAP.txt` defining populations A, B, and C,
run:

```bash
ipyrad2 bpp \
    -d assembly.hdf5 \
    -t '((A,B)AB,C)R;' \
    -i IMAP.txt \
    -g MINMAP.txt \
    -N 1000 \
    -L 100 \
    --analysis MSC-M \
    --migration AB,C \
    --w-prior 2 200 \
    --theta-model linked-mscm \
    --locus-seed 107 \
    --mcmc-seed 207
```

`--migration SOURCE,TARGET` uses BPP's donor/source-to-recipient/target
convention. Thus `--migration AB,C` specifies migration from ancestral
population `AB` into population `C`, and the corresponding BPP parameter is
(w_{AB,C}). In the forward-time biological interpretation, C receives
migrants from AB. In the backward-time coalescent representation, a lineage
sampled in C can trace back into AB; this is the same model direction viewed
backward in time.

Interpret the posterior samples and summary for (w_{AB,C}) rather than reading
the argument as an undirected connection. A posterior concentrated away from
zero supports ongoing migration in the specified AB → C direction. This command
does not estimate C → AB migration: use `--migration C,AB` for the reverse model,
or specify both `--migration AB,C --migration C,AB` to fit two directional bands.
Comparing those explicitly fitted alternatives is necessary when assessing
whether the data support directionality.

The Recipes section will include additional examples demonstrating complete
runs on real data, including dataset-specific filtering, priors, chain lengths,
and interpretation of BPP output.
