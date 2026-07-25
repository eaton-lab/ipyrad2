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
- Optionally download BPP v.4.8.6 or newer. Or `ipyrad2` can auto-download it for you.


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


Here let's run the `--speciestree` mode (model "10" in BPP) to infer a species tree. I increased the
values of some run parameters to run the BPP MCMC chain longer than the default.
This runs for about 20 minutes on my laptop.

```bash
ipyrad2 bpp \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-bpp/ \
  -n msc-infer \
  --tree SRP021469/TREE.nwk \
  --imap SRP021469/IMAP2.txt \
  --minmap SRP021469/MINMAP2.txt \
  --min-length 50 \
  --threads 6 \
  --seed 123 \
  --tauprior 3 0.05 \
  --max-loci 100 \
  --speciestree \
  --samplefreq 5 \
  --nsample 50_000 \
  --burnin 10000
```

??? "ipyrad2 bpp log"

    ```literal
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | ---- ipyrad2 bpp: single-run BPP staging and execution ----
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | -------------------------------------------------------
    2026-07-23 15:56:33 | INFO     | cli_analysis.py      | CMD: ipyrad2 bpp -d SRP021469/OUT/assembly.hdf5 -o SRP021469/output-bpp/ -n msc-infer --imap SRP021469/IMAP2.txt --minmap SRP021469/MINMAP2.txt --min-length 50 --threads 6 --seed 123 --tauprior 3 0.05 -N 100 --tree SRP021469/TREE.nwk --speciestree --samplefreq 5 --nsample 50_000 --burnin 10000
    2026-07-23 15:56:33 | INFO     | sequence_common.py   | imap file doesn't include minmap info, parsing standard imap file format.
    2026-07-23 15:56:33 | INFO     | sequence_windows.py  | No windows specified; selecting the full length of all scaffolds. Use -w to subset scaffold windows and -P to view scaffold names.
    2026-07-23 15:56:33 | INFO     | sequence_windows.py  | selected 45448 windows from 45448 scaffolds
    [####################] 100% | Filtering seqex loci - total jobs: 44686
    2026-07-23 15:56:36 | INFO     | seqex.py             | wrote 100 filtered loci as independent records in one PHYLIP file to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/msc-infer.phy
    2026-07-23 15:56:36 | INFO     | seqex.py             | wrote stats report to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp/msc-infer.stats.txt
    2026-07-23 15:56:37 | INFO     | bpp.py               | bpp v4.8.6 (/tmp/bpp-4.8.6-linux-x86_64/bin/bpp)
    2026-07-23 16:36:24 | INFO     | bpp.py               | wrote BPP outputs to /home/deren/Documents/tools/ipyrad2/SRP021469/output-bpp
    ```


This writes several output files. The first few are generated by seqex running internally,
which writes `{name}.stats.txt` and `{name}.stats.json`. Then ipyrad2 writes the BPP data
and configuration files `{name}.phy` `{name}.imap.txt` and `{name}.ctl.txt`.
Finally, if BPP is run, it generates many result files, which will vary depending on the settings
and algorithm run.

In this example, BPP generated a posterior distribution of gene trees for each locus,
labeled as `{name}.gtree.L{locus_number}`. And it generated a result summary file labeled
`{name}.txt`. This file includes a lot of information, including a summary of the data itself,
the MCMC chain, the convergence diagnostics, and at the very end, the summarized results, which
in this case is a set newick trees.

```bash
tail -n 100 SRP021469/output-bpp/msc-infer.txt
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

The `ipyrad2 bpp` arguments control two stages of the analysis. First, `seqex`
selects and formats loci from an ipyrad2 assembly. Then, ipyrad2 writes a BPP
control file and, unless `--write-only` is used, runs BPP. The tables below
describe both the immediate effect of each argument and the BPP setting it
ultimately controls.

!!! warning "Choose priors and chain lengths for your data"

    The defaults are starting values, not generally appropriate values for every
    dataset. Read the [BPP documentation](https://bpp.github.io/bpp-manual/bpp-4-manual/),
    inspect the generated `{name}.ctl.txt` file, and run independent chains with
    different seeds to assess convergence.

### Inputs and locus selection

| Argument | Required/default | Effect |
| --- | --- | --- |
| `-d, --data PATH` | required | Reads loci from an ipyrad2 assembly HDF5 file. |
| `-o, --out PATH` | `output-bpp/` | Sets the directory for the staged sequence, map, control, statistics, and BPP result files. |
| `-n, --name STR` | `bpp` | Sets the common filename prefix, for example `bpp.phy`, `bpp.ctl.txt`, and `bpp.mcmc.txt`. |
| `--tree STR or PATH` | required | Supplies a guide tree as inline Newick or as a path to a Newick file. Its tip labels must exactly match the population names in the IMAP. MSC-I requires extended Newick annotations, and MSC-M migration involving ancestral populations requires suitable internal-node or clade references. |
| `-i, --imap PATH` | required | Assigns assembly samples to populations using `sample<TAB>population` or `glob<TAB>population`. Only matched samples are included, and the population names define the species used by BPP. |
| `-g, --minmap PATH` | required | Sets the minimum number of samples with data required from each population at a locus using `population<TAB>min`. Its population keys must match the IMAP keys. |
| `-N, --max-loci INT` | required | Caps the number of loci written after filtering. If more loci pass, a random subset is selected; if fewer pass, all passing loci are written and the control file records the smaller, actual `nloci`. |
| `-L, --min-length INT` | required | Rejects loci shorter than this number of base pairs, both before and after sequence filtering. |

The same value supplied to `--seed` controls random locus subsampling, making
the selected set repeatable when a positive seed is used.

### Model selection

The model flags are mutually exclusive. With no model flag, BPP uses the
standard multispecies coalescent with the supplied species tree and species
assignments fixed, and estimates population-size and divergence-time
parameters.

| Argument | Default | Effect |
| --- | --- | --- |
| `--msc-i` | off | Fits an MSC-I model with episodic introgression on a fixed network. The guide tree must use BPP extended Newick containing `&phi=` introgression annotations. This enables gene flow and adds `phiprior` to the control file. |
| `--msc-m SOURCE,TARGET [...]` | off | Fits an MSC-M isolation-with-migration model with continuous, directional migration on a fixed tree. Each token adds one migration band from the source population or clade to the target, for example `A,B` or `(A,B),C`. This enables gene flow and adds `wprior` and the migration pairs to the control file. |
| `--speciestree` | off | Estimates the species-tree topology while keeping the sample-to-species assignments fixed. The tree supplied with `--tree` provides the starting topology. |
| `--speciesdelimitation` | off | Estimates species boundaries on the fixed guide-tree topology by allowing guide-tree nodes to collapse. ipyrad2 writes BPP's delimitation proposal setting `1 0 2`. |

### Priors

Numeric distribution parameters require positive values; `--speciesmodelprior` instead accepts codes `0` through `3`. BPP uses shape-rate
parameterizations for the Gamma distribution and the inverse-Gamma
parameterization described below.

| Argument | Default | Effect |
| --- | --- | --- |
| `--thetaprior A B` | `3 0.03` | Assigns every contemporary and ancestral population-size parameter, theta, the inverse-Gamma prior `IG(A, B)`. Its mean is `B / (A - 1)`; increasing `B` raises the prior scale, while increasing `A` lowers it. BPP requires `A > 2` for a finite variance. ipyrad2 estimates theta rather than integrating it out. |
| `--tauprior A B` | `3 0.03` | Assigns the root divergence time, tau, the inverse-Gamma prior `IG(A, B)`, also with mean `B / (A - 1)`. Other divergence times are distributed conditionally within the root age, so this prior sets the overall time scale in expected substitutions per site. |
| `--speciesmodelprior INT` | `0` | Selects BPP's prior on tree/delimitation models: `0` gives equal probabilities to labeled histories; `1` gives equal probabilities to rooted species trees; `2` gives equal probabilities to numbers of species and apportions probability within each number in proportion to labeled histories; `3` instead apportions it uniformly among compatible models. Codes `2` and `3` apply only to joint species-tree and delimitation inference, which the current ipyrad2 model flags do not expose; use `0` or `1` with the available inference modes. |
| `--phiprior A B` | `1 1` | Under `--msc-i`, assigns each introgression probability phi a Beta prior with mean `A / (A + B)`. Larger `A + B` concentrates the prior more tightly around that mean. It has no effect unless `--msc-i` is selected. |
| `--wprior A B` | `2 200` | Under `--msc-m`, assigns each mutation-scaled migration rate `w = m / mu` a Gamma shape-rate prior with mean `A / B`. The default prior mean is `0.01`. It has no effect unless `--msc-m` is selected. |
| `--alphaprior A B NCAT` | `1 1 4` | Places a Gamma shape-rate prior on the shape parameter alpha for among-site rate variation; its prior mean is `A / B`. Conditional site-rate variance is `1 / alpha`, so smaller alpha permits more variation. `NCAT` is the number of discrete-Gamma rate categories: more categories improve the approximation but increase computation. |

### Rate variation

`--locusrate` and `--clock` pass their values directly to the corresponding BPP
control-file fields. ipyrad2 checks that tokens are present but leaves their
model-specific syntax and compatibility to BPP.

| Argument | Default | Effect |
| --- | --- | --- |
| `--locusrate TOKENS [...]` | `1 2 3 2 iid` | Controls substitution-rate variation among loci. Mode `0` fixes equal rates; mode `1 A B C iid or dir` estimates locus rates, with `A B` defining the Gamma prior on their mean and `C` controlling dispersion among loci; mode `2 FILE` reads fixed relative rates; mode `3 A B` is for tip dating. The ipyrad2 default estimates rates using a conditional-i.i.d. prior. Tip-dating mode is not usable through this wrapper because it does not expose BPP's `datefile`. |
| `--clock TOKENS [...]` | `2 10.0 100.0 5.0 iid LN` | Controls rate variation among species-tree branches. Mode `1` is a strict clock, mode `2` uses independent branch rates, and mode `3` uses autocorrelated rates. The runtime default is an independent-rates relaxed clock with a Gamma prior of mean `10 / 100 = 0.1` on the average rate variance, conditional-i.i.d. variation among loci with concentration `5`, and log-normal branch rates. The current CLI help displays `dir`, but the generated control file uses `iid`. |

When both settings specify `iid` or `dir`, BPP requires them to agree. Some
clock models are also incompatible with particular gene-flow models; consult
the BPP documentation before supplying custom token sequences.

### MCMC and runtime

The total number of MCMC iterations is:

```text
burnin + (nsample * samplefreq)
```

| Argument | Default | Effect |
| --- | --- | --- |
| `--burnin INT` | `1000` | Runs this many initial MCMC iterations without retaining samples. ipyrad2 enables automatic proposal tuning during burn-in, for which BPP requires more than 200 iterations. A longer burn-in provides more time to tune and approach the posterior but increases runtime. |
| `--samplefreq INT` | `2` | Retains one posterior sample after every `INT` post-burn-in iterations. Increasing it spaces retained samples farther apart and, for a fixed `nsample`, increases the total chain length and runtime. |
| `--nsample INT` | `10000` | Sets the number of posterior samples retained after burn-in. Increasing it can improve posterior precision and effective sample sizes, while increasing runtime and output size. |
| `--threads N [START INCREMENT]` | BPP default | With one positive integer, requests `N` BPP threads. Three integers additionally pin them beginning at hardware-thread `START` and stepping by `INCREMENT`. `N` also controls seqex preprocessing, cannot exceed `--max-loci`, and is reduced in the BPP control file if fewer than `N` loci pass filtering. If omitted, no `threads` line is written and BPP chooses its default distribution. |
| `--seed none or INT` | `none` | A positive integer makes both locus sampling and the BPP random-number stream repeatable. Omitting the option or using `none` lets locus sampling initialize nondeterministically and writes `seed = -1`, asking BPP to obtain an automatic seed. Use different automatic or positive seeds for independent convergence checks. |
| `--write-only` | off | Writes the BPP PHYLIP sequence file, IMAP file, control file, and seqex statistics without locating, downloading, or running the BPP executable. |
| `-f, --force` | off | Permits same-name staged files and BPP outputs in the output directory to be overwritten. Without it, existing outputs stop the run. |
| `-l, --log-level STR` | `INFO` | Sets ipyrad2 logging verbosity to `TRACE`, `DEBUG`, `INFO`, `WARNING`, or `ERROR`; it does not change which BPP results are written. |
| `-h, --help` | — | Prints the current command syntax, defaults, and examples, then exits. |


## Examples

See the Recipes section for more examples implementing different BPP algorithms [TODO].