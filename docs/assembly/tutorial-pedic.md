
# Empirical SE denovo assembly tutorial

Here we demonstrate a denovo assembly for an empirical RAD data set using ipyrad2.

The dataset involves 13 samples from the Eaton and Ree (2013) dataset, which is composed of
single-end 75bp reads from a RAD-seq library prepared with the PstI enzyme by Floragenex Inc.
The dataset includes all species within a small monophyletic clade of *Pedicularis*, including
multiple individuals from 5 species and several subspecies, as well as an outgroup species.
The sampling spans from population-level variation where species boundaries are unclear, to
higher-level divergence where species boundaries are quite distinct.
This is a common scale at which RAD-seq data are often very useful.

## Setup

If you haven’t done so yet, start by installing ipyrad2, as well as a few additional packages that will
be used to download the dataset from SRA, and run some downstream analyses.

```bash
conda install ipyrad2 sra-tools raxml-ng -c bioconda -c conda-forge
```

## Download the dataset

First we will download the metadata which includes sample accession IDs, names, and other information.
Here we use some bash commands to fetch the metadata from NCBI using a public URL for the study
accession (SRP021469) and save its metadata to a file (runinfo.csv), and also print a subset of it
to the terminal for viewing. The data could alternatively be fetched manually online.

```bash
# Create directories for the metadata and sequencing data
mkdir -p SRP021469/{sra,fastq,tmp}

# Download the study metadata
wget --no-verbose \
  -O SRP021469/runinfo.csv \
  'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP021469'

# Display the run accession, number of spots, and library name
cut -d',' -f1,4,29,30 SRP021469/runinfo.csv | tr ',' '\t'
```

This displays the following table where each row represents one sequenced sample, its number of reads (spots),
the scientific name, and the sample name that was assigned by the researchers.
```literal
Run          spots     ScientificName                 SampleName
SRR1754715   696994    Pedicularis superba            29154_superba
SRR1754720   1452316   Pedicularis thamnophila        30556_thamno
SRR1754730   1253109   Pedicularis cyathophylla       30686_cyathophylla
SRR1754729   964244    Pedicularis przewalskii        32082_przewalskii
SRR1754728   636625    Pedicularis thamnophila        33413_thamno
SRR1754727   1002923   Pedicularis przewalskii        33588_przewalskii
SRR1754731   1803858   Pedicularis rex                35236_rex
SRR1754726   1409843   Pedicularis rex                35855_rex
SRR1754725   1391175   Pedicularis rex                38362_rex
SRR1754723   822263    Pedicularis rex                39618_rex
SRR1754724   1707942   Pedicularis rex                40578_rex
SRR1754722   2199740   Pedicularis cyathophylloides   41478_cyathophylloides
SRR1754721   2199613   Pedicularis cyathophylloides   41954_cyathophylloides
```

Save the run accession and sample name, excluding the header, which we will use
to select and rename the downloaded sequences.

```bash
cut -d ',' -f1,30 SRP021469/runinfo.csv |
  tail -n +2 |
  tr -d '"' |
  tr ',' '\t' > SRP021469/samples.tsv
```

Then call the following bash script that uses `fasterq-dump` from the sra-tools package that we installed
to download the FASTQ files for each sample. This will probably take a few minutes at most, the total data
size is approximately 5Gb raw, and closer to 1Gb after cleanup and compression below.

```bash
while IFS=$'\t' read -r run sample; do
    echo "Downloading ${sample} (${run})"

    # Download the SRA data
    prefetch "$run" \
      --max-size u \
      --output-directory SRP021469/sra

    # Convert the SRA data to FASTQ
    fasterq-dump "SRP021469/sra/${run}" \
      --split-3 \
      --threads 4 \
      --temp SRP021469/tmp \
      --outdir SRP021469/fastq

    # Replace the run accession with the sample name
    for file in SRP021469/fastq/${run}*.fastq; do
        suffix=${file#SRP021469/fastq/${run}}
        mv "$file" "SRP021469/fastq/${sample}${suffix}"
    done

done < SRP021469/samples.tsv
```

Finally, let's compress the FASTQ files and view them. Note the path to our downloaded FASTQ data files is `SRP021469/fastq/`.
If the data were paired-end they would end in \_1.fastq.gz and \_2.fastq.gz, but here because it is single-end data it
produces just one .fastq.gz file per sample.

```bash
# compress the FASTQ files
gzip SRP021469/fastq/*.fastq

# clean up temp files to save space
rm -r SRP021469/sra/

# show the final files
ls -lh SRP021469/fastq/
```

```literal
total 1.1G
-rw-rw-r-- 1 deren deren  42M Jul 20 12:51 29154_superba.fastq.gz
-rw-rw-r-- 1 deren deren  85M Jul 20 12:51 30556_thamno.fastq.gz
-rw-rw-r-- 1 deren deren  79M Jul 20 12:51 30686_cyathophylla.fastq.gz
-rw-rw-r-- 1 deren deren  58M Jul 20 12:51 32082_przewalskii.fastq.gz
-rw-rw-r-- 1 deren deren  40M Jul 20 12:52 33413_thamno.fastq.gz
-rw-rw-r-- 1 deren deren  61M Jul 20 12:52 33588_przewalskii.fastq.gz
-rw-rw-r-- 1 deren deren 108M Jul 20 12:52 35236_rex.fastq.gz
-rw-rw-r-- 1 deren deren  84M Jul 20 12:52 35855_rex.fastq.gz
-rw-rw-r-- 1 deren deren  82M Jul 20 12:53 38362_rex.fastq.gz
-rw-rw-r-- 1 deren deren  50M Jul 20 12:53 39618_rex.fastq.gz
-rw-rw-r-- 1 deren deren 100M Jul 20 12:53 40578_rex.fastq.gz
-rw-rw-r-- 1 deren deren 127M Jul 20 12:53 41478_cyathophylloides.fastq.gz
-rw-rw-r-- 1 deren deren 134M Jul 20 12:53 41954_cyathophylloides.fastq.gz
```


## Assembly

### trim

#### run
The data that we downloaded is already demultiplexed to individual samples, so we can start by running read trimming.
Here we just enter the path to our data (-d) and the path where we want the trimmed reads to be written (-o). I also
specify the total number of cores to use (-c 12) and how to distribute these resources among threaded jobs (-t 4),
which specified to run 3 4-threaded jobs at a time. All other settings are left at the default.

```bash
ipyrad2 trim \
  -d SRP021469/fastq/*.fastq.gz \
  -o SRP021469/TRIM/ \
  -c 12 -t 4
```

#### logging

This will write a log to the terminal (stdout) describing the steps it is performing:

```literal
2026-07-20 13:07:44 | INFO     | cli_main.py          | ----------------------------------------------------------
2026-07-20 13:07:44 | INFO     | cli_main.py          | ----- ipyrad2 trim: quality, adapter, and cutsite motif trimming -----
2026-07-20 13:07:44 | INFO     | cli_main.py          | ----------------------------------------------------------
2026-07-20 13:07:44 | INFO     | cli_main.py          | CMD: ipyrad2 trim -d SRP021469/fastq/29154_superba.fastq.gz SRP021469/fastq/30556_thamno.fastq.gz SRP021469/fastq/30686_cyathophylla.fastq.gz SRP021469/fastq/32082_przewalskii.fastq.gz SRP021469/fastq/33413_thamno.fastq.gz ...[truncated; 13 total matched paths] -o SRP021469/TRIM/
2026-07-20 13:07:44 | INFO     | names.py             | failed to pair files, assuming data in single-end
2026-07-20 13:07:44 | INFO     | names.py             | parsed names by stripping known file suffixes
2026-07-20 13:07:44 | INFO     | names.py             | showing first 10/13 names parsed from file paths
2026-07-20 13:07:44 | INFO     | names.py             | 29154_superba          <- 29154_superba.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 30556_thamno           <- 30556_thamno.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 30686_cyathophylla     <- 30686_cyathophylla.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 32082_przewalskii      <- 32082_przewalskii.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 33413_thamno           <- 33413_thamno.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 33588_przewalskii      <- 33588_przewalskii.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 35236_rex              <- 35236_rex.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 35855_rex              <- 35855_rex.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 38362_rex              <- 38362_rex.fastq.gz
2026-07-20 13:07:44 | INFO     | names.py             | 39618_rex              <- 39618_rex.fastq.gz
2026-07-20 13:07:44 | INFO     | trim_fastqs.py       | trim input preflight found 13 usable samples and 0 skipped empty samples
[####################] 100% | Counting kmers - total jobs: 13
2026-07-20 13:07:46 | INFO     | trim_fastqs.py       | cutsite motifs set to R1=[TGCAG] at offset 0 R2=[<none>] at offset 0
2026-07-20 13:07:46 | INFO     | trim_fastqs.py       | trimming/filtering 13 samples with 'fastp' and writing to /home/deren/Documents/tools/ipyrad2/SRP021469/TRIM
2026-07-20 13:07:46 | INFO     | trim_fastqs.py       | running up to 3 parallel jobs each using up to 4 threads
[####################] 100% | Trimming - total jobs: 13
2026-07-20 13:11:47 | INFO     | trim_fastqs.py       | trimming stats written to /home/deren/Documents/tools/ipyrad2/SRP021469/TRIM/ipyrad_trim_stats_0.txt and /home/deren/Documents/tools/ipyrad2/SRP021469/TRIM/ipyrad_trim_stats_0.json

```

This log includes a few notable things to pay attention to.

First, it says `failed to pair files, assuming data in single-end`. This is expected, since our data is single-end data.

Second, you can see that it parses sample names from the file names by trimming the `.fastq.gz` suffix. The `trim` command
includes additional options for how to parse names from file names in case you want to further edit these.

Third, it says `cutsite motifs set to R1=[TGCAG] at offset 0 R2=[<none>] at offset 0`. This indicates that ipyrad2
detected the restriction enzyme motif that is left on the sequences from the restriction digestion and ligation used
during library preparation. Here "TGCAG" is the expected motif from using the enzyme PstI. It will be trimmed from all of the reads.

#### stats

When the run is finished it writes a stats file. Here we use `cat` to read it, but you can view this text file using any suitable method.

```bash
cat SRP021469/TRIM/ipyrad_trim_stats_0.txt
```

```literal
CMD: ipyrad2 trim -d SRP021469/fastq/29154_superba.fastq.gz SRP021469/fastq/30556_thamno.fastq.gz SRP021469/fastq/30686_cyathophylla.fastq.gz SRP021469/fastq/32082_przewalskii.fastq.gz SRP021469/fastq/33413_thamno.fastq.gz ...[truncated; 13 total matched paths] -o SRP021469/TRIM/

                       total_reads_before total_bases_before q20_rate_before q30_rate_before read1_mean_length_before total_reads_after total_bases_after q20_rate_after q30_rate_after read1_mean_length_after reads_filtered_by_low_quality reads_filtered_by_too_many_N reads_filtered_by_low_complexity reads_filtered_by_too_short
29154_superba                      696994           51577556        0.978960        0.950942                       74            676794          45763500       0.996413       0.975397                      67                          6354                          107                               18                       13721
30556_thamno                      1452316          107471384        0.981800        0.955849                       74           1415109          96004719       0.996578       0.976690                      67                         12213                          217                               98                       24679
30686_cyathophylla                1253109           92730066        0.967338        0.928985                       74           1173909          78490436       0.995187       0.967386                      66                         14618                          133                               28                       64421
32082_przewalskii                  964244           71354056        0.978022        0.948303                       74            936282          63238739       0.996041       0.973448                      67                          8960                          132                               63                       18807
33413_thamno                       636625           47110250        0.969651        0.934449                       74            609722          40887008       0.995692       0.969977                      67                          7494                           69                               44                       19296
33588_przewalskii                 1002923           74216302        0.979369        0.951517                       74            974319          65926297       0.996426       0.975387                      67                          9346                          122                               57                       19079
35236_rex                         1803858          133485492        0.979061        0.950356                       74           1752242         118525053       0.996323       0.974806                      67                         16711                          222                              113                       34570
35855_rex                         1409843          104328382        0.979161        0.950734                       74           1369600          92700424       0.996353       0.974826                      67                         13835                          167                               74                       26167
38362_rex                         1391175          102946950        0.980388        0.953297                       74           1353667          91632437       0.996562       0.976178                      67                         11966                          176                              125                       25241
39618_rex                          822263           60847462        0.976933        0.947744                       74            796778          53807432       0.996447       0.974718                      67                          7563                           97                               52                       17773
40578_rex                         1707942          126387708        0.981325        0.954125                       74           1667208         112853027       0.996327       0.975202                      67                         13226                          204                              106                       27198
41478_cyathophylloides            2199740          162780760        0.982123        0.956337                       74           2149301         145729619       0.996548       0.976818                      67                         17484                          322                              165                       32468
41954_cyathophylloides            2199613          162771362        0.974913        0.943468                       74           2126924         143521850       0.996085       0.972801                      67                         23268                          289                               56                       49076
```

### denovo

Next, we will assemble a denovo pseudoreference genome from the data. A close reference genome is not available for
this subclade of *Pedicularis*, and this pseudoreference will likely serve better than using a distantly related
reference genome. If you have a reference genome you can skip this step and proceed straight to `map`.

#### run

Here we specify the input data path (-d), output data path (-o), and the clustering thresholds within (-s)
and between samples (-S), and the resources to be used. Because this dataset is pretty small, we also use the
`--use-all-samples` option to build the pseudoreference from all 13 samples, instead of randomly sampling
a subset of samples, which is the default for this step.

Note that another important parameter in this step is the `-m/--min-derep-size` argument, which specifies the
minimum number of times a sequence must be observed to be included in this step. Setting a lower value will retain
more data, but increases run times. I suggest starting with the default `-m=5` setting, but if you do not recover
many thousands of loci in your denovo pseudoreference then you may want to consider rerunning it with a lower setting.

Another way to achieve speed improvements in this step is to use `--no-alignment` to skip alignment, which will select
the first/longest sequence in each component rather than inferring a consensus.

In this example run on my laptop the step takes about 40 minutes.

```bash
ipyrad2 denovo \
  -d SRP021469/TRIM/*.fastq.gz \
  -o SRP021469/DENOVO/ \
  -s 0.94 \
  -S 0.85 \
  -m 5 \
  --use-all-samples \
  -c 12 -t 3
```

#### logging
```literal
2026-07-20 13:53:56 | INFO     | cli_main.py          | ------------------------------------------------------------
2026-07-20 13:53:56 | INFO     | cli_main.py          | ----- ipyrad2 denovo: construct locus reference library -----
2026-07-20 13:53:56 | INFO     | cli_main.py          | ------------------------------------------------------------
2026-07-20 13:53:56 | INFO     | cli_main.py          | CMD: ipyrad2 denovo -d SRP021469/TRIM/29154_superba.trimmed.fastq.gz SRP021469/TRIM/30556_thamno.trimmed.fastq.gz SRP021469/TRIM/30686_cyathophylla.trimmed.fastq.gz SRP021469/TRIM/32082_przewalskii.trimmed.fastq.gz SRP021469/TRIM/33413_thamno.trimmed.fastq.gz ...[truncated; 13 total matched paths] -o SRP021469/DENOVO/ -s 0.94 -S 0.85 --use-all-samples -c 12 -t 3 -f
2026-07-20 13:53:56 | INFO     | denovo.py            | loading FASTQ inputs
2026-07-20 13:53:56 | INFO     | names.py             | failed to pair files, assuming data in single-end
2026-07-20 13:53:56 | INFO     | names.py             | parsed names by stripping known file suffixes
2026-07-20 13:53:56 | INFO     | names.py             | showing first 10/13 names parsed from file paths
2026-07-20 13:53:56 | INFO     | names.py             | 29154_superba.trimmed          <- 29154_superba.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 30556_thamno.trimmed           <- 30556_thamno.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 30686_cyathophylla.trimmed     <- 30686_cyathophylla.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 32082_przewalskii.trimmed      <- 32082_przewalskii.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 33413_thamno.trimmed           <- 33413_thamno.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 33588_przewalskii.trimmed      <- 33588_przewalskii.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 35236_rex.trimmed              <- 35236_rex.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 35855_rex.trimmed              <- 35855_rex.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 38362_rex.trimmed              <- 38362_rex.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | 39618_rex.trimmed              <- 39618_rex.trimmed.fastq.gz
2026-07-20 13:53:56 | INFO     | names.py             | normalized 13 parsed FASTQ sample name(s) by stripping recognized workflow suffixes
2026-07-20 13:53:56 | INFO     | names.py             | 29154_superba.trimmed -> 29154_superba
2026-07-20 13:53:56 | INFO     | names.py             | 30556_thamno.trimmed -> 30556_thamno
2026-07-20 13:53:56 | INFO     | names.py             | 30686_cyathophylla.trimmed -> 30686_cyathophylla
2026-07-20 13:53:56 | INFO     | names.py             | 32082_przewalskii.trimmed -> 32082_przewalskii
2026-07-20 13:53:56 | INFO     | names.py             | 33413_thamno.trimmed -> 33413_thamno
2026-07-20 13:53:56 | INFO     | names.py             | 33588_przewalskii.trimmed -> 33588_przewalskii
2026-07-20 13:53:56 | INFO     | names.py             | 35236_rex.trimmed -> 35236_rex
2026-07-20 13:53:56 | INFO     | names.py             | 35855_rex.trimmed -> 35855_rex
2026-07-20 13:53:56 | INFO     | names.py             | 38362_rex.trimmed -> 38362_rex
2026-07-20 13:53:56 | INFO     | names.py             | 39618_rex.trimmed -> 39618_rex
2026-07-20 13:53:56 | INFO     | names.py             | 40578_rex.trimmed -> 40578_rex
2026-07-20 13:53:56 | INFO     | names.py             | 41478_cyathophylloides.trimmed -> 41478_cyathophylloides
2026-07-20 13:53:56 | INFO     | names.py             | 41954_cyathophylloides.trimmed -> 41954_cyathophylloides
2026-07-20 13:53:56 | INFO     | denovo.py            | loaded 13 denovo input samples
2026-07-20 13:53:56 | INFO     | denovo.py            | selecting denovo samples
2026-07-20 13:53:56 | INFO     | denovo.py            | using all 13 denovo input samples
2026-07-20 13:53:56 | INFO     | denovo.py            | clustering within samples
[####################] 100% | Clustering within samples - total jobs: 13
2026-07-20 13:56:08 | INFO     | denovo.py            | within-sample clustering complete for 13 selected samples
2026-07-20 13:56:08 | INFO     | denovo.py            | combining per-sample summaries
2026-07-20 13:56:13 | INFO     | denovo.py            | combined 644744 consensus records across 13 selected samples
2026-07-20 13:56:13 | INFO     | denovo.py            | clustering consensus sequences across samples
[####################] 100% | Across-sample clustering
2026-07-20 13:58:04 | INFO     | denovo.py            | building denovo locus tables
[####################] 100% | Splitting global clusters - total jobs: 181327
2026-07-20 14:01:16 | WARNING  | graph.py             | retaining 12 raw oversize clusters as unsplit placeholder loci (limit=130 raw nodes; total_oversize_nodes=4255)
2026-07-20 14:01:16 | INFO     | graph.py             | built 203941 loci from 181327 graph components (rescued oversize: 0, raw placeholders: 12, post-contraction placeholders: 0)
2026-07-20 14:01:17 | INFO     | denovo.py            | building denovo reference (MAFFT)
[####################] 100% | Aligning loci - total jobs: 203941
2026-07-20 14:42:54 | INFO     | align.py             | wrote denovo reference
2026-07-20 14:42:54 | INFO     | denovo.py            | collecting final denovo QC
2026-07-20 14:42:59 | INFO     | denovo.py            | writing denovo summary report
2026-07-20 14:42:59 | INFO     | denovo.py            | wrote denovo summary report
2026-07-20 14:42:59 | INFO     | denovo.py            | denovo complete; outputs written to /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO
```

[add description of some of the logged info...]

#### pseudoreference FASTA

Below we peek at the first few lines of the pseudoreference FASTA file that is the main output of the `denovo` step.

```bash
head -n 30 SRP021469/DENOVO/denovo_reference.fa
```

```literal
>locus_1_1
ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTGCAG
>locus_1_2
ATTTGTTCTTTTTCCTAATTAAAGATCAGCCCCTGGCTATGTGTTTTCACATCGAGAATTATTTTTAGA
>locus_1_3
ATCTGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCAGTGATACGCTGCAA
>locus_2_1
AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
>locus_2_2
AAAGAAAAACAGCAAAATCCAATCCAATTTATCGTAATCGATTCCTTAACATGTTGGTTAACCGTATTCT
>locus_3_1
GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACCGAG
>locus_3_2
TCCCAAATGGATTGGCTTATTCCAAAAAGACCTTGTTCTTGAAAGATGTATCTCGTGTCTAGTACTAAG
>locus_4_1
TCGCTGCCCAGAAAGAATGATGTTGGTTTCGATGTTGGCAACATCAGACTCAACATCAATGTCTACTGCTCAT
>locus_4_2
TTGCTGCCCAGAATGATGTTGGTATCGATGTTGGCAACATCAGACTCAACATCGATGTCTGCTGCTCCA
>locus_4_3
TTGCTTCCTAGATTGATGTTGGTATCGATGTTGGCAACATCAGAATCAACATTGATGTCTGCTGCTCCA
>locus_4_4
TTGCTGCCCAGAATGATGTTGGTATCGATGTTGGCAACATCAGACTCAACATCGATGTCTGCAATTTTC
>locus_4_5
CTGCTGCCTATAAAAATGATGTTGGTTAAGATATTGGCAACATCAGACTCAACATCGATGTCTACTGCT
>locus_4_6
TTGCTGCCTAGAATGATGTTGTTATTGATGTTGGCAACATCAGACTCAACATCGATGTTTGTTGCTCCA
>locus_4_7
TCGCTGCGTAGAAAATTGATGTTGGTTTAGATGTTTGCAACATCAGACTCAACATCAATGTCTGTTGCT
>locus_4_8
CTGCTGCCCAAAAAAATGATGTTGGTTCCGATGTTGGCAACATTAGACTCAACATCGTTTTCTGTTGCT
...
```

Here the pseudoreference is composed of >200K loci that were recovered by graph clustering and splitting during
the denovo assembly process. The process involved first clustering within each sample to collapse similar reads into
a consensus sequence (e.g., using 0.95 threshold), and clustering across samples at a lower threshold (0.85) to cluster
homologs across samples. The resulting graph included many clusters that contained duplicated regions of the genome from
a single sample, i.e., paralogs. These are split using a graph splitting algorithm to find graph components that contain
at most one sequence per sample (i.e., no duplications).

You can see the result of this process in the names of the loci in the pseudoreference. For example,
given the consensus sequences that grouped into the first cluster (`locus_1`), the graph splitting algorithm further split
this into three distinct components (`locus_1_1`, `locus_1_2`, and `locus_1_3`).
The sequences of these three loci look pretty similar, but clearly differ in their sequences.
We purposely keep all three paralogous copies in the final pseudoreference as these will later allow
for reads to map best to one locus versus another.

#### stats

<!-- The full stats file for denovo is quite long, you can click below to expand it and see the full thing. -->

This dataset recovered 640,636 consensus sequences across the 13 samples, which were clustered and split
into a final set of 203,941 pseudoreference loci. Graph splitting was applied to 36,668 loci
that contained duplications. After splitting, 3,998 loci contained sequences from all 13 samples;
35,423 contained more than half the samples; 89,687 contained more than two samples; and
114,254 contained only a single sample.

<!-- ??? note "click to expand full stats file for denovo run" -->

```bash
cat SRP021469/DENOVO/denovo.stats.txt
```

```literal
CMD: ipyrad2 denovo -d SRP021469/TRIM/29154_superba.trimmed.fastq.gz SRP021469/TRIM/30556_thamno.trimmed.fastq.gz SRP021469/TRIM/30686_cyathophylla.trimmed.fastq.gz SRP021469/TRIM/32082_przewalskii.trimmed.fastq.gz SRP021469/TRIM/33413_thamno.trimmed.fastq.gz ...[truncated; 13 total matched paths] -o SRP021469/DENOVO/ -s 0.94 -S 0.85 --use-all-samples -c 12 -t 3 -f

# Inputs
FASTQ files            13
Selected samples       13
Total input samples    13
Sample selection mode  all
Read layout            single-end

# Clustering Parameters
Within-sample similarity        0.940000
Across-sample similarity        0.850000
Minimum VSEARCH query coverage  0.750000
Minimum dereplication size      5
Minimum read length             35
Minimum merge overlap           20
Maximum merge differences       4
Allow reverse complement        False

# Denovo Summary
Consensus records                     640,636
Loci written                          203,941
Single-sequence loci                  113,659
Identical-sequence loci               18,679
Loci requiring MAFFT                  71,603
Joined-spacer loci                    0
Mixed reconciled spacer loci          1,836
Spacer-stripped output loci           202,105
Duplicated components seen            14,054
Same-sample reconciliation attempted  14,042
Components reconciled                 1,694
Joined-only reconciled loci           0
Mixed reconciled loci                 1,836
Mixed reconciled groups               2,276

# Locus QC
Singleton loci                           114,254
Singleton locus fraction                 0.560231
Loci with 2+ samples                     89,687
Loci with half or more selected samples  35,423
Loci with all selected samples           3,998
Mean samples per locus                   3.125
Median samples per locus                 1.000
Maximum samples per locus                13
Mean cores per locus                     3.141
Median cores per locus                   1.000
Maximum cores per locus                  39
Multi-core single-sample loci            595
Duplicated-component loci                36,668
Reconciled loci                          6,294

# Component QC
Audited components           14,054
Processed components         14,042
Oversize unsplit components  12
Largest component nodes      1,052

# Component Node Summary
Quantile  Input nodes  Contracted nodes
p50       9            8
p90       19           19
p99       40           37
max       1,052        126

# Selected Sample Summary
Sample                  Consensus records  Read count  Joined records  Merged records  Single records
41954_cyathophylloides  78,525             1,570,430   0               0               78,525
40578_rex               55,687             1,244,939   0               0               55,687
41478_cyathophylloides  54,705             1,762,638   0               0               54,705
35236_rex               54,492             1,059,418   0               0               54,492
35855_rex               53,657             937,397     0               0               53,657
38362_rex               52,475             1,078,200   0               0               52,475
30556_thamno            50,756             976,319     0               0               50,756
30686_cyathophylla      50,687             711,961     0               0               50,687
33588_przewalskii       45,808             631,242     0               0               45,808
39618_rex               43,082             514,105     0               0               43,082
32082_przewalskii       41,656             598,909     0               0               41,656
29154_superba           35,145             427,815     0               0               35,145
33413_thamno            28,069             287,939     0               0               28,069

# Locus Occupancy
Samples with data  Loci     Fraction of final loci
0                  0        0.000000
1                  114,254  0.560231
2                  29,320   0.143767
3                  9,123    0.044734
4                  6,337    0.031073
5                  4,925    0.024149
6                  4,559    0.022355
7                  4,031    0.019766
8                  3,939    0.019314
9                  4,896    0.024007
10                 6,092    0.029871
11                 6,577    0.032250
12                 5,890    0.028881
13                 3,998    0.019604

# Runtime
Cores                     12
VSEARCH threads per job   3
VSEARCH worker processes  4
MAFFT threads per job     1
MAFFT worker processes    12
Alignment mode            mafft
MAFFT timeout (seconds)   900
Keep intermediates        False

# Outputs
Reference FASTA       /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo_reference.fa
Locus mapping table   /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.loci.mapping.tsv
Locus stats table     /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.loci.stats.tsv
Sample graph summary  /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.sample_graph_summary.tsv
Run summary report    /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.stats.txt
Run stats json        /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.stats.json
Audit directory       /home/deren/Documents/tools/ipyrad2/SRP021469/DENOVO/denovo.audit
Intermediate files    cleaned on success
```


### map

Next, we will map reads from each sample to the pseudoreference genome to generate BAM alignment files.

#### run

In the `map` command we specify the trimmed fastq files as the data input (-d), the pseudoreference genome fasta
as the reference input (-r), and specify a path to write the output files (-o). We also specify the resources to
be used, here assigning 12 cores to be distributed multiple 4-threaded jobs.

```bash
ipyrad2 map \
  -d SRP021469/TRIM/*.fastq.gz \
  -r SRP021469/DENOVO/denovo_reference.fa \
  -o SRP021469/MAP/ \
  -c 12 -t 4
```

#### logging

The logging report indicates that it successfully identified the samples and their names, completed mapping,
and calculated stats for the mapped reads.

```literal
2026-07-20 14:47:09 | INFO     | cli_main.py          | --------------------------------------------------------------
2026-07-20 14:47:09 | INFO     | cli_main.py          | ----- ipyrad2 map: map reads and write coordinate-sorted BAMs -----
2026-07-20 14:47:09 | INFO     | cli_main.py          | --------------------------------------------------------------
2026-07-20 14:47:09 | INFO     | cli_main.py          | CMD: ipyrad2 map -d SRP021469/TRIM/29154_superba.trimmed.fastq.gz SRP021469/TRIM/30556_thamno.trimmed.fastq.gz SRP021469/TRIM/30686_cyathophylla.trimmed.fastq.gz SRP021469/TRIM/32082_przewalskii.trimmed.fastq.gz SRP021469/TRIM/33413_thamno.trimmed.fastq.gz ...[truncated; 13 total matched paths] -r SRP021469/DENOVO/denovo_reference.fa -o SRP021469/MAP/ -c 12 -t 4
2026-07-20 14:47:09 | INFO     | names.py             | failed to pair files, assuming data in single-end
2026-07-20 14:47:09 | INFO     | names.py             | parsed names by stripping known file suffixes
2026-07-20 14:47:09 | INFO     | names.py             | showing first 10/13 names parsed from file paths
2026-07-20 14:47:09 | INFO     | names.py             | 29154_superba.trimmed          <- 29154_superba.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 30556_thamno.trimmed           <- 30556_thamno.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 30686_cyathophylla.trimmed     <- 30686_cyathophylla.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 32082_przewalskii.trimmed      <- 32082_przewalskii.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 33413_thamno.trimmed           <- 33413_thamno.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 33588_przewalskii.trimmed      <- 33588_przewalskii.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 35236_rex.trimmed              <- 35236_rex.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 35855_rex.trimmed              <- 35855_rex.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 38362_rex.trimmed              <- 38362_rex.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | 39618_rex.trimmed              <- 39618_rex.trimmed.fastq.gz
2026-07-20 14:47:09 | INFO     | names.py             | normalized 13 parsed FASTQ sample name(s) by stripping recognized workflow suffixes
2026-07-20 14:47:09 | INFO     | names.py             | 29154_superba.trimmed -> 29154_superba
2026-07-20 14:47:09 | INFO     | names.py             | 30556_thamno.trimmed -> 30556_thamno
2026-07-20 14:47:09 | INFO     | names.py             | 30686_cyathophylla.trimmed -> 30686_cyathophylla
2026-07-20 14:47:09 | INFO     | names.py             | 32082_przewalskii.trimmed -> 32082_przewalskii
2026-07-20 14:47:09 | INFO     | names.py             | 33413_thamno.trimmed -> 33413_thamno
2026-07-20 14:47:09 | INFO     | names.py             | 33588_przewalskii.trimmed -> 33588_przewalskii
2026-07-20 14:47:09 | INFO     | names.py             | 35236_rex.trimmed -> 35236_rex
2026-07-20 14:47:09 | INFO     | names.py             | 35855_rex.trimmed -> 35855_rex
2026-07-20 14:47:09 | INFO     | names.py             | 38362_rex.trimmed -> 38362_rex
2026-07-20 14:47:09 | INFO     | names.py             | 39618_rex.trimmed -> 39618_rex
2026-07-20 14:47:09 | INFO     | names.py             | 40578_rex.trimmed -> 40578_rex
2026-07-20 14:47:09 | INFO     | names.py             | 41478_cyathophylloides.trimmed -> 41478_cyathophylloides
2026-07-20 14:47:09 | INFO     | names.py             | 41954_cyathophylloides.trimmed -> 41954_cyathophylloides
2026-07-20 14:47:09 | INFO     | mapper.py            | indexing reference: denovo_reference.fa
2026-07-20 14:47:14 | INFO     | mapper.py            | mapping 13 samples to coordinate-sorted BAMs in /home/deren/Documents/tools/ipyrad2/SRP021469/MAP
2026-07-20 14:47:14 | INFO     | mapper.py            | using up to 12 cores (up to 3 multi-threaded jobs using 4 threads)
[####################] 100% | Mapping - total jobs: 13
[####################] 100% | Gathering mapping stats - total jobs: 13
2026-07-20 14:50:41 | INFO     | mapper.py            | mapping stats written to /home/deren/Documents/tools/ipyrad2/SRP021469/MAP/ipyrad_map_stats_0.txt and /home/deren/Documents/tools/ipyrad2/SRP021469/MAP/ipyrad_map_stats_0.json
```

#### stats

The human-readable stats file from `map` is important to read as it provides some useful guidance on the parameter settings
that should be used in the next step.

First, we can see in the `## Applied mapping summary` section that between 74-92\% of reads were successfully mapped to the pseudoreference.
However, not that this does not yet indicate how well these reads mapped. We will apply a number of filters to only keep reads that mapped uniquely and accurately to loci in the reference. The next section `## Assemble read-filter preview` shows us how the filters that will be applied in the next step (assemble) will affect the mapping. This shows that 88-94\% of reads will pass filtering, suggesting that these filters are not too strict. This is further reinforced by the next section `## Preview metric summaries` which shows that the MAPQ scores of these reads are generally very high (>50), with few soft clipped bases (<5) and low edit distances to the reference (<2).


```bash
cat ./SRP021469/MAP/ipyrad_map_stats_0.txt
```

```literal
CMD: ipyrad2 map -d SRP021469/TRIM/29154_superba.trimmed.fastq.gz SRP021469/TRIM/30556_thamno.trimmed.fastq.gz SRP021469/TRIM/30686_cyathophylla.trimmed.fastq.gz SRP021469/TRIM/32082_przewalskii.trimmed.fastq.gz SRP021469/TRIM/33413_thamno.trimmed.fastq.gz ...[truncated; 13 total matched paths] -r SRP021469/DENOVO/denovo_reference.fa -o SRP021469/MAP/ -c 12 -t 4

# ipyrad2 map stats
# Final BAMs are coordinate sorted and indexed.

## Applied mapping summary
# These counts describe filters already applied during ipyrad2 map.

                       input_reads reads_removed_unmapped_or_nonprimary reads_in_final_bam fraction_input_reads_retained_in_final_bam
sample
29154_superba               686092                               120393             565699                                      0.825
30556_thamno               1421254                               217471            1203783                                      0.847
30686_cyathophylla         1179291                               272866             906425                                      0.769
32082_przewalskii           945670                               183438             762232                                      0.806
33413_thamno                612382                               159210             453172                                      0.740
33588_przewalskii           985540                               193737             791803                                      0.803
35236_rex                  1761279                               521219            1240060                                      0.704
35855_rex                  1377005                               182748            1194257                                      0.867
38362_rex                  1361839                               104941            1256898                                      0.923
39618_rex                   802399                               121029             681370                                      0.849
40578_rex                  1674280                               220551            1453729                                      0.868
41478_cyathophylloides     2175418                               192717            1982701                                      0.911
41954_cyathophylloides     2137180                               321352            1815828                                      0.850

## Assemble read-filter preview (not applied during mapping)
# These preview thresholds were not applied during mapping.
# Use them to guide ipyrad2 assemble read filters: -qm/--min-map-q, -ms/--max-softclip, -me/--max-nm.
# Preview mode: read-level thresholds evaluated on final BAM reads.
# MAPQ threshold: 20
# Soft-clipped bases threshold: 25
# NM threshold: 50

### Preview filter effects
                       reads_failing_min_mapq_20 reads_failing_max_softclip_25 reads_failing_max_nm_50 reads_passing_all_preview_filters fraction_reads_passing_all_preview_filters
sample
29154_superba                              48936                         20797                       0                            505581                                      0.894
30556_thamno                               79936                         22846                       0                           1113404                                      0.925
30686_cyathophylla                         51496                         19823                       0                            845628                                      0.933
32082_przewalskii                          80174                         31606                       0                            664465                                      0.872
33413_thamno                               28442                          9780                       0                            419778                                      0.926
33588_przewalskii                          75783                         31944                       0                            698451                                      0.882
35236_rex                                  69843                         21527                       0                           1159361                                      0.935
35855_rex                                  79813                         21062                       0                           1104175                                      0.925
38362_rex                                  73497                         17725                       0                           1174805                                      0.935
39618_rex                                  47389                         13921                       0                            627166                                      0.920
40578_rex                                  73232                         18845                       0                           1372330                                      0.944
41478_cyathophylloides                    137127                         27238                       0                           1833296                                      0.925
41954_cyathophylloides                    105090                         16614                       0                           1703374                                      0.938

### Preview metric summaries
                       mapq_mean mapq_median mapq_stdev softclip_mean softclip_median softclip_stdev nm_mean nm_median nm_stdev
sample
29154_superba             52.410      60.000     16.373         2.233           0.000          7.403   1.331     1.000    1.851
30556_thamno              53.898      60.000     14.574         1.389           0.000          5.960   1.064     0.000    1.665
30686_cyathophylla        54.543      60.000     13.868         1.524           0.000          6.298   1.149     0.000    1.765
32082_przewalskii         50.736      60.000     17.452         2.980           0.000          8.452   1.583     1.000    1.951
33413_thamno              53.908      60.000     14.431         1.507           0.000          6.236   1.177     1.000    1.764
33588_przewalskii         51.283      60.000     16.984         2.997           0.000          8.434   1.629     1.000    1.968
35236_rex                 54.550      60.000     13.766         1.226           0.000          5.598   1.056     0.000    1.631
35855_rex                 53.500      60.000     14.806         1.306           0.000          5.705   1.175     1.000    1.707
38362_rex                 54.364      60.000     13.909         1.072           0.000          5.098   0.981     0.000    1.582
39618_rex                 53.520      60.000     15.005         1.421           0.000          6.049   1.134     1.000    1.689
40578_rex                 54.799      60.000     13.183         0.997           0.000          4.921   0.933     0.000    1.572
41478_cyathophylloides    53.779      60.000     14.946         1.446           0.000          5.615   1.184     1.000    1.739
41954_cyathophylloides    54.373      60.000     13.999         0.910           0.000          4.532   0.934     0.000    1.503
```

### assemble

Finally, the `assemble` step represents the main step of the ipyrad2 assembly workflow. Here we apply filters to the BAM
alignments to keep only confidently mapped reads, which are then used to delimit RAD loci with sufficiently high coverage
across samples, make variant calls, filter for paralogy, and write the final locus alignments into a database and several
output files. A verbose stats file is also produced.

#### run

Here we indicate one or more paths to BAM files as the data input (-d); a path to the reference FASTA (-r);
a path to store the outputs (-o); and a prefix name for the results (-n). We also show two of the most commonly
changed parameter settings, `-m/--min-locus-sample-coverage`, `-s/--min-sample-depth`, and `-qm/--min-map-q`.
Respectively, these parameters effect the minimum number of samples that must be present in a locus to be kept;
the minimum read coverage that must be present at a site in a sample to make a variant call;
and the minimum mapping score of a read to be retained.

```bash
ipyrad2 assemble \
  -d SRP021469/MAP/*.bam \
  -r SRP021469/DENOVO/denovo_reference.fa \
  -o SRP021469/OUT \
  -n assembly \
  -m 4 \
  -s 5 \
  -qm 40 \
  -c 12 -t 4
```

#### stats

The stats

```literal
CMD: ipyrad2 assemble -d SRP021469/MAP/29154_superba.trimmed.sorted.bam SRP021469/MAP/30556_thamno.trimmed.sorted.bam SRP021469/MAP/30686_cyathophylla.trimmed.sorted.bam SRP021469/MAP/32082_przewalskii.trimmed.sorted.bam ...[truncated; 13 total matched paths] -r SRP021469/DENOVO/denovo_reference.fa -o SRP021469/OUT -n assembly -qm 40 -c 12 -t 4

# Assemble Summary
Samples                                               13
Shared loci before minimum sample coverage filter     199,531
Shared loci after delimiting                          49,417
Shared loci after paralog filtering                   49,378
Final loci written                                    45,448
Final loci retained fraction after paralog filtering  0.920410
Final loci retained fraction after delimiting         0.919684
Assembled sites                                       3,100,317
Final SNP sites written                               204,914
Variable sites                                        183,882
Phylogenetically informative sites                    77,117
Alignment matrix occupancy fraction                   0.670662
Overlapping indel clusters masked                     1,728
Overlapping indel records removed                     4,037
Overlapping indel bases masked                        11,543

# Locus Filtering
Loci filtered by minimum length                 119
Loci filtered by minimum sample coverage        150,114
Loci filtered by maximum variant frequency      0
Loci filtered by maximum shared heterozygosity  3,906
Loci filtered by maximum depth outlier          0

# Sample Masking
Loci with samples masked by minimum observed fraction threshold  4
Sample masks triggered by minimum observed fraction threshold    4
Loci with samples masked by sample heterozygosity threshold      135
Sample masks triggered by sample heterozygosity threshold        153

# Alignment Summary
Mean locus length                           68.217
Median locus length                         69.000
Minimum locus length                        25
Maximum locus length                        89
Mean samples per locus                      8.842
Median samples per locus                    9.000
Sites with sample coverage >= 2             3,083,388
Sites with sample coverage >= 3             3,075,143
Sites with sample coverage >= 4             3,068,712
Sites with sample coverage >= trim minimum  3,068,712

# Sample Summary
Sample                  Sample type  Read layout  Reads before filtering  Reads after filtering  Loci in alignment  Loci fraction in alignment  Shared loci with nonzero depth  Shared-depth loci fraction  Mean depth in shared loci  Median depth in shared loci  Mean depth in nonzero shared loci  Median depth in nonzero shared loci  Masked by minimum observed fraction threshold  Masked by sample heterozygosity threshold
29154_superba           RAD          SE           565,699                 479,661                26,106             0.574415                    26,120                          0.574723                    7.216                      5.000                        12.556                             8.000                                0                                              14
30556_thamno            RAD          SE           1,203,783               1,054,723              35,767             0.786987                    35,771                          0.787075                    16.652                     11.913                       21.157                             13.928                               0                                              4
30686_cyathophylla      RAD          SE           906,425                 805,755                30,766             0.676949                    30,776                          0.677170                    10.366                     8.000                        15.307                             10.725                               0                                              10
32082_przewalskii       RAD          SE           762,232                 613,739                16,664             0.366661                    16,693                          0.367299                    6.118                      0.000                        16.656                             9.000                                0                                              29
33413_thamno            RAD          SE           453,172                 396,835                24,518             0.539474                    24,528                          0.539694                    5.696                      5.000                        10.553                             7.000                                0                                              10
33588_przewalskii       RAD          SE           791,803                 646,001                18,695             0.411349                    18,718                          0.411855                    6.138                      0.000                        14.904                             10.028                               0                                              23
35236_rex               RAD          SE           1,240,060               1,101,371              36,994             0.813985                    37,003                          0.814183                    18.193                     15.725                       22.345                             17.884                               1                                              8
35855_rex               RAD          SE           1,194,257               1,039,563              37,436             0.823711                    37,445                          0.823909                    15.646                     12.493                       18.991                             14.000                               0                                              9
38362_rex               RAD          SE           1,256,898               1,117,068              37,349             0.821796                    37,361                          0.822060                    18.421                     14.000                       22.409                             15.967                               2                                              10
39618_rex               RAD          SE           681,370                 590,414                33,200             0.730505                    33,206                          0.730637                    8.719                      7.768                        11.933                             9.000                                1                                              5
40578_rex               RAD          SE           1,453,729               1,311,633              38,000             0.836120                    38,011                          0.836362                    22.080                     15.000                       26.400                             17.000                               0                                              11
41478_cyathophylloides  RAD          SE           1,982,701               1,731,466              34,228             0.753124                    34,239                          0.753366                    27.202                     23.000                       36.107                             27.000                               0                                              11
41954_cyathophylloides  RAD          SE           1,815,828               1,603,428              32,140             0.707182                    32,149                          0.707380                    11.821                     9.000                        16.711                             11.435                               0                                              9

# Locus Occupancy
Samples with data  RAD loci before min sample coverage  RAD loci after min sample coverage  Final filtered RAD loci with WGS  Cumulative final loci  Fraction of final loci
0                  0                                    0                                   0                                 0                      0.000000
1                  106,886                              0                                   0                                 0                      0.000000
2                  33,836                               0                                   0                                 0                      0.000000
3                  9,363                                0                                   0                                 0                      0.000000
4                  6,260                                6,232                               5,529                             5,529                  0.121656
5                  4,605                                4,604                               3,842                             9,371                  0.084536
6                  4,220                                4,220                               3,462                             12,833                 0.076175
7                  3,643                                3,643                               3,274                             16,107                 0.072038
8                  3,262                                3,262                               2,999                             19,106                 0.065988
9                  4,109                                4,110                               3,859                             22,965                 0.084910
10                 5,383                                5,382                               5,225                             28,190                 0.114967
11                 6,602                                6,602                               6,377                             34,567                 0.140314
12                 5,872                                5,872                               5,726                             40,293                 0.125990
13                 5,490                                5,490                               5,155                             45,448                 0.113426

```

#### loci

```bash
head -n 300 SRP021469/OUT/assembly.stats.txt
```

```literal
assembly_reference_sequence  ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTGCAG
29154_superba                ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
30556_thamno                 ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
30686_cyathophylla           ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
32082_przewalskii            ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
33413_thamno                 ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
33588_przewalskii            ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
35236_rex                    ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
35855_rex                    ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
38362_rex                    ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
39618_rex                    ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
40578_rex                    ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
41478_cyathophylloides       ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
41954_cyathophylloides       ATCCGCTCTTTTTCCTATTCAAAGATCAGCCCCCTGGCTCTGTGTTTTCACATCGAGAATTATTTNCAG
//                                                                                                |0:locus_1_1:1-69
assembly_reference_sequence  AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
29154_superba                AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
30556_thamno                 AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
30686_cyathophylla           AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
32082_przewalskii            AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
33413_thamno                 AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
33588_przewalskii            AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
35236_rex                    AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
35855_rex                    AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
38362_rex                    AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
39618_rex                    AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
40578_rex                    AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
41478_cyathophylloides       AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
41954_cyathophylloides       AAGAAAAAACAGCAAAATCCGATCCAATTTATCGTAATCGATTAGTTAACATGTTGGTTAACCGTATTC
//                                                                                                |1:locus_2_1:1-69
assembly_reference_sequence  GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACC
29154_superba                GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
30556_thamno                 GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
30686_cyathophylla           GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACC
32082_przewalskii            GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
33413_thamno                 GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
33588_przewalskii            GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
35236_rex                    GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACC
35855_rex                    GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
38362_rex                    GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACC
39618_rex                    GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
40578_rex                    GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
41478_cyathophylloides       GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACN
41954_cyathophylloides       GGTCCCAAATGAATTGGCTTATTCGAAAAAGGCCTTGTTCTTTGGAAGATCTATCTCGTGTCTGGTACT
//                                                                                               *|2:locus_3_1:1-69
assembly_reference_sequence  CCTAGAAATGATGTTGGCTTCGATGTTGGCAACATCATACTCCACATCAATGTCTGCTGCTCC
30556_thamno                 CCTAGAAATGATGTTGGCTTCGATNTTGGCAACATCATACTCCACATCAATGTCTNCTGCTCC
35236_rex                    CCTAGAAATGATGTTGGCTTCGATNTTGGCAACATCATACTCCAYATCAATGTCTNCTGCTCC
35855_rex                    CCTAGAAATGATGTTGGCTTCGATNTTGGCAACATCATACTCCACATCAATGTCTNCTGCTCC
38362_rex                    CCTAGAAATGATGTTGGCTTCGATNTTGGMAACATCATACTCCACATCAATGTCWNCTGCTCC
39618_rex                    CCTAGAAATGATGTTGGCTTCGATNTTGGCAACATCATACTCCACATCAATGTCTNCTGCTCC
40578_rex                    CCTAGAAATGATGTTGGCTTCGATNTTGGCAACATCATACTCCACATCAATGTCTNCTGCTCC
//                                                        -              -         -        |3:locus_4_25:9-71
assembly_reference_sequence  TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
29154_superba                TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
30556_thamno                 TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
30686_cyathophylla           TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
32082_przewalskii            TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATCAC
33413_thamno                 TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
33588_przewalskii            TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATCAC
35236_rex                    TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
35855_rex                    TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
38362_rex                    TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
39618_rex                    TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
40578_rex                    TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
41478_cyathophylloides       TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
41954_cyathophylloides       TACCTCGACGTGACATGAGCGTGAAAGGGGTTTAAGAATCAGTTTTCTTTTTATAAGGGCTAAAATTAC
//                                                                                             *  |4:locus_5_1:1-69
assembly_reference_sequence  CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
29154_superba                CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
30556_thamno                 CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
30686_cyathophylla           CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
33413_thamno                 CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
33588_przewalskii            CTGCACGAGCCCTTCCTGCATGCCACAAATGACCTACGAAGAAGAAGAATCCTAGAACAAAATGAGAGG
35236_rex                    CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
35855_rex                    CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
38362_rex                    CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
39618_rex                    CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
40578_rex                    CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
41478_cyathophylloides       CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
41954_cyathophylloides       CTGCACGAGCCCTTCCCGCATGCCACAAATGACCTACGAATAAGAAGAATCCTAGAACAAAATGAGAGG
//                                           -                       -                            |5:locus_6_1:1-69
assembly_reference_sequence  TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
29154_superba                TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGACTTACTAGCCTTGATCGTT
30556_thamno                 TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
30686_cyathophylla           TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
32082_przewalskii            TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGACTTACTAGCCTTGATCGTT
33413_thamno                 TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
33588_przewalskii            TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGACTTACTAGCCTTGATCGTT
35236_rex                    TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
35855_rex                    TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
38362_rex                    TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
39618_rex                    TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
40578_rex                    TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGGCTTACTAGCCTTGATCGTT
41478_cyathophylloides       TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGACTTACTAGCCTTGATCGTT
41954_cyathophylloides       TAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACCGATGGACTTACTAGCCTTGATCGTT
//                                                                            *                   |6:locus_7_1:1-69
assembly_reference_sequence  CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
29154_superba                CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
30556_thamno                 CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
30686_cyathophylla           CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
32082_przewalskii            CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
33413_thamno                 CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
33588_przewalskii            CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
35236_rex                    CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
35855_rex                    CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
38362_rex                    CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
39618_rex                    CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
40578_rex                    CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
41478_cyathophylloides       CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
41954_cyathophylloides       CCGGATTTGAAAAAGGAATTGATCGCGATTTTGAACCTGTTCTTTCCATGACCCCTCTTAATTGAGATG
//                                                                                                |7:locus_8_1:1-69
assembly_reference_sequence  CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
29154_superba                CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
30556_thamno                 CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
30686_cyathophylla           CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
32082_przewalskii            CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
33413_thamno                 CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
33588_przewalskii            CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
35236_rex                    CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
35855_rex                    CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
38362_rex                    CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
39618_rex                    CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
40578_rex                    CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
41478_cyathophylloides       CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
41954_cyathophylloides       CCCCTGCTTCTTCAGGCGGAACTCCAGGTTGAGGAGTTACTCGGAATGCTGCCAAGATATCA
//                                                                                         |8:locus_9_1:1-62
assembly_reference_sequence  TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
29154_superba                TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
30556_thamno                 TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
30686_cyathophylla           TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
32082_przewalskii            TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
33413_thamno                 TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
33588_przewalskii            TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
35236_rex                    TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
35855_rex                    TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
38362_rex                    TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
39618_rex                    TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
40578_rex                    TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
41478_cyathophylloides       TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
41954_cyathophylloides       TTGTGATTGATCAAGAAGGAAATCCAAAAGGAACTCGCATTTTTGGTGCAATCCCGCGGGAATTGCGAC
//                                                                                                |9:locus_10_1:1-69
assembly_reference_sequence  AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
29154_superba                AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCGTAAAGTAGAGATCCAGAAACAGGTTCACGAA
30556_thamno                 AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
30686_cyathophylla           AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
32082_przewalskii            AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCGTAAAGTAGAGATCCAGAAACAGGTTCACGAA
33413_thamno                 AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
33588_przewalskii            AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCGTAAAGTAGAGATCCAGAAACAGGTTCACGAA
35236_rex                    AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
35855_rex                    AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
38362_rex                    AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
39618_rex                    AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
40578_rex                    AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCATAAAGTAGAGATCCAGAAACAGGTTCACGAA
41478_cyathophylloides       AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCGTAAAGTAGAGATCCAGAAACAGGTTCACGAA
41954_cyathophylloides       AAGTAGGAATAATGGCACCCGAGATAATATTGTTTCCGTAAAGTAGAGATCCAGAAACAGGTTCACGAA
//                                                                *                               |10:locus_11_1:1-69
assembly_reference_sequence  CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
29154_superba                CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
30556_thamno                 CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
30686_cyathophylla           CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
32082_przewalskii            CATTGTCATCATATCGTATTATCATGCCGCTGTTACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
33413_thamno                 CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
33588_przewalskii            CATTGTCATCATATCGTATTATCATGCCGCTGTTACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
35236_rex                    CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
35855_rex                    CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
38362_rex                    CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
39618_rex                    CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
40578_rex                    CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
41478_cyathophylloides       CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
41954_cyathophylloides       CATTGTCATCATATCGTATTATCATGCCGCTGTCACGTTTAAGTTCTTTACAGGTACGAACAATGACAG
//                                                            *                                   |11:locus_13_1:1-69
assembly_reference_sequence  CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
29154_superba                CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
30556_thamno                 CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
30686_cyathophylla           CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
32082_przewalskii            CAGAATAAACCAATTTAAAAATGGGATAACATGCTTTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
33413_thamno                 CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
33588_przewalskii            CAGAATAAACCAATTTAAAAATGGGATAACATGCTTTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
35236_rex                    CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
35855_rex                    CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
38362_rex                    CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
39618_rex                    CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
40578_rex                    CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
41478_cyathophylloides       CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
41954_cyathophylloides       CAGAATAAACCAATTTAAAAATGGGATAACATGCTCTATAGGGCATGAGCTCGAGTATCATAAGTGTTT
//                                                              *                                 |12:locus_14_1:1-69
assembly_reference_sequence  GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
29154_superba                GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
30556_thamno                 GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
30686_cyathophylla           GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
32082_przewalskii            GCCCTAGAACCAGAAATAAATAAGCTTATTCTTTGTTCACTTGAATCAGAATTCTAACCCGAACTCAAA
33413_thamno                 GCCTTAGAACTAGAAATAAATAAGCTTATTTTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
33588_przewalskii            GCCCTAGAACCAGAAATAAATAAGCTTATTCTTTGTTCACTTGAATCAGAATTCTAACCCGAACTCAAA
35236_rex                    GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
35855_rex                    GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
38362_rex                    GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
39618_rex                    GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
40578_rex                    GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
41478_cyathophylloides       GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
41954_cyathophylloides       GCCTTAGAACTAGAAATAAATAAGCTTATTCTTTGTTCATTTGAATCAGAATTCCAACCCGAACTCAAA
//                              *      *                   -        *              *              |13:locus_15_1:1-69
assembly_reference_sequence  AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
29154_superba                AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
30556_thamno                 AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
30686_cyathophylla           AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
32082_przewalskii            AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGGTCTGGTTCACCCGACACGTACAC
33413_thamno                 AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
33588_przewalskii            AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGGTCTGGTTCACCCGACACGTACAC
35236_rex                    AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
35855_rex                    AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
38362_rex                    AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
39618_rex                    AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
40578_rex                    AACACGTGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
41478_cyathophylloides       AACACGCGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
41954_cyathophylloides       AACACGCGCGAGCCCCTTCGAATGGAAAAATAAAATTTAATGAGGATCTGGTTCACCCGACACGTACAC
//                                 *                                      *                       |14:locus_16_1:1-69
assembly_reference_sequence  CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
29154_superba                CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
30556_thamno                 CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
30686_cyathophylla           CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
32082_przewalskii            CTCTAGTGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
33413_thamno                 CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
33588_przewalskii            CTCTAGTGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
35236_rex                    CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
35855_rex                    CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
38362_rex                    CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
39618_rex                    CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
40578_rex                    CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
41478_cyathophylloides       CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
41954_cyathophylloides       CTCTAGCGGTCGGCTCGGTTCTTTCAAATTGTTTCTCATTATTGAGAAAAGGTAACAAAGATAAAATAC
//                                 *                                                              |15:locus_17_1:1-69
assembly_reference_sequence  TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
29154_superba                TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
30556_thamno                 TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
30686_cyathophylla           TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
32082_przewalskii            TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
33413_thamno                 TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
33588_przewalskii            TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
35236_rex                    TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
35855_rex                    TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
38362_rex                    TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
39618_rex                    TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
40578_rex                    TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
41478_cyathophylloides       TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
41954_cyathophylloides       TTTGAGCAGCAAAGGGTGTCCCTCTTCTCGTACCCTTGAATCCACAAGTACCGGCCGAGGACC
//                                                                                          |16:locus_18_1:1-63
assembly_reference_sequence  TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
29154_superba                TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
30556_thamno                 TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
30686_cyathophylla           TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
32082_przewalskii            TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
33413_thamno                 TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
33588_przewalskii            TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
35236_rex                    TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
35855_rex                    TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
38362_rex                    TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
39618_rex                    TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
40578_rex                    TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
41478_cyathophylloides       TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
41954_cyathophylloides       TACCCCCCGTGAATACTCCGCCGGTATGAAAAGTTCTTAATGTTAATTGAGTGCCCGGTTCTCCAATTG
//                                                                                                |17:locus_19_1:1-69
assembly_reference_sequence  CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
29154_superba                CAAACGCAACTTTCAGTATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
30556_thamno                 CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
30686_cyathophylla           CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
32082_przewalskii            CAAACGCAGCTTTCAATATGGGTTCCGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
33413_thamno                 CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
33588_przewalskii            CAAACGCAGCTTTCAATATGGGTTCCGACGAAGCCAATTTAGTAATTAGTAAAGMTGAGGTTAATGAGG
35236_rex                    CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
35855_rex                    CAAACGCAACTTTCAATATRGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGYTGAGGTTAATGAGG
38362_rex                    CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
39618_rex                    CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
40578_rex                    CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
41478_cyathophylloides       CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
41954_cyathophylloides       CAAACGCAACTTTCAATATGGGTTCTGACGAAGCCAATTTAGTAATTAGTAAAGCTGAGGTTAATGAGG
//                                   *      -   -     *                            -              |18:locus_20_1:1-69
assembly_reference_sequence  CAAACGCAGCTTACAATATGGGTTCCGACGAGGCTAATTTAGTAATTAGTAAAGCTGAGGTTAATGAAG
30686_cyathophylla           CAAACGCAGCTTACAATATGGGTTCCGACGAGGCTAATTTAGTAATTAGTAAAGCTGAGGTTAATGAAG
32082_przewalskii            CAAACGCAGCTTACAATATGGGTTCCRACGAGGCTAATTTAGTAATTAGTAAAGCTGARGTTAATGAAG
33588_przewalskii            CAAACGCAGCTTACAATATGGGTTCCAACGAGGCTAATTTAGTAATTAGTAAAGCTGAGGTTAATGAAG
35855_rex                    CAAACACAACTTACAATATGGGTTCCGACGAGGCCAATTCAGTAATTAGTAAAGTTGAGGTTAATGAGG
41478_cyathophylloides       CAAACGCAGCTTACAATATGGGTTCCGACGAGGCTAATTTAGTAATTAGTAAAGCTGAGGTTAATGAAG
41954_cyathophylloides       CAAACGCAGCTTACAATATGGGTTCCGACGAGGCTAATTTAGTAATTAGTAAAGCTGAGGTTAATGAAG
//                                -  -                 *       -    -              -   -        - |19:locus_20_3:1-69
assembly_reference_sequence  CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
29154_superba                CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
30556_thamno                 CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
30686_cyathophylla           CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
32082_przewalskii            CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
33413_thamno                 CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
33588_przewalskii            CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
35236_rex                    CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
35855_rex                    CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
38362_rex                    CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
39618_rex                    CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
40578_rex                    CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
41478_cyathophylloides       CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
41954_cyathophylloides       CTATCGGTTTGCACTTTTACCCAATCTGGGAAGCAGCATCCGTTGATGAATGGTTATACAATGG
//                                                                                           |20:locus_21_1:1-64
```


#### output files

```bash
ls -lh SRP021469/OUT/
```

```literal
total 37M
-rw-rw-r-- 1 deren deren 948K Jul 20 15:11 assembly.bed
-rw-rw-r-- 1 deren deren 5.7M Jul 20 15:11 assembly.loci.gz
-rw-rw-r-- 1 deren deren 6.1M Jul 20 15:12 assembly.vcf.gz
-rw-rw-r-- 1 deren deren 521K Jul 20 15:12 assembly.vcf.gz.csi
-rw-rw-r-- 1 deren deren  24M Jul 20 15:12 assembly.hdf5
-rw-rw-r-- 1 deren deren  11K Jul 20 15:12 assembly.stats.txt
-rw-rw-r-- 1 deren deren  15K Jul 20 15:12 assembly.stats.json
```

## Analysis

This dataset is used as an example for each tool in the Analysis section. Below I show just one example
of using the window extracter (wex) tool filter and write a concatenated alignment, followed by raxml-ng
to infer a ML phylogenetic tree.

### window-extracter (wex)

```bash
ipyrad2 wex \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/output-wex \
  -n assembly_min8 \
  -m 8 \
  -r 0.9
```

```literal
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | -------------------------------------------------------
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | ----- ipyrad2 wex: extract alignments from windows -----
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | -------------------------------------------------------
2026-07-20 16:59:25 | INFO     | cli_analysis.py      | CMD: ipyrad2 wex -d SRP021469/OUT/assembly.hdf5 -m 8 -o SRP021469/output-wex -n assembly_min8 -r 0.9
2026-07-20 16:59:25 | INFO     | window_extracter.py  | No windows specified; selecting the full length of all scaffolds. Use -w to subset scaffold windows and -P to view scaffold names.
2026-07-20 16:59:25 | INFO     | window_extracter.py  | selected 45448 windows from 45448 scaffolds
2026-07-20 16:59:43 | INFO     | window_extracter.py  | wrote alignment (13, 1965012) to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.phy
2026-07-20 16:59:43 | INFO     | window_extracter.py  | wrote stats/log to: /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.stats.txt
```

This generated an alignment that is 13 taxa x 1.96M sites. Examining the stats file we can see additional information.

```bash
cat SRP021469/output-wex/assembly_min8.stats.txt
```

```literal
CMD: ipyrad2 wex -d SRP021469/OUT/assembly.hdf5 -m 8 -o SRP021469/output-wex -n assembly_min8 -r 0.9

# Extract Summary
infile                    SRP021469/OUT/assembly.hdf5
outfile                   /home/deren/Documents/tools/ipyrad2/SRP021469/output-wex/assembly_min8.phy
out_format                phy
windows_selected          45,448
selected_windows_preview  locus_1_1:1-69, locus_2_1:1-69, locus_3_1:1-72, locus_4_25:1-71, locus_5_1:1-69, locus_6_1:1-69, locus_7_1:1-69, locus_8_1:1-69, locus_9_1:1-62, locus_10_1:1-69, ... (45448 total)

# Filtering Summary
populations                     all
min_sample_coverage_filter      all=8
max_sample_missing              0.900000
samples_selected_initial        13
samples_dropped_by_max_missing  0
samples_final                   13

# Alignment Summary
nsamples_before_filtering              13
nsites_in_windows_before_filtering     3,100,317
nvariants_in_windows_before_filtering  183,379
nsamples_after_filtering               13
nsites_in_windows_after_filtering      1,965,012
nvariants_in_windows_after_filtering   126,655

# Sample Summary
sample                  population  percent_missing  dropped_by_max_missing
29154_superba           all         26.241           no
30556_thamno            all         6.652            no
30686_cyathophylla      all         14.528           no
32082_przewalskii       all         49.630           no
33413_thamno            all         33.084           no
33588_przewalskii       all         43.547           no
35236_rex               all         4.721            no
35855_rex               all         4.861            no
38362_rex               all         4.950            no
39618_rex               all         13.964           no
40578_rex               all         3.844            no
41478_cyathophylloides  all         6.613            no
41954_cyathophylloides  all         10.684           no
```

Now that have written a concatenated phylip file with the supermatrix alignment we can run a tree inference tool on it, such as `raxml-ng` below.

### raxml-ng concatenation tree

Here we use `raxml-ng` to infer a phylogenetic tree from the supermatrix alignment generated by `wex`. We specify the `--all` option which will perform a tree search and bootstrap analysis to calculate support values. We specify the input (`--msa`) as the phy file produced in the previous step, indicate the substitution model choice, here using the common default `GTR+G`, and tell it to perform 100 non-parametric bootstrap replicate searches (`--bs-trees`). This will likely take 20 minutes or more to run.

```bash
raxml-ng \
    --all \
    --msa SRP021469/output-wex/assembly_min8.phy \
    --model GTR+G \
    --bs-trees 100 \
    --workers 2
```

When it finishes you an examine the results by plotting a tree. The default result files are saved to the same folder as the input file.
One easy way to plot a tree is using [`toytree`](https://eaton-lab.org/toytree/), which can plot trees either in the terminal, or as high quality vector graphics in formats like PDF. It can also be used to perform operations on a tree, such as re-rooting on an outgroup, as shown below.

Let's start by rooting the tree on an outgroup, which in this case is the taxa labeled 'przewalksii'.

```bash
toytree root \
    -i SRP021469/output-wex/assembly_min8.phy.raxml.support \
    -o SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted \
    -n "~prz" \
    --mad
```

Then we can print a tree visualization to the terminal:

```bash
toytree view \
    -i SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted \
    --ladderize
```

```literal
                               ┌───32082_przewalskii
┌──────────────────────────────┤
│                              └───33588_przewalskii
│
│                       ┌─────────29154_superba
│                    ┌──┤
│                    │  └─────────30686_cyathophylla
│                ┌───┤
│                │   │            ┌─41954_cyathophylloides
│                │   └────────────┤
└────────────────┤                └─41478_cyathophylloides
                 │
                 │      ┌───────33413_thamno
                 │      │
                 └──────┤┌─────────30556_thamno
                        ││
                        └┤   ┌────35855_rex
                         │┌──┤
                         ││  └────40578_rex
                         └┤
                          │┌────────35236_rex
                          └┤
                           │      ┌──38362_rex
                           └──────┤
                                  └──39618_rex
```

Or generate a high quality PDF tree visualization and open it externally:

```bash
toytree draw \
    -i SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted \
    -o SRP021469/output-wex/assembly_min8.phy.raxml.support.rooted.pdf \
    --node-labels 'support' \
    --ladderize
```

Or you can open a python or jupyter session and use toytree interactively in Python to generate
a tree figure with many more styling options. See the toytree docs.


## EXIT


### treeslider tree set

Another useful phylogenetic analysis is to infer a *species tree* using ASTRAL. Here, you must first infer a gene tree
for each locus individually, and then analyze the distribution of gene trees to find the best species tree that can explain
the variation among trees under assumptions of multi-species coalescent model. This takes into account the expectation that
incomplete lineage sorting will cause gene tree variation among close relatives.

The `treeslider` tool in ipyrad2 makes it easy to infer a gene tree for each locus, while also filtering the dataset to only
consider loci that meet some minimum filtering requirement. For example, you can require that all samples are present, or that
at least a certain subset of samples are present using the `-m`, `--imap`, and `--minmap` arguments.
Here we just use `-m 10` to require that a locus has data for at least 13 samples. After filtering the set of loci, it will run
a `raxml-ng` analysis for each one to generate a tree, with the final set of trees saved to the result table.

```bash
ipyrad2 treeslider \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/OUT/output-treeslider \
  -m 10
```

```literal
...
```

### astral species tree

Create an IMAP file mapping each sample to a group/population name.

```literal
32082_przewalskii       przewalskii
33588_przewalskii       przewalskii
29154_superba           superba
30686_cyathophylla      cyathophylla
41954_cyathophylloides  cyathophylloides
41478_cyathophylloides  cyathophylloides
33413_thamno            thamnophila_subsp_cupuliformis
30556_thamno            thamnophila_subsp_thamnophila
35855_rex               rex_subsp_rex
40578_rex               rex_subsp_rex
35236_rex               rex_subsp_rockii
38362_rex               rex_subsp_lipskyana
39618_rex               rex_subsp_lipskyana
```

Infer a species tree in astral4 by providing the input set of trees and mapping file.

```bash
astral4 \
    --input ... \
    --mapping ... \
    --output ... \
    --thread 6 \
    --root przewalskii \
    ...
```

```bash
toytree
```

### PCA

```bash
ipyrad2 pca \
  -d SRP021469/OUT/assembly.hdf5 \
  -o SRP021469/OUT/output-pca \
  ...
```


### BPP model fit

The `ipyrad2 bpp` tool makes it easy to setup and run BPP analyses from an ipyrad assembly. This includes
fitting MSC, MSC-i, and MSC-m models, inferring species trees, and performing species delimitation.

```bash
ipyrad2 bpp \
  -d ... \
  -o ... \

```
