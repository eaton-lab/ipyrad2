

# Empirical PE reference assembly tutorial

This tutorial demonstrates a reference-based assembly for an empirical 2RAD dataset in ipyrad2
with additional WGS samples included. All reads are paired-end 2 x 150bp in length. The RAD
dataset was prepared using the 3RAD protocol.
 <!-- (enzyme set 1: Cla1 and HindIII restriction enzymes). -->


## Download the data

The data in this tutorial are not yet publicly available.

<!--

mkdir -p TUTORIAL/REFERENCE/
mkdir -p TUTORIAL/FASTQS/{RAD,WGS}
(
  cd TUTORIAL/FASTQS/RAD/
  ln -s ../../../examples/Ama-PE-ddRAD/*.gz .
)
(
  cd TUTORIAL/FASTQS/WGS/
  ln -s ../../../examples/Ama-WGS/*.gz .
)
(
  cd TUTORIAL/REFERENCE/
  ln -s ../../examples/Atub-genome/* .
)

 -->

## Assembly

### trim RAD data

The data are already demultiplexed to individual samples, so we can start by running read trimming. Here we
will run `trim` twice, first on the RAD-seq samples, and then on the WGS samples. This is will allow the
`trim` method to detect and remove the restriction motif from the beginning of each RAD-seq read. These
patterns are not present on the WGS sequences, so we will run trimming on the RAD and WGS samples separately.

Here we use `-d` to specify input fastq paths, and `-o` to specify output directory paths, and also specify
the number of cores available (`-c`) and how to distribute these resources among multi-threaded jobs (`-t`).

<!-- (Cla1 motif) -->
<!-- (HindIII motif).  -->
The logging report shows that it detected 19 paired samples from our input data paths, and
that these samples consistely start with the restriction cutsite motifs "ATCGG" on R1s
and "TAGCTT" on R2s.
It reports the progress of the run until completed, and
then prints the path to the output fastqs and stats file.

```bash
ipyrad2 trim \
  -d TUTORIAL/FASTQS/RAD/*.fastq.gz \
  -o TUTORIAL/TRIM/RAD/ \
  -c 12 -t 4
```

??? note "ipyrad2 trim rad log"

    ```literal
    2026-07-22 15:28:04 | INFO     | cli_main.py          | ----------------------------------------------------------
    2026-07-22 15:28:04 | INFO     | cli_main.py          | ----- ipyrad2 trim: quality, adapter, and cutsite motif trimming -----
    2026-07-22 15:28:04 | INFO     | cli_main.py          | ----------------------------------------------------------
    2026-07-22 15:28:04 | INFO     | cli_main.py          | CMD: ipyrad2 trim -d TUTORIAL/FASTQS/RAD/SLH_AL_0012_R1.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0012_R2.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0013_R1.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0013_R2.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0014_R1.fastq.gz ...[truncated; 38 total matched paths] -o TUTORIAL/TRIM/RAD/ -c 8 -t 4 -f
    2026-07-22 15:28:04 | INFO     | names.py             | paired files by auto-detecting mate tokens in filenames
    2026-07-22 15:28:04 | INFO     | names.py             | showing first 10/19 names parsed from file paths
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0012 <- ('SLH_AL_0012_R1.fastq.gz', 'SLH_AL_0012_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0013 <- ('SLH_AL_0013_R1.fastq.gz', 'SLH_AL_0013_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0014 <- ('SLH_AL_0014_R1.fastq.gz', 'SLH_AL_0014_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0018 <- ('SLH_AL_0018_R1.fastq.gz', 'SLH_AL_0018_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0030 <- ('SLH_AL_0030_R1.fastq.gz', 'SLH_AL_0030_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0034 <- ('SLH_AL_0034_R1.fastq.gz', 'SLH_AL_0034_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0036 <- ('SLH_AL_0036_R1.fastq.gz', 'SLH_AL_0036_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0042 <- ('SLH_AL_0042_R1.fastq.gz', 'SLH_AL_0042_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0048 <- ('SLH_AL_0048_R1.fastq.gz', 'SLH_AL_0048_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | names.py             | SLH_AL_0063 <- ('SLH_AL_0063_R1.fastq.gz', 'SLH_AL_0063_R2.fastq.gz')
    2026-07-22 15:28:04 | INFO     | trim_fastqs.py       | trim input preflight found 19 usable samples and 0 skipped empty samples
    [####################] 100% | Counting kmers - total jobs: 19
    [####################] 100% | Counting kmers - total jobs: 19
    2026-07-22 15:28:08 | INFO     | trim_fastqs.py       | cutsite motifs set to R1=[ATCGG] at offset 0 R2=[TAGCTT] at offset 0
    2026-07-22 15:28:08 | INFO     | trim_fastqs.py       | trimming/filtering 19 samples with 'fastp' and writing to /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM
    2026-07-22 15:28:08 | INFO     | trim_fastqs.py       | running up to 2 parallel jobs each using up to 4 threads
    [####################] 100% | Trimming - total jobs: 19
    2026-07-22 15:38:26 | INFO     | trim_fastqs.py       | trimming stats written to /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM/RADipyrad_trim_stats_0.txt and /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM/RAD/ipyrad_trim_stats_0.json
    ```

The stats file report shows the following. There are ~500K read pairs per sample. The proportion of bases with quality scores >20 or >30 increases slightly after trimming. About ~20K reads were filtered from each sample, and ~100K bases were trimmed from each sample.

```bash
cat TUTORIAL/TRIM/RAD/ipyrad_trim_stats_0.txt
```

```literal
CMD: ipyrad2 trim -d TUTORIAL/FASTQS/RAD/SLH_AL_0012_R1.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0012_R2.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0013_R1.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0013_R2.fastq.gz TUTORIAL/FASTQS/RAD/SLH_AL_0014_R1.fastq.gz ...[truncated; 38 total matched paths] -o TUTORIAL/TRIM/RAD/ -c 8 -t 4 -f

            total_reads_before total_bases_before q20_rate_before q30_rate_before read1_mean_length_before read2_mean_length_before total_reads_after total_bases_after q20_rate_after q30_rate_after read1_mean_length_after read2_mean_length_after reads_filtered_by_low_quality reads_filtered_by_too_many_N reads_filtered_by_low_complexity reads_filtered_by_too_short adapter_trimmed_reads adapter_trimmed_bases
SLH_AL_0012             563176           80001289        0.954069        0.911111                      141                      142            497944          66873966       0.981349       0.948682                     134                     134                         63800                           32                                2                        1398                 16966                105734
SLH_AL_0013             146094           20281405        0.954693        0.910566                      139                      138            129446          16969956       0.981259       0.947499                     132                     129                         16214                           10                                2                         422                  8521                 61117
SLH_AL_0014             443004           61969524        0.956256        0.912269                      139                      140            395550          52340611       0.983580       0.950033                     132                     132                         46132                           52                                4                        1266                 14011                100130
SLH_AL_0018             399974           56367849        0.955489        0.912426                      141                      140            354654          47231153       0.981967       0.949155                     134                     131                         44304                           34                                4                         978                 12822                 82471
SLH_AL_0030             549442           77490308        0.949871        0.904072                      140                      141            479712          63927960       0.979807       0.945329                     133                     133                         68412                           62                               16                        1240                 16358                106977
SLH_AL_0034             919528          129257044        0.957052        0.913076                      139                      141            824964         109686676       0.983515       0.949865                     132                     133                         92498                           86                                6                        1974                 25575                172774
SLH_AL_0036             560770           79099258        0.956150        0.912363                      139                      142            501304          66891797       0.983463       0.950077                     132                     134                         58250                           66                               16                        1134                 16225                105836
SLH_AL_0042             541960           76126373        0.948882        0.902434                      140                      140            471154          62526059       0.979608       0.944598                     133                     132                         69448                           44                                6                        1308                 16587                115752
SLH_AL_0048             381132           53734169        0.956733        0.913493                      139                      142            341132          45510921       0.983839       0.950918                     132                     134                         39024                           44                                2                         930                 11480                 73544
SLH_AL_0063             322578           45163232        0.956649        0.915023                      141                      138            287666          38081673       0.981918       0.949984                     133                     130                         33886                           34                                0                         992                 13174                 96015
SLH_AL_0064             706128          100558899        0.952658        0.908034                      142                      142            619146          83410659       0.980672       0.947017                     135                     134                         85520                           62                                2                        1398                 21361                114197
SLH_AL_0084             249296           35235114        0.956800        0.914646                      140                      141            221634          29619760       0.982823       0.950826                     133                     133                         26972                           30                                2                         658                 11371                 66049
SLH_AL_0100             630198           88883272        0.951891        0.907221                      142                      139            553546          73808801       0.980239       0.946311                     135                     131                         75416                           74                                2                        1160                 17966                118649
SLH_AL_0101             894498          127625723        0.950925        0.905928                      142                      142            782744         105652472       0.980218       0.946168                     135                     134                        110264                           70                                2                        1418                 21758                111997
SLH_AL_0104             728720          103161918        0.952420        0.908515                      140                      142            641002          85810622       0.980896       0.947689                     133                     134                         86212                           78                                8                        1420                 21701                137332
SLH_AL_0105             529698           73862769        0.958234        0.915753                      139                      139            476570          62883408       0.983888       0.951379                     132                     131                         51740                           52                                8                        1328                 16232                126928
SLH_AL_0106             676214           96331896        0.953310        0.909506                      142                      142            594574          80161167       0.980783       0.947751                     135                     134                         80058                           64                                2                        1516                 20532                115082
SLH_AL_3065             200288           28146435        0.953608        0.910337                      139                      141            176412          23447184       0.981538       0.948949                     132                     133                         23316                           10                                4                         546                 10433                 64039
SLH_AL_3066             461010           65195993        0.956640        0.914975                      141                      141            411400          55046995       0.981874       0.949931                     134                     133                         48518                           36                                0                        1056                 14140                 93647
```

### trim WGS data

Next we run `trim` on the WGS samples. You can skip this step if you do not wish to add any
WGS samples. Because these samples contain many more reads they take a bit longer to run.
To improve runtimes, and normalize inputs among samples, you can optionally use the
`-x/--max-reads` flag here to keep only the first N number of reads from any sample.

```bash
ipyrad2 trim \
  -d TUTORIAL/FASTQS/WGS/*.fastq.gz \
  -o TUTORIAL/TRIM/WGS/ \
  -x 5_000_000 \
  -E \
  -c 12 -t 4
```

??? note "ipyrad2 trim wgs log"

    ```literal
    2026-07-22 18:38:59 | INFO     | cli_main.py          | ----------------------------------------------------------
    2026-07-22 18:38:59 | INFO     | cli_main.py          | ----- ipyrad2 trim: quality, adapter, and cutsite motif trimming -----
    2026-07-22 18:38:59 | INFO     | cli_main.py          | ----------------------------------------------------------
    2026-07-22 18:38:59 | INFO     | cli_main.py          | CMD: ipyrad2 trim -d TUTORIAL/FASTQS/WGS/21040XD-01-07_S39_L002_R1_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-07_S39_L002_R2_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-08_S40_L002_R1_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-08_S40_L002_R2_001.fastq.gz ...[truncated; 8 total matched paths] -o TUTORIAL/TRIM/WGS/ -x 5_000_000 -E -c 12 -t 4
    2026-07-22 18:38:59 | INFO     | names.py             | paired files by auto-detecting mate tokens in filenames
    2026-07-22 18:38:59 | INFO     | names.py             | showing first 4/4 names parsed from file paths
    2026-07-22 18:38:59 | INFO     | names.py             | 21040XD-01-07_S39_L002 <- ('21040XD-01-07_S39_L002_R1_001.fastq.gz', '21040XD-01-07_S39_L002_R2_001.fastq.gz')
    2026-07-22 18:38:59 | INFO     | names.py             | 21040XD-01-08_S40_L002 <- ('21040XD-01-08_S40_L002_R1_001.fastq.gz', '21040XD-01-08_S40_L002_R2_001.fastq.gz')
    2026-07-22 18:38:59 | INFO     | names.py             | 21040XD-01-09_S41_L002 <- ('21040XD-01-09_S41_L002_R1_001.fastq.gz', '21040XD-01-09_S41_L002_R2_001.fastq.gz')
    2026-07-22 18:38:59 | INFO     | names.py             | SRR15412865            <- ('SRR15412865_1.fastq.gz', 'SRR15412865_2.fastq.gz')
    2026-07-22 18:38:59 | INFO     | trim_fastqs.py       | trim input preflight found 4 usable samples and 0 skipped empty samples
    2026-07-22 18:38:59 | INFO     | trim_fastqs.py       | cutsite motifs set to R1=[<none>] at offset 0 R2=[<none>] at offset 0
    2026-07-22 18:38:59 | INFO     | trim_fastqs.py       | trimming/filtering 4 samples with 'fastp' and writing to /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM/WGS
    2026-07-22 18:38:59 | INFO     | trim_fastqs.py       | running up to 3 parallel jobs each using up to 4 threads
    [####################] 100% | Trimming - total jobs: 4
    2026-07-22 18:54:09 | INFO     | trim_fastqs.py       | trimming stats written to /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM/WGS/ipyrad_trim_stats_0.txt and /home/deren/Documents/ipyrad-tests/TUTORIAL/TRIM/WGS/ipyrad_trim_stats_0.json
    ```

The stats report shows a slightly higher adapter contamination in this dataset, with one sample trimming adapters in nearly 10% of reads. The base quality is lower in the one sample with 300 bp reads than in the other samples with 150 bp reads, reflecting differences in the technologies used to sequence these samples.

```bash
cat TUTORIAL/TRIM/WGS/ipyrad_trim_stats_0.txt
```

```literal
CMD: ipyrad2 trim -d TUTORIAL/FASTQS/WGS/21040XD-01-07_S39_L002_R1_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-07_S39_L002_R2_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-08_S40_L002_R1_001.fastq.gz TUTORIAL/FASTQS/WGS/21040XD-01-08_S40_L002_R2_001.fastq.gz ...[truncated; 8 total matched paths] -o TUTORIAL/TRIM/WGS/ -x 5_000_000 -E -c 12 -t 4

                       total_reads_before total_bases_before q20_rate_before q30_rate_before read1_mean_length_before read2_mean_length_before total_reads_after total_bases_after q20_rate_after q30_rate_after read1_mean_length_after read2_mean_length_after reads_filtered_by_low_quality reads_filtered_by_too_many_N reads_filtered_by_low_complexity reads_filtered_by_too_short adapter_trimmed_reads adapter_trimmed_bases
21040XD-01-07_S39_L002           10000000         1510000000        0.961961        0.917460                      151                      151           8718158        1298721695       0.987971       0.962181                     149                     148                       1248882                           26                              398                       32536                370126               7662366
21040XD-01-08_S40_L002           10000000         1510000000        0.969503        0.931227                      151                      151           9090450        1334383169       0.989339       0.965528                     147                     146                        842934                           10                              476                       66130               1116069              31804344
21040XD-01-09_S41_L002           10000000         1510000000        0.960757        0.916889                      151                      151           8805934        1304337361       0.987865       0.962097                     148                     147                       1067662                           22                              246                      126136                505151              21563173
SRR15412865                      10000000         2494480397        0.936768        0.896631                      249                      249           8258378        2020158740       0.981638       0.967588                     246                     242                       1721446                          818                              398                       18960                 15326               1531004
```

### map RAD data

Next map the trimmed RAD fastqs to the reference genome.

```bash
ipyrad2 map \
  -d TUTORIAL/TRIM/RAD/*.fastq.gz \
  -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa \
  -o TUTORIAL/MAP/RAD/ \
  -c 12 -t 4
```

??? note "ipyrad2 map rad log"

    ```literal
    2026-07-22 17:17:26 | INFO     | cli_main.py          | CMD: ipyrad2 map -d TUTORIAL/TRIM/RAD/SLH_AL_0012.R1.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0012.R2.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0013.R1.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0013.R2.trimmed.fastq.gz ...[truncated; 38 total matched paths] -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa -o TUTORIAL/MAP/RAD/ -c 8 -t 4 -f
    2026-07-22 17:17:26 | INFO     | names.py             | paired files by auto-detecting mate tokens in filenames
    2026-07-22 17:17:26 | INFO     | names.py             | showing first 10/19 names parsed from file paths
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0012 <- ('SLH_AL_0012.R1.trimmed.fastq.gz', 'SLH_AL_0012.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0013 <- ('SLH_AL_0013.R1.trimmed.fastq.gz', 'SLH_AL_0013.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0014 <- ('SLH_AL_0014.R1.trimmed.fastq.gz', 'SLH_AL_0014.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0018 <- ('SLH_AL_0018.R1.trimmed.fastq.gz', 'SLH_AL_0018.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0030 <- ('SLH_AL_0030.R1.trimmed.fastq.gz', 'SLH_AL_0030.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0034 <- ('SLH_AL_0034.R1.trimmed.fastq.gz', 'SLH_AL_0034.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0036 <- ('SLH_AL_0036.R1.trimmed.fastq.gz', 'SLH_AL_0036.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0042 <- ('SLH_AL_0042.R1.trimmed.fastq.gz', 'SLH_AL_0042.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0048 <- ('SLH_AL_0048.R1.trimmed.fastq.gz', 'SLH_AL_0048.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | names.py             | SLH_AL_0063 <- ('SLH_AL_0063.R1.trimmed.fastq.gz', 'SLH_AL_0063.R2.trimmed.fastq.gz')
    2026-07-22 17:17:26 | INFO     | mapper.py            | using existing bwa-mem2 reference index: AmaTu_v01_no00_renamed.fa
    2026-07-22 17:17:26 | INFO     | mapper.py            | mapping 19 samples to coordinate-sorted BAMs in /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/RAD
    2026-07-22 17:17:26 | INFO     | mapper.py            | using up to 8 cores (up to 2 multi-threaded jobs using 4 threads)
    [####################] 100% | Mapping - total jobs: 19
    [####################] 100% | Gathering mapping stats - total jobs: 19
    2026-07-22 17:33:08 | INFO     | mapper.py            | mapping stats written to /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/RAD/ipyrad_map_stats_0.txt and /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/RAD/ipyrad_map_stats_0.json
    ```

The stats file shows the following: Approximately 60% of read pairs were retained in each sample, after removing reads that were either unmapped or not primary alignments, or that did not map to the same scaffold. This file also includes a table below the main stats with a "preview", showing how typical filters that will be applied in the next stage (assemble) will affect these reads. The mapping scores (MAPQ) are quite high (mean 60), and the number of reads that are soft-clipped (reflecting that only part of the read mapped) is
generally low. This indicates accurate and unique mapping of reads.

```bash
cat TUTORIAL/MAP/RAD/ipyrad_map_stats_0.txt
```

```literal
CMD: ipyrad2 map -d TUTORIAL/TRIM/RAD/SLH_AL_0012.R1.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0012.R2.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0013.R1.trimmed.fastq.gz TUTORIAL/TRIM/RAD/SLH_AL_0013.R2.trimmed.fastq.gz ...[truncated; 38 total matched paths] -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa -o TUTORIAL/MAP/RAD/ -c 8 -t 4 -f

# ipyrad2 map stats
# Final BAMs are coordinate sorted and indexed.
# Paired-end final BAMs keep only mapped mates on the same scaffold.

## Applied mapping summary
# These counts describe filters already applied during ipyrad2 map.

            input_templates reads_removed_unmapped_or_nonprimary reads_removed_same_scaffold_pairing duplicate_records_removed templates_in_final_bam fraction_input_templates_retained_in_final_bam
sample
SLH_AL_0012          290153                               107494                               83133                         0                 194840                                          0.672
SLH_AL_0013           75545                                30878                               22685                         0                  48764                                          0.645
SLH_AL_0014          226714                                78384                               61470                         0                 156787                                          0.692
SLH_AL_0018          206146                                75981                               56979                         0                 139666                                          0.678
SLH_AL_0030          274253                                94098                               74132                         0                 190138                                          0.693
SLH_AL_0034          471707                               166676                              130923                         0                 322908                                          0.685
SLH_AL_0036          286987                                99038                               80581                         0                 197178                                          0.687
SLH_AL_0042          269676                                89208                               71256                         0                 189444                                          0.702
SLH_AL_0048          195509                                72579                               53924                         0                 132258                                          0.676
SLH_AL_0063          175750                                83593                               61027                         0                 103440                                          0.589
SLH_AL_0064          385654                               191051                              144909                         0                 217674                                          0.564
SLH_AL_0084          136286                                69609                               47522                         0                  77721                                          0.570
SLH_AL_0100          342401                               163709                              126075                         0                 197509                                          0.577
SLH_AL_0101          487776                               238841                              184651                         0                 276030                                          0.566
SLH_AL_0104          396203                               190541                              137683                         0                 232091                                          0.586
SLH_AL_0105          288765                               130067                               99683                         0                 173890                                          0.602
SLH_AL_0106          365910                               170894                              128175                         0                 216376                                          0.591
SLH_AL_3065          106854                                47885                               33100                         0                  66362                                          0.621
SLH_AL_3066          251002                               112340                               85804                         0                 151930                                          0.605

## Assemble read-filter preview (not applied during mapping)
# These preview thresholds were not applied during mapping.
# Use them to guide ipyrad2 assemble read filters: -qm/--min-map-q, -ms/--max-softclip, -me/--max-nm, -mt/--max-tlen.
# Preview mode: pair-level thresholds evaluated on final BAM templates.
# MAPQ threshold: 20
# Soft-clipped bases threshold: 25
# NM threshold: 50
# Absolute TLEN threshold: 2000

### Preview filter effects
            templates_failing_min_mapq_20 templates_failing_max_softclip_25 templates_failing_max_nm_50 templates_failing_max_abs_tlen_2000 templates_passing_all_preview_filters fraction_templates_passing_all_preview_filters
sample
SLH_AL_0012                         18856                             35798                          24                                5254                                146946                                          0.754
SLH_AL_0013                          3552                             11407                           1                                1228                                 35581                                          0.730
SLH_AL_0014                         13796                             27875                          18                                3773                                120088                                          0.766
SLH_AL_0018                         12597                             26094                          15                                3459                                105733                                          0.757
SLH_AL_0030                         17494                             33001                          23                                4410                                145537                                          0.765
SLH_AL_0034                         26239                             56760                          57                                6696                                249280                                          0.772
SLH_AL_0036                         21345                             34426                          31                                4615                                147784                                          0.749
SLH_AL_0042                         18850                             32395                          38                                4392                                144192                                          0.761
SLH_AL_0048                          9749                             24783                          16                                3057                                101933                                          0.771
SLH_AL_0063                          9155                             31789                          13                                1953                                 68071                                          0.658
SLH_AL_0064                         19200                             69365                          72                                4459                                141478                                          0.650
SLH_AL_0084                          7132                             25495                          32                                1603                                 50155                                          0.645
SLH_AL_0100                         18499                             60753                          71                                4125                                129355                                          0.655
SLH_AL_0101                         27357                             86969                         113                                6173                                178390                                          0.646
SLH_AL_0104                         24521                             74828                          48                                7214                                148878                                          0.641
SLH_AL_0105                         16023                             51910                          46                                3462                                115442                                          0.664
SLH_AL_0106                         19855                             66010                          79                                4528                                142986                                          0.661
SLH_AL_3065                          6075                             20527                          40                                1614                                 43804                                          0.660
SLH_AL_3066                         13328                             45124                          29                                2870                                101385                                          0.667

### Preview metric summaries
            min_mapq_mean min_mapq_median min_mapq_stdev max_softclip_mean max_softclip_median max_softclip_stdev max_nm_mean max_nm_median max_nm_stdev abs_tlen_mean abs_tlen_median abs_tlen_stdev
sample
SLH_AL_0012        50.970          60.000         17.228            14.637               0.000             29.928       4.848         3.000        5.213    815719.101         262.000    6146044.205
SLH_AL_0013        51.509          60.000         15.184            19.401               0.000             34.062       4.088         2.000        4.780    985999.414         224.000    6883261.900
SLH_AL_0014        51.535          60.000         16.697            14.028               0.000             28.978       4.767         3.000        5.145    690880.530         254.000    5657727.339
SLH_AL_0018        51.456          60.000         16.695            14.884               0.000             29.886       4.777         3.000        5.199    824328.422         253.000    6229009.670
SLH_AL_0030        51.486          60.000         16.991            13.604               0.000             28.629       4.827         3.000        5.233    647552.637         259.000    5442898.687
SLH_AL_0034        52.095          60.000         16.261            13.824               0.000             28.882       4.742         3.000        5.157    596535.918         255.000    5248795.868
SLH_AL_0036        50.212          60.000         18.138            13.714               0.000             28.759       4.810         3.000        5.178    647023.665         259.000    5409961.391
SLH_AL_0042        50.963          60.000         17.533            13.275               0.000             27.991       4.900         3.000        5.207    633754.676         257.000    5401961.524
SLH_AL_0048        52.728          60.000         15.483            15.092               0.000             30.403       4.782         3.000        5.201    726915.591         255.000    5844245.247
SLH_AL_0063        51.002          60.000         16.744            22.992               2.000             33.322       8.009         7.000        6.073    164400.919         238.000    1984768.225
SLH_AL_0064        50.521          60.000         16.744            24.195               3.000             34.304       8.006         7.000        6.234    203219.342         248.000    2246946.658
SLH_AL_0084        50.532          60.000         16.940            25.440               4.000             35.474       7.920         7.000        6.347    199785.962         232.000    2189641.413
SLH_AL_0100        50.335          60.000         17.136            22.938               2.000             33.201       8.187         7.000        6.140    200412.323         250.000    2223153.129
SLH_AL_0101        50.098          60.000         17.467            23.637               3.000             33.793       8.439         8.000        6.228    214730.168         257.000    2282749.335
SLH_AL_0104        50.252          60.000         17.688            24.738               4.000             34.836       8.352         8.000        6.171    755339.140         248.000    5696474.467
SLH_AL_0105        50.918          60.000         17.062            22.015               2.000             32.360       8.266         7.000        6.087    185289.978         247.000    2116470.481
SLH_AL_0106        50.898          60.000         16.997            22.937               2.000             33.400       8.534         8.000        6.251    193010.700         257.000    2161833.714
SLH_AL_3065        50.888          60.000         16.822            24.089               3.000             34.795       7.922         7.000        6.360    614217.356         228.000    5205936.882
SLH_AL_3066        51.208          60.000         16.740            22.099               2.000             32.683       8.322         8.000        6.083    182257.623         252.000    2137144.878
```

### map WGS data

Next map the trimmed WGS fastqs to the reference genome. Here `-m` will remove PCR duplicates based on mapping coordinates,
a feature we can apply for WGS reads but not for RAD reads.

```bash
ipyrad2 map \
  -d TUTORIAL/TRIM/WGS/*.fastq.gz \
  -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa \
  -o TUTORIAL/MAP/WGS/ \
  -m \
  -c 8 -t 4
```

??? "ipyrad2 map wgs log"

    ```parsed-literal
    2026-07-22 18:58:06 | INFO     | cli_main.py          | --------------------------------------------------------------
    2026-07-22 18:58:06 | INFO     | cli_main.py          | ----- ipyrad2 map: map reads and write coordinate-sorted BAMs -----
    2026-07-22 18:58:06 | INFO     | cli_main.py          | --------------------------------------------------------------
    2026-07-22 18:58:06 | INFO     | cli_main.py          | CMD: ipyrad2 map -d TUTORIAL/TRIM/WGS/21040XD-01-07_S39_L002.R1.trimmed.fastq.gz TUTORIAL/TRIM/WGS/21040XD-01-07_S39_L002.R2.trimmed.fastq.gz TUTORIAL/TRIM/WGS/21040XD-01-08_S40_L002.R1.trimmed.fastq.gz ...[truncated; 8 total matched paths] -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa -o TUTORIAL/MAP/WGS/ -m -c 8 -t 4
    2026-07-22 18:58:06 | INFO     | names.py             | paired files by auto-detecting mate tokens in filenames
    2026-07-22 18:58:06 | INFO     | names.py             | showing first 4/4 names parsed from file paths
    2026-07-22 18:58:06 | INFO     | names.py             | 21040XD-01-07_S39_L002 <- ('21040XD-01-07_S39_L002.R1.trimmed.fastq.gz', '21040XD-01-07_S39_L002.R2.trimmed.fastq.gz')
    2026-07-22 18:58:06 | INFO     | names.py             | 21040XD-01-08_S40_L002 <- ('21040XD-01-08_S40_L002.R1.trimmed.fastq.gz', '21040XD-01-08_S40_L002.R2.trimmed.fastq.gz')
    2026-07-22 18:58:06 | INFO     | names.py             | 21040XD-01-09_S41_L002 <- ('21040XD-01-09_S41_L002.R1.trimmed.fastq.gz', '21040XD-01-09_S41_L002.R2.trimmed.fastq.gz')
    2026-07-22 18:58:06 | INFO     | names.py             | SRR15412865            <- ('SRR15412865.R1.trimmed.fastq.gz', 'SRR15412865.R2.trimmed.fastq.gz')
    2026-07-22 18:58:06 | WARNING  | mapper.py            | removing PCR duplicates by coordinates. Be sure this run includes only WGS samples, not RAD
    2026-07-22 18:58:06 | INFO     | mapper.py            | using existing bwa-mem2 reference index: AmaTu_v01_no00_renamed.fa
    2026-07-22 18:58:06 | INFO     | mapper.py            | mapping 4 samples to coordinate-sorted BAMs in /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/WGS
    2026-07-22 18:58:06 | INFO     | mapper.py            | using up to 8 cores (up to 2 multi-threaded jobs using 4 threads)
    [####################] 100% | Mapping - total jobs: 4
    [####################] 100% | Gathering mapping stats - total jobs: 4
    2026-07-23 09:26:10 | INFO     | mapper.py            | mapping stats written to /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/WGS/ipyrad_map_stats_0.txt and /home/deren/Documents/ipyrad-tests/TUTORIAL/MAP/WGS/ipyrad_map_stats_0.json
    ```

The stats report shows that the mapping rate for these samples was approximately ... Many read pairs were excluded becaues they
did not map, were not primary alignments, or paired to different scaffolds. Note, this could reflect these samples
being more distantly related to the reference genome, rather than being an artifact of being WGS versus RAD.


```bash
cat TUTORIAL/MAP/WGS/ipyrad_map_stats_0.txt
```
```literal
CMD: ipyrad2 map -d TUTORIAL/TRIM/WGS/21040XD-01-07_S39_L002.R1.trimmed.fastq.gz TUTORIAL/TRIM/WGS/21040XD-01-07_S39_L002.R2.trimmed.fastq.gz TUTORIAL/TRIM/WGS/21040XD-01-08_S40_L002.R1.trimmed.fastq.gz ...[truncated; 8 total matched paths] -r TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa -o TUTORIAL/MAP/WGS/ -m -c 8 -t 4

# ipyrad2 map stats
# Final BAMs are coordinate sorted and indexed.
# Paired-end final BAMs keep only mapped mates on the same scaffold.

## Applied mapping summary
# These counts describe filters already applied during ipyrad2 map.

                       input_templates reads_removed_unmapped_or_nonprimary reads_removed_same_scaffold_pairing duplicate_records_removed templates_in_final_bam fraction_input_templates_retained_in_final_bam
sample
21040XD-01-07_S39_L002         6221847                              3966472                             2741413                    347115                2694348                                          0.433
21040XD-01-08_S40_L002         6120214                              3373345                             2131772                    496787                3119263                                          0.510
21040XD-01-09_S41_L002         6188210                              3921824                             2829802                    762609                2431094                                          0.393
SRR15412865                    7548622                              8714761                             4892654                      3938                 742946                                          0.098

## Assemble read-filter preview (not applied during mapping)
# These preview thresholds were not applied during mapping.
# Use them to guide ipyrad2 assemble read filters: -qm/--min-map-q, -ms/--max-softclip, -me/--max-nm, -mt/--max-tlen.
# Preview mode: pair-level thresholds evaluated on final BAM templates.
# MAPQ threshold: 20
# Soft-clipped bases threshold: 25
# NM threshold: 50
# Absolute TLEN threshold: 2000

### Preview filter effects
                       templates_failing_min_mapq_20 templates_failing_max_softclip_25 templates_failing_max_nm_50 templates_failing_max_abs_tlen_2000 templates_passing_all_preview_filters fraction_templates_passing_all_preview_filters
sample
21040XD-01-07_S39_L002                        931985                           1089614                        1605                              150000                               1146875                                          0.426
21040XD-01-08_S40_L002                       1124354                           1111925                        1492                              135900                               1392662                                          0.446
21040XD-01-09_S41_L002                        846003                           1105404                        1482                              158592                               1005254                                          0.413
SRR15412865                                   414198                            618744                        1217                              141652                                 91999                                          0.124

### Preview metric summaries
                       min_mapq_mean min_mapq_median min_mapq_stdev max_softclip_mean max_softclip_median max_softclip_stdev max_nm_mean max_nm_median max_nm_stdev abs_tlen_mean abs_tlen_median abs_tlen_stdev
sample
21040XD-01-07_S39_L002        34.851          40.000         25.218            33.562              11.000             41.737       9.957         9.000        6.962    454316.974         241.000    3231238.069
21040XD-01-08_S40_L002        33.990          40.000         25.243            29.223               7.000             39.495       9.147         8.000        6.711    297615.529         199.000    2573490.135
21040XD-01-09_S41_L002        34.873          40.000         24.900            37.417              19.000             42.354      10.700        10.000        6.916    415525.864         237.000    2880333.931
SRR15412865                   23.402          12.000         24.111           161.495             209.000             85.971       7.367         2.000       10.350   2505987.176         346.000    7198230.542
```


### assemble

We are now ready to assemble the dataset. Here we specify separate paths to the RAD and WGS BAM alignments
using `--rad-bams` and `--wgs-bams`, respectively, and specify the path to the reference genome fasta (`-r`).
We saw in the mapping stats above that the mean MAPQ score was generally >50. Based on this we can confidently
exclude reads with much lower mapping scores, so I raised the `-qm` parameter to 40.
Most other options are left at their defaults.

Let's also rename our samples at this point so that the final assembled data contains easily interpretable
names, rather than obscure accession IDs, for downstream analyses. To do this we provide an 2-column file
using the `--rename` flag to map current names to new names.

```bash
ipyrad2 assemble \
  --rad-bams TUTORIAL/MAP/RAD/*.bam \
  --wgs-bams TUTORIAL/MAP/WGS/*.bam \
  --reference TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa \
  --out TUTORIAL/OUT/ \
  --rename TUTORIAL/NAMES.txt \
  --name assembly \
  -qm 40 \
  -c 8 -t 4
```
??? note "ipyrad2 assemble log"

    ```parsed-literal
    2026-07-23 09:32:26 | INFO     | cli_main.py          | -----------------------------------------------------------
    2026-07-23 09:32:26 | INFO     | cli_main.py          | ----- ipyrad2 assemble: delimit loci and call variants -----
    2026-07-23 09:32:26 | INFO     | cli_main.py          | -----------------------------------------------------------
    2026-07-23 09:32:26 | INFO     | cli_main.py          | CMD: ipyrad2 assemble --rad-bams TUTORIAL/MAP/RAD/SLH_AL_0012.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0013.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0014.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0018.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0030.trimmed.sorted.bam ...[truncated; 19 total matched paths] --wgs-bams TUTORIAL/MAP/WGS/21040XD-01-07_S39_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-08_S40_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-09_S41_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/SRR15412865.trimmed.sorted.bam --reference TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa --out TUTORIAL/OUT/ --name assembly -qm 40 -c 8 -t 4
    2026-07-23 09:32:26 | ERROR    | cli_main.py          | Error: outfiles with prefix assembly already exist in /home/deren/Documents/ipyrad-tests/TUTORIAL/OUT. Use --force to overwrite.
    2026-07-23 09:32:26 | ERROR    | cli_main.py          | see error message above. Shutting down.
    (ipyrad2) deren@rex ~/Documents/ipyrad-tests $ ipyrad2 assemble   --rad-bams TUTORIAL/MAP/RAD/*.bam   --wgs-bams TUTORIAL/MAP/WGS/*.bam   --reference TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa   --out TUTORIAL/OUT/   --name assembly   -qm 40   -c 8 -t 4 -f
    2026-07-23 09:32:29 | INFO     | cli_main.py          | -----------------------------------------------------------
    2026-07-23 09:32:29 | INFO     | cli_main.py          | ----- ipyrad2 assemble: delimit loci and call variants -----
    2026-07-23 09:32:29 | INFO     | cli_main.py          | -----------------------------------------------------------
    2026-07-23 09:32:29 | INFO     | cli_main.py          | CMD: ipyrad2 assemble --rad-bams TUTORIAL/MAP/RAD/SLH_AL_0012.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0013.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0014.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0018.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0030.trimmed.sorted.bam ...[truncated; 19 total matched paths] --wgs-bams TUTORIAL/MAP/WGS/21040XD-01-07_S39_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-08_S40_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-09_S41_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/SRR15412865.trimmed.sorted.bam --reference TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa --out TUTORIAL/OUT/ --name assembly -qm 40 -c 8 -t 4 -f
    2026-07-23 09:32:29 | INFO     | assemble.py          | loading BAM inputs
    2026-07-23 09:32:29 | INFO     | assemble.py          | loaded 19 RAD samples
    2026-07-23 09:32:29 | INFO     | assemble.py          | loaded 4 WGS samples
    [####################] 100% | Scanning BAM headers - total jobs: 23
    2026-07-23 09:32:31 | INFO     | assemble.py          | BAM layout: 23 paired-end, 0 single-end
    2026-07-23 09:32:34 | INFO     | assemble.py          | validating BAM headers against the reference
    2026-07-23 09:32:34 | INFO     | assemble.py          | filtering mapped reads before assembly: MAPQ>=40, same scaffold pairs only, abs(TLEN)<=1000
    [####################] 100% | Filtering mapped reads - total jobs: 23
    2026-07-23 09:33:41 | INFO     | assemble.py          | filtered analysis BAMs ready for 23 samples
    2026-07-23 09:33:41 | INFO     | assemble.py          | using up to 8 cores (2 concurrent jobs, 4 threads per job)
    2026-07-23 09:33:41 | INFO     | assemble.py          | building per-sample coverage BEDs
    [####################] 100% | Building per-sample coverage BEDs - total jobs: 23
    2026-07-23 09:36:11 | INFO     | assemble.py          | building loci from shared sample coverage BEDs
    2026-07-23 09:36:13 | INFO     | assemble.py          | scoring paralog evidence
    2026-07-23 09:36:19 | INFO     | assemble.py          | mixed RAD/WGS assembly detected; skipping softclip-based paralog failure for WGS samples
    2026-07-23 09:36:19 | INFO     | assemble.py          | paralog scoring uses the shared loci BED for all samples
    2026-07-23 09:36:19 | INFO     | assemble.py          | preparing loci-restricted BAMs for paralog scoring
    [####################] 100% | Preparing loci-restricted paralog BAMs - total jobs: 23
    2026-07-23 09:36:34 | INFO     | assemble.py          | loci-restricted paralog BAMs ready for 23 samples
    [####################] 100% | Scoring paralog evidence - total jobs: 23
    2026-07-23 09:40:34 | INFO     | assemble.py          | aggregating paralog filters across samples
    2026-07-23 09:40:34 | INFO     | assemble.py          | mixed RAD/WGS assembly detected; using RAD samples to control shared paralog locus retention and WGS samples for QC only
    2026-07-23 09:40:37 | INFO     | assemble.py          | paralog filtering retained 10693/12008 shared loci
    2026-07-23 09:40:37 | INFO     | assemble.py          | preparing cleaned BAMs for joint calling
    [####################] 100% | Preparing cleaned calling BAMs - total jobs: 23
    2026-07-23 09:40:50 | INFO     | assemble.py          | cleaned calling BAMs ready for 23 samples
    [####################] 100% | Calling variants - total jobs: 8
    2026-07-23 09:42:52 | INFO     | assemble.py          | filtering variant calls
    2026-07-23 09:42:56 | INFO     | assemble.py          | masking WGS heterozygous genotypes by allele balance
    2026-07-23 09:43:00 | INFO     | assemble.py          | masked 708 / 3008 WGS heterozygous genotypes outside allele-balance range [0.20, 0.80]
    2026-07-23 09:43:00 | INFO     | assemble.py          | resolving indels and SNPs
    2026-07-23 09:43:15 | INFO     | variants.py          | masked 598 overlapping-indel clusters (1317 records, 6249 bp)
    [####################] 100% | Building low-depth masks - total jobs: 23
    [####################] 100% | Building sample-specific paralog masks - total jobs: 23
    [####################] 100% | Merging sample masks - total jobs: 23
    2026-07-23 09:43:20 | INFO     | assemble.py          | preparing locus reference sequence
    2026-07-23 09:43:20 | INFO     | assemble.py          | building consensus sequences
    [####################] 100% | Building consensus sequences - total jobs: 23
    2026-07-23 09:43:30 | INFO     | assemble.py          | building locus database
    2026-07-23 09:43:30 | INFO     | assemble.py          | built locus database from 24 FASTA inputs
    2026-07-23 09:43:30 | INFO     | assemble.py          | writing final loci and summary files
    [####################] 100% | Resolving and writing final loci - total jobs: 84
    2026-07-23 09:44:02 | INFO     | assemble.py          | wrote final loci: 8699 loci, 2797293 sites
    2026-07-23 09:44:02 | INFO     | assemble.py          | writing final VCF
    [####################] 100% | Building final VCF masks - total jobs: 23
    2026-07-23 09:44:05 | INFO     | variants.py          | masking final VCF chunks
    [####################] 100% | Masking final VCF chunks - total jobs: 32
    2026-07-23 09:44:13 | INFO     | variants.py          | concatenating masked final VCF chunks
    2026-07-23 09:44:14 | INFO     | assemble.py          | wrote final VCF
    2026-07-23 09:44:17 | INFO     | assemble.py          | writing SNP database
    [####################] 100% | Building SNP database chunks - total jobs: 16
    2026-07-23 09:44:30 | INFO     | assemble.py          | wrote SNP database with 292466 SNP sites
    2026-07-23 09:44:30 | INFO     | assemble.py          | preparing final sample depth summaries
    [####################] 100% | Preparing final depth summaries - total jobs: 23
    2026-07-23 09:44:31 | INFO     | assemble.py          | summarizing final sample depth
    [####################] 100% | Summarizing final sample depth - total jobs: 23
    2026-07-23 09:44:33 | INFO     | assemble.py          | final sample depth summary ready for 23 samples
    2026-07-23 09:44:33 | INFO     | loci.py              | wrote assemble summary report
    2026-07-23 09:44:33 | INFO     | assemble.py          | removed assemble tmpdir /home/deren/Documents/ipyrad-tests/TUTORIAL/OUT/assembly_tmpdir
    2026-07-23 09:44:33 | INFO     | assemble.py          | assemble complete; outputs written to /home/deren/Documents/ipyrad-tests/TUTORIAL/OUT
    ```

The stats file report shows:


## Assembly stats

We now have a finished assembly stored in `TUTORIAL/OUT/`. The first thing to do is to look at the
human-readable stats file.

```bash
cat TUTORIAL/OUT/assembly.stats.txt
```
```literal
    CMD: ipyrad2 assemble --rad-bams TUTORIAL/MAP/RAD/SLH_AL_0012.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0013.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0014.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0018.trimmed.sorted.bam TUTORIAL/MAP/RAD/SLH_AL_0030.trimmed.sorted.bam ...[truncated; 19 total matched paths] --wgs-bams TUTORIAL/MAP/WGS/21040XD-01-07_S39_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-08_S40_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/21040XD-01-09_S41_L002.trimmed.sorted.bam TUTORIAL/MAP/WGS/SRR15412865.trimmed.sorted.bam --reference TUTORIAL/REFERENCE/AmaTu_v01_no00_renamed.fa --out TUTORIAL/OUT/ --name assembly -qm 40 -c 8 -t 4 -f

    # Assemble Summary
    Samples                                               23
    Shared loci before minimum sample coverage filter     26,905
    Shared loci after delimiting                          12,008
    Shared loci after paralog filtering                   10,693
    Final loci written                                    8,699
    Final loci retained fraction after paralog filtering  0.813523
    Final loci retained fraction after delimiting         0.724434
    Assembled sites                                       2,797,293
    Final SNP sites written                               292,466
    Variable sites                                        156,941
    Phylogenetically informative sites                    64,231
    Alignment matrix occupancy fraction                   0.360839
    Overlapping indel clusters masked                     598
    Overlapping indel records removed                     1,317
    Overlapping indel bases masked                        6,249

    # Mixed RAD/WGS Diagnostics
    RAD samples                                          19
    WGS samples                                          4
    Loci failed by RAD paralog QC                        1,315
    Loci failed by WGS paralog QC                        704
    Loci failed by RAD and WGS paralog QC                157
    Loci kept by RAD but failed by WGS QC                547
    WGS heterozygous genotypes masked by allele balance  708
    Sites supported by RAD only                          140,754
    Sites supported by WGS only                          7,477
    Sites supported by RAD and WGS                       8,608
    Sites supported by neither RAD nor WGS               135,627

    # Locus Filtering
    Loci filtered by minimum length                 7
    Loci filtered by minimum sample coverage        14,897
    Loci filtered by maximum variant frequency      0
    Loci filtered by maximum shared heterozygosity  1,992
    Loci filtered by maximum depth outlier          0

    # Sample Masking
    Loci with samples masked by minimum observed fraction threshold  1,384
    Sample masks triggered by minimum observed fraction threshold    1,547
    Loci with samples masked by sample heterozygosity threshold      5
    Sample masks triggered by sample heterozygosity threshold        5

    # Alignment Summary
    Mean locus length                           321.565
    Median locus length                         289.000
    Minimum locus length                        27
    Maximum locus length                        1,631
    Mean samples per locus                      9.522
    Median samples per locus                    8.000
    Sites with sample coverage >= 2             2,715,376
    Sites with sample coverage >= 3             2,706,445
    Sites with sample coverage >= 4             2,695,768
    Sites with sample coverage >= trim minimum  2,695,768

    # Sample Summary
    Sample                  Sample type  Read layout  Reads before filtering  Reads after filtering  Loci in alignment  Loci fraction in alignment  Shared loci with nonzero depth  Shared-depth loci fraction  Mean depth in shared loci  Median depth in shared loci  Mean depth in nonzero shared loci  Median depth in nonzero shared loci  Masked by minimum observed fraction threshold  Masked by sample heterozygosity threshold
    21040XD-01-07_S39_L002  WGS          PE           5,388,695               3,073,781              2,340              0.268996                    2,852                           0.327854                    0.874                      0.000                        2.666                              1.737                                509                                            2
    21040XD-01-08_S40_L002  WGS          PE           6,238,525               3,453,715              1,708              0.196344                    2,242                           0.257731                    0.606                      0.000                        2.352                              1.336                                533                                            1
    21040XD-01-09_S41_L002  WGS          PE           4,862,185               2,765,049              1,435              0.164961                    1,836                           0.211059                    0.551                      0.000                        2.609                              1.450                                401                                            0
    SLH_AL_0012             RAD          PE           389,680                 339,937                4,094              0.470629                    4,096                           0.470859                    7.787                      0.000                        16.538                             9.000                                5                                              0
    SLH_AL_0013             RAD          PE           97,528                  88,684                 388                0.044603                    389                             0.044718                    0.801                      0.000                        17.915                             5.000                                1                                              0
    SLH_AL_0014             RAD          PE           313,574                 276,308                3,752              0.431314                    3,758                           0.432004                    6.209                      0.000                        14.373                             8.000                                9                                              0
    SLH_AL_0018             RAD          PE           279,332                 246,319                3,358              0.386021                    3,358                           0.386021                    5.465                      0.000                        14.157                             7.291                                0                                              0
    SLH_AL_0030             RAD          PE           380,276                 333,838                4,126              0.474307                    4,134                           0.475227                    7.547                      0.000                        15.882                             9.374                                6                                              1
    SLH_AL_0034             RAD          PE           645,816                 575,161                4,824              0.554546                    4,826                           0.554776                    13.640                     6.000                        24.587                             14.396                               3                                              1
    SLH_AL_0036             RAD          PE           394,356                 337,531                4,106              0.472008                    4,107                           0.472123                    7.279                      0.000                        15.418                             9.000                                2                                              0
    SLH_AL_0042             RAD          PE           378,888                 328,846                4,109              0.472353                    4,111                           0.472583                    7.184                      0.000                        15.202                             9.000                                4                                              0
    SLH_AL_0048             RAD          PE           264,516                 238,431                3,306              0.380044                    3,311                           0.380618                    5.517                      0.000                        14.494                             8.000                                6                                              0
    SLH_AL_0063             RAD          PE           206,880                 182,070                3,750              0.431084                    3,751                           0.431199                    5.248                      0.000                        12.170                             8.000                                2                                              0
    SLH_AL_0064             RAD          PE           435,348                 382,907                5,934              0.682147                    5,935                           0.682262                    13.198                     9.000                        19.345                             13.313                               1                                              0
    SLH_AL_0084             RAD          PE           155,442                 136,433                1,950              0.224164                    1,957                           0.224968                    2.942                      0.000                        13.077                             6.928                                7                                              0
    SLH_AL_0100             RAD          PE           395,018                 343,862                5,894              0.677549                    5,895                           0.677664                    11.896                     9.000                        17.554                             13.000                               1                                              0
    SLH_AL_0101             RAD          PE           552,060                 477,004                6,174              0.709737                    6,174                           0.709737                    16.592                     13.000                       23.378                             18.000                               2                                              0
    SLH_AL_0104             RAD          PE           464,182                 400,608                3,973              0.456719                    3,976                           0.457064                    9.799                      0.000                        21.439                             13.000                               5                                              0
    SLH_AL_0105             RAD          PE           347,780                 304,528                4,986              0.573169                    4,988                           0.573399                    9.468                      6.000                        16.513                             12.000                               3                                              0
    SLH_AL_0106             RAD          PE           432,752                 378,310                6,009              0.690769                    6,009                           0.690769                    13.149                     10.410                       19.036                             15.000                               1                                              0
    SLH_AL_3065             RAD          PE           132,724                 116,456                1,805              0.207495                    1,811                           0.208185                    2.473                      0.000                        11.878                             6.000                                7                                              0
    SLH_AL_3066             RAD          PE           303,860                 267,303                4,669              0.536728                    4,671                           0.536958                    8.188                      5.000                        15.249                             11.000                               2                                              0
    SRR15412865             WGS          PE           1,485,892               580,822                141                0.016209                    178                             0.020462                    0.078                      0.000                        3.798                              1.535                                37                                             0

    # Locus Occupancy
    Samples with data  RAD loci before min sample coverage  RAD loci after min sample coverage  Final filtered RAD loci with WGS  Cumulative final loci  Fraction of final loci
    0                  0                                    0                                   0                                 0                      0.000000
    1                  9,138                                0                                   0                                 0                      0.000000
    2                  3,313                                0                                   0                                 0                      0.000000
    3                  2,238                                0                                   0                                 0                      0.000000
    4                  1,829                                1,798                               940                               940                    0.108058
    5                  1,380                                1,370                               981                               1,921                  0.112772
    6                  1,220                                1,203                               784                               2,705                  0.090125
    7                  1,407                                1,406                               1,009                             3,714                  0.115990
    8                  1,346                                1,349                               883                               4,597                  0.101506
    9                  778                                  749                                 627                               5,224                  0.072077
    10                 587                                  545                                 548                               5,772                  0.062996
    11                 417                                  389                                 394                               6,166                  0.045293
    12                 356                                  340                                 299                               6,465                  0.034372
    13                 376                                  366                                 259                               6,724                  0.029774
    14                 448                                  448                                 286                               7,010                  0.032877
    15                 461                                  470                                 343                               7,353                  0.039430
    16                 473                                  460                                 330                               7,683                  0.037935
    17                 500                                  493                                 356                               8,039                  0.040924
    18                 488                                  475                                 337                               8,376                  0.038740
    19                 150                                  147                                 183                               8,559                  0.021037
    20                 0                                    0                                   82                                8,641                  0.009426
    21                 0                                    0                                   24                                8,665                  0.002759
    22                 0                                    0                                   27                                8,692                  0.003104
    23                 0                                    0                                   7                                 8,699                  0.000805
```

We should be pretty satisifed with this result. We recovered many loci and without too much missing data --
i.e., most loci have data from nearly all samples, not just a subset of them. These loci seem to be highly
variable, with nearly X variant sites recovered.

## Assembled loci

After examining the stats file the next step should be to look at the loci file, which provides
a human-readable format for examining aligned loci to ensure that they look reasonable. This is
a large file so bash command-line tools like `less` or `head` are often useful for viewing it.

[describe loci format]

```bash
zcat TUTORIAL/OUT/assembly.loci.gz | head -n 202
```
```literal
    assembly_reference_sequence  TGTGAGTAATATTTGATGTTGAATAGTAAATGAAATATATTCCTATTTAAGGGGGGGTGTCCAAGACACCATTAAAACTTATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAATACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACCATAATATACACAATTTTAATCAAGTAAACTCTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAATCAATTCAAAACCCTTAAAAATAAAATAAGAAACAAC
    SLH_AL_0013                  TGTGAGTAATATTTGATGTTGAATAGTAAATGAAATATATTCCTATTTANGGGGNGNTGTCCAAGACNCCATTAAAACNTATGTTCAAATCCTTTCATTTCATAAAGCTAANAATGAANACCNCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACCATAATATACACAATTTTAATCAAGTAAAMTCTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGSAAAACACAAAATCAATTCAAAACNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0042                  TGTGAGTAATATTTGATGTTGAATAGTAAATGAAATATATTCCTATTTAGGGGGGGGTGTCCAAGACACCATTAAAACTTATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAATACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACCATAATATACACAATTTTAATCAAGTAAACTCTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAATCAATTCAAAACCCTTTAAAATAAAATAAGAAACAAC
    SLH_AL_0063                  TGTGAATAATATTTCATGTTGAATAGTAAATGAAATACATTCCTATTTAAGGGGAGCTGTCCAAGACTCCATTAAAACATATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAACACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACTATAATATCCACAATTTTAATCAAGTAAACACTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAATCAATTCAAAACCCTTTAAAATAAAATAAGAAACAAC
    SLH_AL_0084                  TGTGAGTAATATTTGATGTTAAATAGTAAATGAAATAAATTCCTATTTAAGGGGGGGTGTCCAAGACACCATTAAAACTTATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAACACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACCATAATATACACAATTTTAATCAAGTAAACTCTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAATCAATTCAAAACTTTTTAAAATAAAATAAGAAACAAN
    SLH_AL_0104                  TGTGAGTAA-ATTTGATGTTGAATAGTAAATGAAATATATTCCTATTTAAGGGGGGGTGTCCAAGACACCATTAAAACTTATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAATACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACCATAATATACACAATTTTAATCAAGTAAACTCTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAACCAATTCAAAACCCTTTAAAATAAAATAAGAAACAAC
    SLH_AL_0105                  TGTGAATAATATTTCATGTTGAATAGTAAATGAAATACATTCCTATTTAAGGGGAGCTGTCCAAGACTCCATTAAAACATATGTTCAAATCCTTTCATTTCATAAAGCTAAGAATGAACACCTCAGTTACTAACATAACAATATAACATGCATATAGCCCACATCTATGCCCAGAATCAGACTATAATATCCACAATTTTAATCAAGTAAACACTAAAAATTCATTCATAAGAATACTATCAAGCAAGCAATAGCAACGCAAAACACAAAATCAATTCAAAACCCTTTAAAATAAAATAAGAAACAAC
    //                                *        *     -                *           -    * *          *          *                                       *                                                               *       *                    -*                                              -           -           --  -                    |0:A_tuberculatus_Chr01:17883-18190
    assembly_reference_sequence  GCCTTGCCTATAATTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGGTTGATAAAACCTCTTGCTCTAGTCAAATTGCTCACATGTGCAATGATGTGTGTCAAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCCTTCTTGCTGCCATATATAGTTTTCTAGGTCATAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAAATTTTGAAACACTGCAA
    SLH_AL_0063                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGCGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_0064                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGTGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_0100                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGTGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_0101                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGTGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_0104                  GCCTTTCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTTCTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGTGTGTCNAATATGACAGACAAAACAAGGCTCTT-------ATTATGTGGCCTTCTTGCTGCCATGTATAGTTTTCTAGGTCNTAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAAATTTAGATACACTGCAA
    SLH_AL_0105                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGCGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_0106                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGTGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    SLH_AL_3066                  GCCTTGCCTATAGTTCTCCTAAACTTTTACACTCTTTGTGGGTACTATCATGTTATTTGGCTTGATAAAACCTCTTGGTCTAGTCAAACTGCTCACATATGCAATGATGCGTGTCTAATATGACAGACAAAACAAGGCTCTTTTAGCTCATTATGTGGCTTTCTTGCTGCCATGTATAGTTTTCTAGGTCNAAGTGACCCAAATGAAGATTTCTATGTTGCATTCTCACATTGCCATAAAAGATTTAGAAACACTGCAA
    //                                -      -                              -                -                -          -         -          *     -                                           *             -                 *                                                 *    -  -         |1:A_tuberculatus_Chr01:98725-98983
    assembly_reference_sequence  ATTCATGTTTCTTGTTTTTCTCTGATAAACCCTGTTTCTTCTAGTAAACCCACAAGCAAAGCTCCATTTCCTTAACCTCCTCCCAAACAAACACTCCCAACTTATTGTTTCTTCATCATATGTGCCAATCAATTCTTTGTCCAAAAACCTATTTAAAAACCATTACCATATTCATTACGTATGCCCATTGATATTCGATTC
    SLH_AL_0063                  NTTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_0064                  ATTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_0084                  ATTCATGTTTCTCGTTTTMCTCTGATAANCCCTGTTTCTTCRAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCAWATGTGCCAATCAATTCTTTGTCCNAGAACCCATATAAAAACCATTACCATATTCATTACRTATGCCCATTGATATTCAATTC
    SLH_AL_0100                  ATTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_0101                  ATTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_0104                  NNNNNNNNNNNNNGTTTTCCTCTGATAANCCCTGTTTCTTCAAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAAMACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCNANAACCNATNTAAAAACCATTACCATATTCATTACGTATGCCCATTGATATTCAATTC
    SLH_AL_0105                  NTTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_0106                  ATTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    SLH_AL_3065                  NNNNNNNNNNNNNGTTTTCCTCTGATAANCCCTGTTTCTTCAAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCNANAACCNATNTAAAAACCATTACCATATTCATTACGTATGCCCATTCATATTCAATTC
    SLH_AL_3066                  NTTCATGTTTCTCGTTTTCCTCTGATAANCCCTGTTTCTTCGAGTAAACCCACAAGCAAAGCTCCATTTCCTTANCCTCCTCCCAAACAAACACTCCCANCTTATTGTTTGTTCATCATATGTGCCAATCAATTCTTTGTCCAAGAACCCATATAAAAACCATTACCATATTCATTACATATGCCCATTGATATTCAATTC
    //                                       -     -                      *                                                 -                  -       -                         -    -  -                         *          -      -    |2:A_tuberculatus_Chr01:101256-101456
    assembly_reference_sequence  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACCGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGAGCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCGGTGATGGAAGTATTGATG
    SLH_AL_0012                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0014                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0018                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGYTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGYGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0030                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0034                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0036                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0042                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATRTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0048                  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACCGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATAGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0063                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_0064                  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACAGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGTTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_0100                  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACAGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGTTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_0101                  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACAGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGTTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_0104                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTTGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_0105                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_0106                  ATGTGTTCAGGTACAATATTTTGATATTTAGGCTTTTAAAGGAAAACAGGCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGTTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    SLH_AL_3065                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATYTGTTGTTTGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGCNGTGATGGAAGTATTGATG
    SLH_AL_3066                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTAAAGAGGTTTGTGATTTGTTGTTGGATATGAAGGAGANCCAAATCGTGCCCGATAAGGTGACCATGAATGCTGTATTGTGCTTCTTTTGTAAGGCTGGAATGATGGATGTTGTAGTTGACTTGTATGAAGACAAGGCTGAATTTGGGCTTACTCCACATGGTTTGGTGTTTAATTTTGTTATGAACACTCTATGTNGTGATGGAAGTATTGATG
    //                                                                          *            -       -       *                                  -                                                 -                                      *           *                                  *                   |3:A_tuberculatus_Chr01:102530-102796
    assembly_reference_sequence  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTGGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGAAGTGCAATGCTTGGCCAGGTTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTACTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGATG
    21040XD-01-08_S40_L002       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTTATAATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0012                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0014                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0018                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0030                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0036                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAMGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0042                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0048                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTKGGGCGTTCGATTTCATCGTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCACAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACTTGGAAGAACTTGTTCAAGAATGCTTTCCACTTGATCTCTANTCTTACAACATGATAATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0063                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACTTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACGTCGANG
    SLH_AL_0064                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATGATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACTTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0100                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATGATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACTTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0101                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATGATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACTTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0105                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACCTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_0106                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATGATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACTTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    SLH_AL_3066                  AAGTATGATGAGGCTGTGCTATTAATCAAAGAAATGGAGAAAGTTGGGCGTTCGATTTCATCTTTTGTTGGTAATACACTTTTGTTGAACTCACTTTATAATGACAACTTGCAGAAGGCTTGGCTTCGTTTGAGAGATCCATCTAATGAAACACCTAGNAGTGCAATGCTTGGCCAGGNTATTGCTGTTTTTTCTGGCCAAATAAATAGTGATGTTAAGGCGGAGGACATGGAAGAACTTGTTCAAAAATGCTTTCCACCTGATCTGTANTCTTACAACATGATTATGAGAAAATTGAGCATCAATAACATCGANG
    //                                                                       *                 *                                    *             *                                                                                                       -          *                 *            *      *                 *                        -      |4:A_tuberculatus_Chr01:103624-103939
    assembly_reference_sequence  GTGGTGAGAAAGAAAGGATCCTAGTCTTATGATATCCTACATGGCATCATACTATTAATACACAAAACCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTATCAACACCAACTATTCTGTTCATGTCACTTTCAGTTGTGATCAAGTTATAAATTCACTTTTTAATTTGA
    21040XD-01-07_S39_L002       GTGGTGNNAAAGAAAGGATCCTAGTCTTATGNTATCCTACATGGCATCATANTATNNATACACNAAANCTCTCTTGTCCAAGTAATTATCTGACTNNNNNNNNNNNNNNNNNNNAACTANTNTGTTNATGTCACTTTCANTTGTGNTCAAGTNATAAATTCACTTTTNNNNNNNN
    SLH_AL_0014                  GTGGTGAGAAAGAAAGGATCCTAGTCTTATGRTATCCTACATGGCATCATACTATTAATACACAAAACCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTATCAACACCAACTAATCTGTTCATGTCACTTTCAGTTGTGATCAAGTTATAAATTCACTTTTTAWTTTGA
    SLH_AL_0034                  GTGGTGAGAAAGAAAGGATCCTAGTCTTATGGTATCCTACATGGCATCATACTATTAATACACCAAACCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTATCWACACCAACTAATCTGTTCATGTCACTTTCAGTTGTGATCAAGTTATAAATTCACTTTTTATTTTGA
    SLH_AL_0042                  GTGGTGAGAAAGAAAGGATCCTAGTCTTATGGTATCCTACATGGCATCATACTATTAATACACAAAACCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTATCAACACCAACTAATCTGTTCATGTCACTTTCAGTTGTGATCAAGTTATAAATTCACTTTTTATTTTGA
    SLH_AL_0104                  GTGGTGGGAAAGAAAGGATCCTAGTCTTATGGTATCCTACATGGCATCATANTATTGATACACAAAACCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTNTCAACACCAACTAATTTGTTTATGTCACTTTCAATTGTGTTCAAGTCATAAATTCACTTTTTATTTTGA
    SLH_AL_0105                  GTGGTGGGAAAGAAAGGATCCTAGTCTTATGGTATCCTACATGGCATCATANTATTAATACACAAAAGCTCTCTTGTCCAAGTAATTATCTGACTCAGCAAATCTNTCNACACCAACTAATTTGTTCATGTCACTTTCAATTGTGTTCAAGTCATAAATTCACTTTTTANNNNNN
    //                                 *                        *                        -      -   -                                        -          - *    -            *     *      *                *     |5:A_tuberculatus_Chr01:110679-110853
    assembly_reference_sequence  CTAGAGAGAGAAAGTGATGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGGAAGTGTAGTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTATGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTAGAAAATCTCATCACAACAACAACAACAACAACAATGGAAACGGAAAACCGAAAACCGATTTTTCCAACAAAGGAAGTCAAAACGACATGAAAATTGTCGTAATTATGG
    21040XD-01-07_S39_L002       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAACANNAATGGAAANNNNAAACCNAAAACNGANNTTNNCNACAAANNAANTCAAANCGACNTNAAANTNGTCGTAATTATGG
    SLH_AL_0012                  CTAGAGAGAGAAAG---TGATGATGATGAAGGAAGTAATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAAATGGAAGT---GTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTAGAAAATCTCATC------ACAACAACAACAACAATGGAAACGGAAAACCAAAAACCGATTTTTCCAACAAAGGAAGTCAAANCGACATGAAANTTGTCGTAATTATGG
    SLH_AL_0014                  CTAGAGAGAGAAAG---TGATGATGATGAAGGAAGTRATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAAATGGAAGT---GTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTRCTTTGCTTATTCTWGCTTGTTCTTTTASAAAATCTCATCACA---ACAACAACAACAACAATGGAAACGGAAAACCAAAAACCGATTTTTCCAACAAAGGAAGTCAAANCGACATGAAANTTGTCGTAATTATGG
    SLH_AL_0018                  CTAGAGAGNGAAAG---TGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTANATGGAAGTGTAGTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTNGCTTGTTCTTTTAGAAAATCTCATC------ACAACAACAACAACAATGGAAACGGAAAACCNAAAACCGATTTTTCCAACAAAGGAAGTCAAANCGACATGAAANTTGTCGTAATTATGG
    SLH_AL_0030                  CTAGAGAGAGAAAG---TGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGGAAGT---GTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTNGCTTGTTCTTTTAGAAAATCTCATCACA---ACAACAACAACAAYAATGGAAACGGAAAACCAAAAACCGATTTTTCCNACAAAGGAAGTCAAANCGACATGAAANTTGTCGTAATTATGG
    SLH_AL_0034                  CTAGAGAGAGAAAG---TGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTARATGGAAGT---GTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTWGCTTGTTCTTTTAGAAAATCTCATCACA---ACAACAACAACAACAATGGAAACGGAAAACCRAAAACCGATTTTTCCAACAAAGGAAGTCAAANCGACWTGAAANTTGTCGTAATTATGG
    SLH_AL_0048                  CTAGAGAGAGAAAG---TGATGATGATGAAGGAAGTNATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTANATGGAAGT---GTTCAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTNGCTTGTTCTTTTAGAAAATCTCATC------ACAACAACAACAACAATGGAAACGGAAAACCAAAAACCGATTTTTCCAACAAAGGAAGTCAAANCGACATGAAANTTGTCGTAATTATGG
    SLH_AL_0063                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACT---------------AACAACACAATGTCAATGGTAGATGGAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTATGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTCTTCTTTTACAAAATCTCATNACAACAACAACAACAACAACAATGGAAACCCAAAACCAAAAACGGATTTTTCCAACAAACAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_0064                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGAAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACAAAATCTCATNACAACAACAACACCAACAACAATGGAAACCCAAAACCTAAAACGGATTTTTCCAACAAAGAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_0084                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGGAAGT---GTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTATGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACTAAATCTCATNACAACAACAACNNCAACAACAATGGAAACCCAAAACCAAAAACGGATNTTTNCAACAAAGGAACTCAAANCGACATNAAANTTGTCGTAATTATGG
    SLH_AL_0100                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGAAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACAAAATCTCATNACAACAACAACACCAACAACAATGGAAACCCAAAACCTAAAACGGATTTTTCCAACAAAGAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_0101                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGAAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACAAAATCTCATNACAACAACAACACCAACAACAATGGAAACCCAAAACCTAAAACGGATTTTTCCAACAAAGAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_0104                  CTAGAGAGAGAAAGTGATGATGATGATGAAGGAAGTGATGASTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGGAAGT---GTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACWAAATCTCATNACAACAACAACMACAACAACAATGGAAACCCWAAACCAAAAACGGATCTTTCCAACAAAGGAASTCAAANCGACATSAAANTTGTCGTAATTATGG
    SLH_AL_0105                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACT---------------AACAACACAATGTCAATGGTAGATGGAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTATGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTCTTCTTTTACAAAATCTCATNACAACAACAACAACAACAACAATGGAAACCCAAAACCAAAAACGGATTTTTCCAACAAACAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_0106                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACTAACACCACAAAAACAAACAACACAATGTCAATGGTAGATGAAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTNTGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTGTTCTTTTACAAAATCTCATNACAACAACAACACCAACAACAATGGAAACCCAAAACCTAAAACGGATTTTTCCAACAAAGAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    SLH_AL_3066                  NTAGAGAGAGAAANTGATGATGATGATGAAGGAAGTGATCACT---------------AACAACACAATGTCAATGGTAGATGGAAGNGTAGTTNAACAAGGATGGAAAACACCAATACCATATCTCTTTGGTGGACTTGCACTTATGCTATGTCTTATTGCTATTGCTTTGCTTATTCTTGCTTCTTCTTTTACAAAATCTCATNACAACAACAACAACAACAACAATGGAAACCCAAAACCAAAAACGGATTTTTCCAACAAACAAACTCAAANCGACATGAAANTCGTCGTAATTATGG
    //                                                               *  - -                                     *   *                                                                                  -             *    *        **                     -*      -         **-     *     *   -           **  *          - -     *             |6:A_tuberculatus_Chr01:116394-116695
    assembly_reference_sequence  CTTTTCCTCTTAATTATTAGTTTTCGAAGTACATAATGAGTGATCACTAAGTCAGCACTAACAAGTTATTGGCATATCTTCATAAATATAACCCTATTTGGACTATAACTTATATTTTAAAGTGCTTTTATATATTGAGTATATAATTTTGTACATGATATGTGTTAATGATTGCTAGCTTAAGTCATTTAACTAAAATAGAACAATGCATATACCATTACAGGAAGGGTTAGCCATACGATATTTA
    21040XD-01-08_S40_L002       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNTAATGAGTGATCACTAAGTNANCACTAACAANTTATTGGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTATNACTTATNTTTTNAAGTGNTTTTATATATTNAGTATATNATTTTNTANANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTANAACAATNCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0063                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTAAAGTGNTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAGAACAATCCATATACCATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_0064                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTGAAGTGNTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAAAACAATCCATATACGATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_0100                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTGAAGTGCTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAAAACAATCCATATACGATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_0101                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATNTTTTGAAGTGNTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAAAACAATCCATATACGATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_0105                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTAAAGTGCTTTTATATATTGAGTATATAATTTTGTAC-TGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAGAACAATCCATATACCATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_0106                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTGAAGTGNTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAAAACAATCCATATACGATTACATGAAGGGTTAGCCATACAATATTTA
    SLH_AL_3066                  CTTTTCCTCTTAATTATTAGTATTCGANGTACATAATGAGTGATCACTAAGTCACCACTAACAAGTTATTGGCATATCTTNATAAATATAACCCTATTTGGACTATAACTTATATTTTAAAGTGNTTTTATATATTGAGTATATAATTTTGTANATGATATGTGTTAATGATTGCTAGCTTAAATCATTTAACTAAAATAGAACAATCCATATACGATTACATGAAGGGTTAGCCATACAATATTTA
    //                                                -                                -                                                               *                                                                -                *      -       *      -                -       |7:A_tuberculatus_Chr01:164410-164656
    assembly_reference_sequence  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGAAAATGAAGTCAGAAATTAGAATAGTAGTACAAGAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0012                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTACGTGGAAAGAWGTAGAAGAAAATGAAGTCAGAAATTAGAATAGTAGTACAANAANATTAGATAAATTAACAGGGTTT
    SLH_AL_0014                  TTAGAGAGATAAACATGGAGATTTTTGTAAWTAGAAAANAGATGGG--TTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAARCGTACGTGGAAAGAWGTAGAAGAWAATGAAGTCAGAAATTAGAATAGTAGTTCAANAATATTAGATAAATWAACAGGGTTT
    SLH_AL_0018                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTRCGTGGAAAGATGTAGAAGATAATGAAGTCAGAAATTAGAATAGTANNNNAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0030                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTACGTGGAAWGAWGTAGAAGAWAATGAAGTCAGAAATTAGAATAGTAGTACAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0034                  TTAGAGAGATAANCATGGAGATNTTTGTAANTAGAAAANAGATGGNTTTTTTTGTNGGGNAAGACNAAATTTNGANTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTACGTGGAAAGAWGTAGAAGAWAATGAAGTCAGAAATTAGAATAGTANNNNAANAANATTRGATAAATTAACAGGGTTT
    SLH_AL_0036                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAGCGTACGTGGAAAGATGTAGAAGATAATGAAGTCAGAAATTAGAATAGTANTNNAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0042                  TTAGAGAGATAAACATGGAGATTTTTGTAATTAGAAAANAGATGGG--TTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGAGAAACGTACGTGGAANGAAGTAGAAGAAAATGAAGTCAGAAATTAGAATAGTANTNNAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0048                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATWKTGKRRAATTAGGCAAAATTTAGGGAAAATATGGAAAGAGANAAGCGTACGTGGAAAGAWGTAGAAGNNAATGAAGTCAGAAATTAGAATAGTANNNNAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0063                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGNTTTTTTTGTTGGGGAAGACCAAATTTAGATTATTGTGGAAAATTAGGCAAAATTTGGAGAAAATATGGAAAGAGAGAAGCGTATGTGGAAAGAAGTAGACGAGAATGAAGTCAGAAATTCGAATAGTAGTAGATNAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0064                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAATAGATGGNTTTTTTTGTGGGGGAAGACAAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGANAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGGAAATGAAGTCAGAAATTAGAATAGTACTAGAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0084                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGT-TTTTTGTTGGGKAAGACMAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGRGAAAAWATGGANAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGRAAATGAMGTCAGAAATTAGAATAGTASTASAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0100                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAATAGATGGNTTTTTTTGTGGGGGAAGACAAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGANAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGGAAATGAAGTCAGAAATTAGAATAGTACTAGAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0101                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAATAGATGGNTTTTTTTGTGGGGGAAGACAAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGANAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGGAAATGAAGTCAGAAATTAGAATAGTACTAGAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0104                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGAT------------TTAGGCAAAATTTAGAGAAAAAATGGAAAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGAAAATGACGTCAGAAATTAGAATAGTAGTAGAANAATATTAGATAWATTAACAGGGTTT
    SLH_AL_0105                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGNTTTTTTTGTTGGGGAAGACCAAATTTAGATTATTGTGGAAAATTAGGCAAAATTTGGAGAAAATATGGAAAGAGAGAAGCGTATGTGGAAAGAAGTAGACGAGAATGAAGTCAGAAATTCGAATAGTAGTAGATNAATATTAGATAAATTAACAGGGTTT
    SLH_AL_0106                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAATAGATGGNTTTTTTTGTGGGGGAAGACAAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGGGAAAATATGGANAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGGAAATGAAGTCAGAAATTAGAATAGTACTAGAANAATATTAGATAAATTAACAGGGTTT
    SLH_AL_3065                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGGTTTTTTTGTTGGGGAAGACCAAATTTTGATTATTGTGGAAAATTAGGCAAAATTTAGAGAAAAAATGGAAAGAGAGAAGCGTNCGTGGAAAGAAGTAGAAGAAAATGACGTCAGAAATTAGAATAGTAGTNNANNAANATTAGATAAATTAACAGGGTTT
    SLH_AL_3066                  TTAGAGAGATAAACATGGAGATTTTTGTAAATAGAAAAAAGATGGNTTTTTTTGTTGGGGAAGACCAAATTTAGATTATTGTGGAAAATTAGGCAAAATTTGGAGAAAATATGGAAAGAGAGAAGCGTACGTGGAAAGAAGTAGAAGAGAATGAAGTCAGAAATTCGAATAGTAGTAGATNAATATTAGATAAATTAACAGGGTTT
    //                                                         *       *                *   -     *      *      --  ---               * *     *              *   -*      -  *     * **     *          *        * -* *       -    -  -          |8:A_tuberculatus_Chr01:179827-180032
    assembly_reference_sequence  TCGCTTTTCATTCATGCGGTAAGAAAATGATAGTGGGACCCACCAATGAGTCAATAACCAATCTACTAGTAGAGATAGAGTAACCTTAGACCAATTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCTACCATTTTCATATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCAGTACAACAAGCTTCCACATCAGCTTCAACACAGGCCCGGCCTATGTGCAGGAGCTCCCTAAGCAACAGCTCAGGGCCTCACTTTGTTGGGGGCCCAGTTTTTCAAGGGGCCCTTATTTAAATTTTTGCTCCGGGCCTCGATATTGTCCAAGACGGCTCTGCTTGAACATCGCAAACAATCATAAAACTCTCCAAGGTTTCATCCACACATCACATTGGTCTCCATTTGCTACACATTGGTAAGGTCTAAAATACACATAGGAAGTTTAAACCCTTTATTTCATTTTTGATGATATAGGGATTAGGGAGAGTATTATAATTCTTTCAATAGTGTAATATATGTTGATATATCACACTAATTAATTGTCTATAACATATTTGTTGGTGGTACTATT
    SLH_AL_0012                  TCGCTTTTCATTCGTGCGGTAAGAAAATGAGAGTGGGACCCACCAATGAGTCAATAACCAATCTACTAGTAGAGATNGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACRAGCTGCTAGAACCTANAACACACACCAGATCTTCTACCATTTTCATATANATAAAAAAATGTCATAAATCTTANCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCAGTACAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0014                  NCGCTTTTCATTCGTGCGGTAAGAAWATGAGAGTGGGACCCACCAATGAGTCAATAACCAATCTACTAGTAGAGATAGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACGAGCTNCTAGAACCTANAACACACACCAGATCTTCWACCATTTTCATATANATAAAAAAATGTCATAAATCTTANCTTTTGCTACGYTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCASTACAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCTATGTGCAGGAGCTCTCTAAGCAGCAGCTCAGGGCCTCACTTTGTTGGGGGCCCAGTTTTTCAAGGGGCCCNTATTTAAACTTTTGCTCCGACCCTCGATATTGTCCAAGACGGCACTGCTTGAACATCGCAAACAATCATAAAACTCTNCAAGNTTTCNTCCACACATCACATTGGTCTCCATTTGCTACACATTGGTAAGGTCNAAAATACACATAGGAAGTTTAAACCCTTTATTTCATATTTGATTATATAGGGATTAGNNNGAGTATTATAATTCTTTCAATAGTGTAATATATGTTGATATNTCACACTAATTAATTGTCTATNAC--ATTTGTTGGTGTTACTATT
    SLH_AL_0018                  NCGCTTTTCATTCGTGCGGTAAGAAAATGAGAGTGGGACCCACCAATGAGTCAATAACCAATCTACTAGTAGAGATAGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACAAGCTNCTAGAACCTANAACACACACCAGATCTTCTACCATTTTCATATAA----------------------ANCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCAGTTCAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0030                  TCGCTTTTCATTCGTGCGGTAAGAAAATGAGAGTGGGACCCACCNNNNNNNNNNNNNNNAATCTACTAGTAGAGATWGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTANAACACACACCAGATCTTCTACCATTTTCATATAA----------------------ANCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCAGTACAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTCAAGGGGCCCNTATTTAAACTTTTGCTCCGACCCTCGATATTGTCCAAGACGGCTCTGCTTNAACATCGCAAACAATCATAAAACTCTTCAAGATTTCGTCCACACATCACATTGGTCTCCATTTGCTACACATTGGTAAGGTCNAAAATACACATAGGAAGTTTAAACCCTTTATTTCATATTTGATGATATAGGGATTAGNNNGAGTATTATAATTCTTTCAATAGTGTAATATATGTTGATATNTCACACTAATTAATTGTCTATNACATATTTGTTGGTGTTACTATT
    SLH_AL_0036                  TCGCTTTTCATTCGTGCGGTAAGAAAATGAGAGTGGGACCCACCAATGAGTCAATAACCAATCTACTAGTAGAGATNGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACRAGCTGCTAGAACCTANAACACACACCAGATCTTCTACCATTTTCATATANATAAAAAAATGTCATAAATCTTANCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCAGTACAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0042                  TCGCTTTTCATTCGTGCGGTAAGAAAATGAGAGTGGGACCCACCNNNNNNNNNNNNNNNAATCTACTAGTAGAGATAGAGTAACCTTAGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTANAACACACACCAGATCTTCTACCATTTTCATATANATAAAAAAATGTCATAAATCTTANCTTTTGCTACGCTTTTCCCACCTTAAATTTCAATACTCTTTTCTCTTAAATCASTACAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTCAAGGGGCCCNTATTTAAACTTTTGCTCCGACCCTCGATATTGTCCAAGACGGCTCTGCTTNAACATCGCAAACAATCATAAAACTCTTCAAGATTTCGTCCACACATCACATTGGTCTCCATTTGCTACACATTGGTAAGGTCNAAAATACACATAGGAAGTTTAAACCCTTTATTTCATATTTGATGATATAGGGATTAGNNNGAGTATTATAATTCTTTCAATAGTGTAATATATGTTGATATNTCACACTAATTAATTGTCTATNACATATTTGTTGGTGTTACTATT
    SLH_AL_0048                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCCTCGATATTGTCCAAGACGGCTCTGCTTG------------AATCATAAAACTCTTCAAGATTTCGTCCACACATCACATTGGTCTCCATTTGCTACACATTGGTAAGGTCGAAAATACACATAGGAAGTTTAAACCCTTTATTTCATATTTGATGATATAGGGATTAGNNNGAGTATTATAATTCTTTCAATAGTGTAATATATGTTGATATNTCACACTAATTAATTGTCTATNACATATTTGTTGGTGTTACTATT
    SLH_AL_0064                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCNACCATTTTCNNATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCATAAATTTCAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0084                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCNACCATTTTCNNATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCATAAATTTCAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0100                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCNACCATTTTCNNATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCATAAATTTCAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0101                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCNACCATTTTCNNATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCATAAATTTCAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0106                  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACCAANTAAATTGGTTGCTTTTAATATACAAGCTGCTAGAACCTAGAACACACACCAGATCTTCNACCATTTTCNNATAAATAAAAAAATGTCATAAATCTTATCTTTTGCTACGCTTTTCCCACCATAAATTTCAATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    //                                        -           -    -                                             -                                         *                                  -                                                  -          *                              * -                                                     -       -                                                        -          -*                      -                                 -    -    -                                             -                                    -      -                                                                                               -       |9:A_tuberculatus_Chr01:190356-190995
    assembly_reference_sequence  CCTGATGAGGAGCAGAAGAATTCAGTGTTCCCTCACCTAGACGATAACTCCCAGGCCAGTAATCAACATAGCAGCTGGAGTTCATATAAGGTAAAACAGGTTGCCCTACAATGATCTGGTGTGATGCAAGACCATGGTTATTAACATAGGCAAAAGCTGGCATCATCGTTAATCCAGCCATGATCATCATTCTTGATTCTTGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    SLH_AL_0012                  CCTGATGAGGAGCAGAAGAATNCAGTGTNCCCTCACNNAGANGATAACTCCCAGGCCANTAATNANNATAGCANCTGGNGTTCATATAAGGTAAAACAGGTTGNCCTACAATGATCTGGTGNGNTGCAAGANCATGGTTATTAACATAGGCAAAAGCTGGCATCATCGTTAATCCAGCCATSATCATCATTYTTGATTCATGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    SLH_AL_0018                  CCTGATGAGGAGCAGAAGAATCCAGTGTTCCCTCACCTAGACGATAACTCCCAGGCCAGTAATTAACATAGCAGCTGGTGTTCATATAAGGTAAAACAGGTTGCCCTACAATGATCTGGTGTGATGCAAGACCATGGTTATTAACATAGGCAAAAGCTGGCATCACNGTTAATCCAGCCATCATCATCATTCTTGATTCATGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    SLH_AL_0034                  CCTGATGAGGAGCAGAAGAATNCAGTGTNCCCTCACNNAGANGATAACTCCCAGGCCANTAATNANNATAGCANCTGGNGTTCATATAAGGTAAAACAGGTTGNCCTACAATGATCTGGTGNGNTGCAAGACCATGGTTATWAACATAGGCAAAAGCTGGCATCATCGTTAAWCCAGCCATCATCATCATTCTTGATTCWTGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    SLH_AL_0042                  CCTGATGAGGAGCAGAAGAATNCAGTGTNCCCTCACNNAGANGATAACTCCCAGGCCANTAATNANNATAGCANCTGGNGTTCATATAAGGTAAAACAGGTTGNCCTACAATGATCTGGTGNGNTGCAAGACCATGGTTATTAACATAGGCAAAAGCTGGCATCATCGTTAATCCAGCCATGATCATCATTCTTGATTCATGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    SLH_AL_0048                  CCTGATGAGGAGCAGAAGAATCCAGTGTTCCCTCACCTAGACGATAACTCCCAGGCCAGTAATCAATATAGCAGCTGGTGTTCATATAAGGTAAAACAGGTTGCCCTACAATGATCTGGTGTGATGCAAGACCATGGTTATTAACATAGGCAAAAGCTGGCATCATCGTTAATCCAGCCATGATCATCATTCTTGATTCTTGATTTCTTATATTTGCTTCTCTTCTCTCTTGTTTATGAGCATTTTGATGCCCTCCTA
    //                                                -                                         -  -           -                                                              -                       -      -        *         -       *                                                          |10:A_tuberculatus_Chr01:202773-203030
    assembly_reference_sequence  ACCAAGCTGAGATGAGAGACCAGACATAAAGGCTTCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGTAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    21040XD-01-08_S40_L002       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTCTCAACAGATTGTGTATCAGCTGATNNNNNNNNNNNNNNNNNNNNNNN
    21040XD-01-09_S41_L002       ACCNAGCTGAGATGAGAGACCAGACATANANGCNTCNAGNGCATCTCCAGNNTCTTNTTCAANCNTAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCNTTCTTATCCAGAAGAGAATCGGCNGTCTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    SLH_AL_0014                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCTTCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGTAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    SLH_AL_0034                  ACCRAGCTGAGATGAGAGACCAGACATAAAGGCTTCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGTAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    SLH_AL_0036                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCTTCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGTAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    SLH_AL_0042                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCTTCCAGTGCATCTCCAGAA---TCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAKAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGTAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    SLH_AL_0064                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCATGTTTCTGTTTCGT
    SLH_AL_0100                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCATGTTTCTGTTTCGT
    SLH_AL_0101                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCATGTTTCTGTTTCGT
    SLH_AL_0104                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTYTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCAGGTTTCTGTTTCGT
    SLH_AL_0105                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTTGCAGGTTTCTGTTTCGT
    SLH_AL_0106                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCCGTCTCAACAGATTGTGTATCAGCTGATTTGTTAGCATGTTTCTGTTTCGT
    SLH_AL_3066                  ACCGAGCTGAGATGAGAGACCAGACATAAAGGCATCCAGTGCATCTCCAGAATCTTCTTCAACCTTAGGCTCTGATGTCTTGCTCTTCTCTTCCAGAATCAATCTTTTTTTCTCTTCCATTTCATTCAGCAACGCCTCCTTCTTATCCAGAAGAGAATCGGCAGTCTCAACAGATTGTGTATCAGCTGATTTGTTTGCAGGTTTCTGTTTCGT
    //                              *                             *                                                       -     -                                                      *           -                                *   *             |11:A_tuberculatus_Chr01:211536-211748
    assembly_reference_sequence  TTTTTTACAATTGCTACGCTTGCATACTAGCCATTATCCATTAAGCTATATTGTCTTACATCTTGTTCTTCATGTATAAAGAAAACAGAGGAAGTACTAATTTGCGAATTGATAACAAGACACACTAATGATGAAAACAC
    SLH_AL_0063                  TTTTTTACAATTGCTACACTTGCATACTAGCCATTACCCATCAAGCTATATTGTCTAACATCTTATTCTTCATGTATAAAGAAAACAGAGGAAGTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAG
    SLH_AL_0084                  TTTTTGACAATTGCTACACTTGCATACTAGCCATTATCCATCAAGCTATATTGTCTTACATCTTATTCTTCATGTATAAAGAAAACAGAGGAAGTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAT
    SLH_AL_0104                  TTTTTTACAATTGCTACACTTGCATACCAGCCATTATCCATCAAGCTATATTGTYTTACATTTTATTCTTCAYGTATAAAGAAAACAGAGGAAGTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAT
    SLH_AL_0105                  TTTTTTACAATTGCTACACTTGCATACTAGCCATTACCCATCAAGCTATATTGTCTAACATCTTATTCTTCATGTATAAAGAAAACAGAGGAAGTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAN
    SLH_AL_3065                  TTTTTTACAATTGCTACACTTGCATACTAGCCATTATCCATCAAGCTATATTGTCTTACATTTTATTCTTCATGTATNNNNNNNNNNNNNNNNNTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAN
    SLH_AL_3066                  TTTTTTACAATTGCTACACTTGCATACTTGCCATTACCCATCAAGCTATATTGTCTTACATCTTATTCTTCATGTATNNNNNNNNNNNNNNNNNTACTACTTTGCGAATTGATAACAAGACACACTAATGATGAAAACAG
    //                                -           -         --       *    -            - *    *  -       -                          -                                       *|12:A_tuberculatus_Chr01:216391-216530
    assembly_reference_sequence  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAATTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAATGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACAACTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTCCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    21040XD-01-09_S41_L002       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGCNNANAGNANNAGNNGGGAATAAATTNCNAGTTACCTCCGTAATCTACAACNGTANNNNNNNNNNNNNNNNNN
    SLH_AL_0012                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAATGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTCCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0014                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAAYGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGYATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTYCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0018                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAAYGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGYATCAATATGAGGTTAGCTASATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTYCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0030                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAATTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAATGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCAGCATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTTCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0034                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAAYGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCRGCATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGYGGGGAATAAATTTCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0036                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAATTTATCCCTCTCAAGATATGCTTCCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCGAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGTATCAATACGAGGTTAGCTANNTGAACCAAAACANCTGAACGTNATTCAAAAGCTAAAAGCATGAGTCGGGAATAAATTTCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0042                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAGTCGGCATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCWWAAAGMAAAAGTGGGGAATAAATTTCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0048                  ACATAAGTTCTCCCTGTTTCTTCAAGTATCCCTGCCGCTGCTGATGCATAGTTTATCCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAAYGGCAATCCCAGATGTCGGGCTTCAATTTTGCACCAAAAAGAAATCAAKTCGGCATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCAAAAGCTAAAAGCAAAAGTGGGGAATAAATTTCAAGTTACCTCCGTAATCTACAACAGTAAGACCATTGCAAAACCTC
    SLH_AL_0064                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATGCCTCTCAATATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCAAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_0100                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATGCCTCTCAATATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCAAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_0101                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATGCCTCTCAATATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCAAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_0104                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATTCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCGAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACATGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_0105                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATGCCTCTCAATATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCGAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCATAAGCTAACAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_0106                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATGCCTCTCAATATATGCTTTCCGAACGATGTTATTGAGAAGTATGGTGGTATCAACGGCAATCCAAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACGTGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    SLH_AL_3065                  ACGTAAGTTCTCCCTGTTTCTTCAAGTATTCCTGCCGCTGCTGATGCATAATTTATTCCTCTCAAGATATGCTTTCCGAACGATGTTATTGAGAAATATGGTGGTATCAACGGCAATCCGAGATGTCGTGCTTCAATTTTGCACCAAAAAGAAATCGAGTCGGTATCAATATGAGGTTAGCTACATGAACCAAAACANCTGAACATGATTCATAAGCTAAAAGCAAAAGTGGGGAATAAATTACAAGTTACCTCCGTAATCTACAACCGTAAGACCATTGCAAAACCTC
    //                             *                          *                    *     *        *        -                    -              *        *        *                           * -  * *       -           -                    *       *    -- -  - --  --           *                        *                     |13:A_tuberculatus_Chr01:274089-274377
    assembly_reference_sequence  TTAAGCTGATGGTTAAAGCCCTAGGACATGTTATATAATCTATCATTTCTCGGACCAACTAGTAACCAAAGTAGTAGTTAATATCGATAGTAAAGTTAGATGGACATAAATCTATGTTGTCTGATTTCGTAGTGTCATTGACAGCAGACTAGAAGCTGTCTCGCCATTGGCGAGCTTGATCTGATAAATCGAAATGAGTACACCTGGTTCTAAAACCCATTAATATGTGTAATCTTAATAGATCGAATTTCATTTCGTGCATCTTCTCTTATTTGGATTATTGTTCTGTACTGATGGAAGTGTTTAAGGATGTTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATATATGGCATATGCATATTAAACCTATATGC
    SLH_AL_0012                  TTAAGCTGATGGTTAAARCTCTAGGACATGTTATATAATCTATCGTTTCTCGGAGCAGCTATTAACCAAAGTAGTAGTTRATATCGATAGAAAAGTTAGATTGACATAAATCTATGTTGTCTGATTTCGTAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATNGAATTTCATTTCGTGCATCTT---TTATTTGGATTATTGTTCTGTACTGATGGAAGTGTTTAAGGATGTTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATATATGGCATATGCATATTAAACYTATATGC
    SLH_AL_0014                  TTAAGCTGATGGTTAAAGCTCTAGGACATGTTATATAATCTATCGTTTCTCGGAGCAGCTATTAACCAAAGTAGTAGTTAATATCGATAGAAAAGTTAGATTGACANAAATCTATGNTGNCNGNTTTCGNAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATNGAATTTCATTTCGTGCATCTT---TTATTTGGATTATTGTTCTGTACTGATGGAAGTGTTTAAGGATGTTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATATATGGCATATGCATATTAAACCTATATGC
    SLH_AL_0018                  TTAAGCTGATGGTCAAAGCCTTAGGACATGTTATATAATTTATCATTTCTCGGACCAACTAGTAACCAAAGTAGTAGTTAATATCGATAGTAAAGTTAGATGGACATAAATCTATGTTGTCTGATTTCGT--TGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATNGAATTTCATTTCGTGCATCTTCTCTTATTTGGATTATTGTTCTGTACTGATGGAAGTGTTTAAGGATGATAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATATATGGCATATGCATAATAATCCTATATGC
    SLH_AL_0036                  TTAAGCTGATGGTYAAAGCCYTAGGACATGTTATATAATYTATCATTTCTCGGAYCAACTAGTAACCAAAGTAGTAGTTAATATCGATAGTAAAGTTAGATGGACATAAATCTATGTTGTCWGATTTCGT--TGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATNGAATTTCATTTCGTGCATCTT---TTATTTGGATTATTGTTCTGTACTGATGGAAGTGTTTAAGGATGWTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATATATGGCATATGCATAWTAAWCCTATATGC
    SLH_AL_0100                  TTAAGTAGATGGTTAAAGCCCTAGGACATGTTATATAATCTATCATTTCTNGNACCAACTTGTNACCAAAGTAGTAGTTAATATCGATCGTAAAGTTAGATGGACATAAATGTATGTTGTCTGATTTCGTAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATCGAATTTCATTTCGTGCATCTTCTCTTGTTTGGATTATTGTTCTGTACTAATGGAAGTCTTTAAGGATGTTAGCTGCAAAAATGTATTTAAACCTTGAAGTTATTCATGTACGGCATATGCATATTAAACCTATATGC
    SLH_AL_0101                  TTAAGTAGATGGTTAAAGCCCTAGGACATGTTATATAATCTATCATTTCTNGNACCAACTTGTNACCAAAGTAGTAGTTAATATCGATCGTAAAGTTAGATGGACATAAATGTATGTTGTCTGATTTCGTAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATCGAATTTCATTTCGTGCATCTTCTCTTGTTTGGATTATTGTTCTGTACTAATGGAAGTCTTTAAGGATGTTAGCTGCAAAAATGTATTTAAACCTTGAAGTTATTCATGTACGGCATATGCATATTAAACCTATATGC
    SLH_AL_0105                  TTAAGTAGATGGTTAAAGCCCTAGGACATGTTATATGATCTATCATTTCTNGCACCAACTTGTNACCAAAGTAGTAGTTGATATCGATAGTAAAGTTAGATGGACATAAATGTATGTTGTCTGNTTTCGNAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATNGAATTTCATTTCGTGCATCTTCTCTTATTTGGATTATTGTTCTGTACTAATGGAAGTCTTTAAGGATGTTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATGTACGGCATATGCATATTAAACCTATATGC
    SLH_AL_0106                  TTAAGTAGATGGTTAAAGCCCTAGGACATGTTATATAATCTATCATTTCTNGNACCAACTTGTNACCAAAGTAGTAGTTAATATCGATCGTAAAGTTAGATGGACATAAATGTATGTTGTCTGATTTCGTAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATCGAATTTCATTTCGTGCATCTTCTCTTGTTTGGATTATTGTTCTGTACTAATGGAAGTCTTTAAGGATGTTAGCTGCAAAAATGTATTTAAACCTTGAAGTTATTCATGTACGGCATATGCATATTAAACCTATATGC
    SLH_AL_3066                  TTAAGTAGATGGTTAAAGCCCTAGGACATGTTATATGATCTATCATTTCTNGCACCAACTTGTNACCAAAGTAGTAGTTGATATCGATAGTAAAGTTAGATGGACATAAATGTATGTTGTCTGNTTTCGNAGTGTCATTGACAGCAGANTAGAAGCTGTCTCGCCANTGGCGANCTTGATCNGATAAATCNAAATGAGTANACCTGGTTNTAAAACCCANNAATATGTGTAATCTTAATAGATCGAATTTCATTTCGTGCATCTTCTCTTATTTGGATTATTGTTCTGTACTAATGGAAGTCTTTAAGGATGTTAGTTGCAAAAATGTATTTAAACCTTGAAGTTATTCATGTACGGCATATGCATATTAAACCTATATGC
    //                                **      *   - **               *  *    *       * *  *  **                 *        * *          *         *         -                                                                                                                                                    *                     *        *          *   *                                  *  *            *   * -       |14:A_tuberculatus_Chr01:276432-276812
    assembly_reference_sequence  TTAATTGAATATAGTATTAAAAAGTCAACATGATGGGAGGTTTCACTACAAGATTGAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAACTGTTTTCATATGAAGAGGTTATTTAATTATATAACATATCTCATAACCTCAATTATCAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTTGAATAACTATTTCCACAAATCTAAAAAAATAAAATTTTAAAAGGGAAAATATTGAAAAAATGGATTCATTTAGGTTGA
    SLH_AL_0012                  TTAATTGAATATAGTATTAAAAAGTCAACATGACGGGAGGTTTC---ACAAGATTAAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAACTGTTTTCANATGAAGAGGTTATTTAATTATATAACATANCTCATAACNTCAATNATNAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTTGAATAACTATTTCCACAAATCTAAAAAAATAAAATTTTAAAAGGGGAAATATTGAAAAAATGGATTCATTTAGGTTGA
    SLH_AL_0018                  TTAATTGAATATAGTATTAAAAAGTCAACATGACGGGAGGTTTCACTACAAGATTGAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAANTGTTTTCANATGAAGAGGTTATTTAATTATATAACATANCTCATAACNTCAATNATNAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTTGAATAACTATTTCCACAAATCTAAAAAAATAAAATTTTAAAAGGGRAAATATTGAAAAAATGGATTCATTTAGGTTGA
    SLH_AL_0034                  TTAATTGAATATAGTATTAAAAAGTCAACATGACGGKAGGTTTCACTACAAGATTGAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAANTGTTTTCANATGAAGAGGTTATTTAATTATATAACATANCTCATAACNTCAATNATNAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTTGGATAACTATTTCCACAAATYTAAAAAAATAAAATTTTAAAAGGGGAAATATTGAAAAAATGGATTCATTTAGGTTGA
    SLH_AL_0036                  TTAATTGAATATAGTATTAAAAAGTCAACATGACGGGAGGTTTCACTACAAGATTGAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAANTGTTTTCANATGAAGAGGTTATTTAATTATATAACATANCTCATAACNTCAATNATNAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTAGAATAACTATTTCCACAAATCTAAAAAAATAAAATTTTAAAAGGGAAAATATTGAAAAAATGGATTCATTTAGGTTGA
    SLH_AL_0042                  TTAATTGAATATAGTATTAAAAAGTCAACATGATGGGAGGTTTCACTACAAGATTGAATCTCATCGACCTCTCATTTCAATTGAAGTTTATAACGCCAAATTGATGAGAGTCTATGTATATCCACAATTCTAGGCCAAANTGTTTTCANATGAAGAGGTTATTTAATTATATAACATANCTCATAACNTCAATNATNAACTTATTGACTTTATTTTTAGTTATCACTATATATTGGCTTTAAATTTCTTATTTTTCTTCCCATTCCATTGTTGAATAACTATTTCCACAAATCTAAAAAAATAAAATTTTAAAAGGGAAAATATTGAAAAAATGGATTCATTTAGGTTGA
    //                                                            *  -                  -                                                                                                                                                                                                                       - -                  -                        *                                |15:A_tuberculatus_Chr01:280635-280984
```

## Exploratory analyses

### PCA

### Phylogenetic tree inference

Here we use the `ipyrad2 seqex` tool to extract loci from the database file (HDF5)
and write a concatenated sequence alignment of a filtered subset of loci. The argument
`-C` specifies to concatenate (combine end-to-end) the retained loci. You can either
apply this to all loci (i.e., genome-wide), or you can limit the genomic scope by selecting
a subset of scaffolds/chromosomes that loci can be on. This is most relevant when working
with a chromosome-scale reference genome.


hile `-m 10` says
to keep only loci with data from at least 10 samples, and `-r 0.5` says to exclude samples
that have missing data for more than 50% of a locus.

```bash
ipyrad2 seqex \
    -d TUTORIAL/OUT/assembly.hdf5 \
    -P
```
??? "scaffolds report"

    ```literal
    ...
    ```

...
```bash
ipyrad2 seqex \
    -d TUTORIAL/OUT/assembly.hdf5 \
    -o TUTORIAL/OUT/output-seqex/ \
    -n assembly_concat_m10_r0.5 \
    -w Chr01 \
    -C \
    -m 10 \
    -r 0.5
```
??? "scaffolds report"

    ```literal
    ...
    ```

The stats file shows which loci were included in the alignment, the size of the final concatenated alignment, and its
amount of missing sequences (Ns). You can toggle parameters of the `seqex` tool to generate alignments with different
numbers of samples included, and allowing different proportions of missing data to investigate the consistency of your
results to these settings.

This alignment has ...

[insert stats]

We can now use this concatenated sequence alignment in a downstream analysis tool. Here we will use `raxml-ng`, a tool
for inferring a phylogenetic tree
