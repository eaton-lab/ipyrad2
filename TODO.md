# WORKING
- [x] refine bedgraph method
- [x] get min-cov masks from bedgraph
- [] support individual-level variant calling by default; optional group called variants
- [] consensus writing exclude if >X% variants are lowQual Ns.
- [] should soft clip counting be used to filter mapped reads before variant calling, or also for paralog detection?
- [] should/can we get/keep phased results?

# TODO
- trim: collect and summarize trim stats files into one file.
- denovo: test on real data.
- map: add checks in concat_tech_reps_into_tmpdir for bad names.
- map: collect stats from stats files into write_stats
- assemble: need to bring back a max-H consensus filter.
- assemble: max_snps filter is not working. Consider variant quality with this...?
- assemble: allow imap to group samples for joint variant calls.
- assemble: consensus do not write empty seqs, use ordered iterators by chr:pos
- assemble: accept a loci.bed file and skip locus delimiting. Allow WGS only data in this case.
- assemble: instead of warning about UMI tags on -U it could check for them.
- assemble: check filters again, why so few shared het hitting?
- assemble: imap for combining technical-reps in assemble [easy for variants, harder for beds?]
- assemble: organize and write assemble stats
- variants: check that snps are being removed from indel regions, and that we want that (in func: get_vcf_with_indels_resolved)
- variants: try bcf is faster than vcf for intermediates.
- parallel: stderr writes temporarily to $TMPDIR. Instead set to outdir/tmpdir?
- wex: check that alignments are correct when using -r or other exclude methods.
- wex: better format for stats. Write missing per sample to stats.
- analysis:
    - [ ] test api with logger
    - [x] pca
    - [ ] structure
    - [ ] bpp - Done but not tested
    - [ ] raxml-ng
    - [ ] treeslider
    - [ ] more converters like wex?


# TODO Low priority
- demux: update to parallelize using run_with_pool like all other code.
- trim: expose more fastp parameters to user (maybe just a generic --kwargs option)
- trim: allow users to enter their own (additional) adapters.
- variants: try bcf is faster than vcf for intermediates.


# ideas to try/consider for V2
```bash
# write 500 loci randomly sampled from chroms 1-5 using random seed 123 and requiring loci to contain 10 samples
ipyrad2 lex -d DATA.hdf5 -s Chr[1-5] -n 500 -m 10 -x 123 --format fasta
# write 1000 loci randomly sampled from chr1 or chr2 requiring minlen, and that loci have MINMAP cov across IMAP population
ipyrad2 lex -d DATA.hdf5 -s Chr1 Chr2 -n 1000 --min-len 100 --imap IMAP.tsv --minmap MINMAP.tsv --format fasta
# write 500 loci from chroms 1-5 requiring data to be across at least 4 individuals in each pop.
ipyrad2 lex -d DATA.hdf5 -s Chr[1-5] -n 500 --imap IMAP.tsv --minmap 4 --format bpp
...


ipyrad2 analysis wex -d HDF5 ... -O phy
ipyrad2 analysis lex -d HDF5 --imap ... --minmap ... -O bpp
ipyrad2 analysis tex -d HDF5 ... --scaffolds chr1 -c 20 -t 10 --wsize 1e5 --ssize 1e5 --binary raxml-ng --kwargs 'model=GTR+G,boots=100,trees=pars{10}'
ipyrad2 analysis tex -d HDF5 ... --scaffolds chr1 -c 20 -t 10 --wsize 1e5 --ssize 1e5 --binary veryfasttree --kwargs '-gtr,-gamma'
ipyrad2 analysis pca -d HDF5 ... --scaffolds [defaults to all] -n name-prefix -o dir --imap ... --minmap ... --nreplicates 50 -c 10

```
- assemble: option to keep genotype likelihoods, instead of discarding (easy). Should be easy, just don't toss them.
- assemble: use -q and -Q for alignment and base qualities? Then what about other two qualities? Besides -q/-Q already has different meaning in trim.
- assemble: add option --min-loci-per-sample (or similar?) to exclude samples from outfiles if < nloci
- cli_main: log only a subset of the CMD because it looks ugly when too many file paths are expanded by the shell. It could look like this: `CMD: ipyrad assemble -d map_hq2/anas-DE110-plate_8.filtered.bam ...[X sample paths; not shown] -r ... `
- assemble: should we include insertions and do a MSA on final loci? Optional? Downsides: much slower,
  and unaligns loci with the VCF positions, so relabel coordinates needed.
- assemble: add a max-len arg to assemble. It will help to identify when users accidentally
  input WGS bams as rad data, and the assembled beds are waaay too big.
- denovo: WHY NOT use spades/assembly for denovo step instead of vsearch pipeline?
  - ideas: assemble all samples data together, or per-sample and then find consensus?
  - if former: does spades do poorly when samples are highly divergent?
  - if former: can we implement graph splitting similar to current vsearch pipeline?
  - if latter: the need to align and get consensus makes this not so different from current vsearch pipeline.
- assemble: how to handle if many samples have >0 coverage of a locus, but none of its variant sites pass filtering? Currently it gets included but looks invariant. [Discussed that it is not worth addressing right now, can recommmend using better parameters].
- [] write bi-allelic SNPs as binary to nexus for MLE_BiMarkers in phylonet.


# ideas explored and abandoned (maybe don't try these again)
- [x] map: if imap is supplied then *only* process the samples from -d with matching names in imap. Print one warning about others being ignored.


# Completed
- [x] trim: increase min q, add poly-x trimming.
- [x] trim: collect and summarize trim stats files into one file.
- [x] map: skip existing bams, warn about -f, and process any unfinished inputs.
- [x] map: collect stats from stats files into write_stats
- [x] map: allow combining fastqs using --imap
- [x] add `ipyrad-classic` mode.
- [x] assemble: organize and write assemble stats with depths.
- [x] assemble: increase the hdf5 chunk size in Mb.
- [x] assemble: write_snps: check and cleanup
- [x] names: log only the first N delimited names as examples.
- [x] assemble: use same -q -Q for beds and pileup.
- [x] assemble: mask and filter final VCF
- [x] assemble: write both dbs into one h5
- [x] assemble: write snps hdf5
- [x] name: do not accept empty names as success.
- [x] name: if splitting on _ from back fails, try from front.
- [x] denovo: test on real data.
- [x] denovo: develop SE pipeline.
- [x] analysis: add analysis tools (wex/lex/treeslider) as 'analysis' subcommand
- [x] map: handle SE
- [x] assemble: replace numba code to be pure numpy and remove numba as a dependency.
- [x] cli: subcommand splash messages should show the ipyrad version.
