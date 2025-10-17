
# TODO
- wex: check that alignments are correct when using -r or other exclude methods.
- assemble: how to handle if many samples have >0 coverage of a locus, but none of its variant sites pass filtering? Currently it gets included but looks invariant.
- assemble: add option --min-loci-per-sample (or similar?) to exclude samples from outfiles if < nloci
- parallel: stderr writes temporarily to $TMPDIR. Instead set to outdir/tmpdir?
- trim: would be nice/convenient to have rename/combine option in trim similar to the one in map.
- trim: maybe a simple option along with rename is to add a prefix or suffix to out file names.
- trim: allow users to enter their own (additional) adapters.
- demux: update to parallelize using run_with_pool like all other code.
- assemble: replace numba code to be pure numpy and remove numba as a dependency.
- assemble: write_snps: check and cleanup
- assemble: consensus do not write empty seqs, use ordered iterators by chr:pos
- map: collect stats from stats files into write_stats
- map: handle SE
- assemble: accept a loci.bed file and skip locus delimiting. Allow WGS only data in this case.
- denovo: test on real data.
- denovo: develop SE pipeline.
- trim: collect and summarize trim stats files into one file.
- variants: try bcf is faster than vcf for intermediates.
- variants: check that snps are being removed from indel regions, and that we want that.
- pops: imap for population variant calls in assemble.
- pops: imap for combining technical-reps in assemble [easy for variants, harder for beds?]
- wex: better format for stats. Write missing per sample to stats.
- assemble: instead of warning about UMI tags on -u it could check for them.
- assemble: check filters again, why no shared het hitting?
- assemble: organize and write assemble stats
- add analysis tools (wex/lex/treeslider) as 'analysis' subcommand
- analysis:
    - test api with logger
    - pca
    - structure
    - bpp
    - raxml-ng
    - treeslider
    - more converters like wex?

# ideas to try/consider
```bash
# write 500 loci randomly sampled from chroms 1-5 using random seed 123 and requiring loci to contain 10 samples
ipyrad lex -d DATA.hdf5 -s Chr[1-5] -n 500 -m 10 -x 123 --format fasta
# write 1000 loci randomly sampled from chr1 or chr2 requiring minlen, and that loci have MINMAP cov across IMAP population
ipyrad lex -d DATA.hdf5 -s Chr1 Chr2 -n 1000 --min-len 100 --imap IMAP.tsv --minmap MINMAP.tsv --format fasta
# write 500 loci from chroms 1-5 requiring data to be across at least 4 individuals in each pop.
ipyrad lex -d DATA.hdf5 -s Chr[1-5] -n 500 --imap IMAP.tsv --minmap 4 --format bpp
...
```
- checkpointing: don't bail out if one result exists, skip finished ones and do unfinished ones?
- option to keep genotype likelihoods, instead of discarding (easy).
- denovo: WHY NOT use spades/assembly for denovo step instead of vsearch pipeline?
  - ideas: assemble all samples data together, or per-sample and then find consensus?
  - if former: does spades do poorly when samples are highly divergent?
  - if former: can we implement graph splitting similar to current vsearch pipeline?
  - if latter: the need to align and get consensus makes this not so different from current vsearch pipeline.
- cli_main: log only a subset of the CMD because it looks ugly when too many file paths are expanded by the shell.
- assemble: should we include insertions and do a MSA on final loci? Optional? Downsides: much slower,
  and unaligns loci with the VCF positions, so relabel coordinates needed.
- assemble: add a max-len arg to assemble. It will help to identify when users accidentally
  input WGS bams as rad data, and the assembled beds are waaay too big.

# ideas explored and abandoned
...


# Completed
- [x] add `ipyrad-classic` mode.
- [x] assemble: increase the hdf5 chunk size in Mb.
- [x] map: allow combining fastqs using --imap
- [x] names: log only the first N delimited names as examples.
- [x] assemble: use same -q -Q for beds and pileup.
- [x] assemble: mask and filter final VCF
- [x] assemble: write both dbs into one h5
- [x] assemble: write snps hdf5
- [x] name: do not accept empty names as success.
- [x] name: if splitting on _ from back fails, try from front.

