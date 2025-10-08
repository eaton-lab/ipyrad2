
# TODO
- assemble: add a max-len arg to assemble. It will help to identify when users accidentally
  input WGS bams as rad data, and the assembled beds are waaay too big.
- assemble: replace numba code to be pure numpy and remove numba as a dependency.
- assemble: write_snps: check and cleanup
- map: test mapping rates on edges of short loci.
- map: test mapping exclusion of merged pairs.
- map: handle SE
- assemble: accept a loci.bed file and skip locus delimiting. Allow WGS only data in this case.
- denovo: test on real data.
- denovo: develop SE pipeline.
- trim: collect and summarize trim stats files into one file.
- variants: try bcf is faster than vcf for intermediates.
- variants: check that snps are being removed from indel regions, and that we want that.
- pops: imap for population variant calls in assemble.
- pops: imap for combining technical-reps in assemble [easy for variants, harder for beds?]
- assemble: instead of warning about UMI tags on -u it could check for them.
- assemble: check filters again, why no shared het hitting?
- assemble: organize and write assemble stats
- analysis:
    - test api with logger
    - pca
    - structure
    - bpp
    - raxml-ng
    - treeslider
    - more converters like wex?

# ideas to try/consider
- denovo: WHY NOT use spades/assembly for denovo step instead of vsearch pipeline?
- cli: develop `ipyrad {subcommand}` easy one-liner for full assembly from raw fastqs to outputs.
- cli_main: log only a subset of the CMD because it looks ugly when too many file paths are expanded by the shell.
- assemble: should we include insertions and do a MSA on final loci? Optional? Downsides: much slower,
  and unaligns loci with the VCF.

# ideas explored and abandoned


# Completed
- [x] names: log only the first N delimited names as examples.
- [x] assemble: use same -q -Q for beds and pileup.
- [x] assemble: mask and filter final VCF
- [x] assemble: write both dbs into one h5
- [x] assemble: write snps hdf5
- [x] name: do not accept empty names as success.
- [x] name: if splitting on _ from back fails, try from front.

