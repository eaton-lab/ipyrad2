
# PLAN

- update demux (use parse_names)
- update assembler
- imap for populations in assemble.
- imap for technical-reps in assemble [easy for variants, harder for beds]
- write snps hdf5
- mask and filter final VCF
- ipa.analysis
    - test api with logger
    - pca
    - structure
    - bpp
    - raxml-ng
    - treeslider
    - window-extracter
- ipyrad window-extracter -d H5 --scaffolds ... --scaffold-idxs ... --...
- ipyrad convert-seqs -d H5 -o fasta nexus ... --seed 123
- ipyrad convert-snps -d H5 -o treemix structure --seed 123


