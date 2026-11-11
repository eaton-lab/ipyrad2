# HDF5 format 3.0

Lean and SNP-evidence assembly databases use root `version=3.0` and `format_version="3.0"`.
This is a storage format version, independent of the ipyrad2 package version.
Root `evidence_scope` is `none` by default (`assemble --evidence none`), or
`snps` with `assemble --evidence snps`. Both retain the ordered SNP allele
catalogue; only `snps` writes the `evidence` group. `--evidence loci` adds
the [all-locus group](rad-site-evidence.md) and uses format 4.0.
No external file links are required. The public modes replace the former
`--no-genotype-evidence` and `--all-site-evidence` flags.

The existing `phy`, `phymap`, `genos`, `snpsmap`, `reference`, `sample_dp` and
`site_qual` arrays retain their meanings. Sequence extraction does not read
the evidence arrays. Native DP is stored separately from `sample_dp` because
canonical masking may have changed or removed that value.

| Dataset | Shape | Representation |
|---|---|---|
| `alleles` | SNP × 4 | uint8 ASCII REF followed by ordered ALT; zero padding |
| `evidence/DP` | sample × SNP | int32 native depth |
| `evidence/AD` | sample × SNP × 4 | int32 native allele depths |
| `evidence/GQ` | sample × SNP | int32 native caller quality |
| `evidence/PL` | sample × SNP × 10 | int32 native diploid genotype likelihoods |
| `evidence/status` | sample × SNP | uint8 status code |

Sample order is `genos.attrs['names']`; SNP order is `snpsmap`. Diploid SNPs
have at most four A/C/G/T alleles, so padded AD/PL slots cover multiallelic
records without variable-length objects. For K alleles only the first K AD
entries and K(K+1)/2 PL entries are meaningful. PL uses VCF Number=G order:
00, 01, 11, 02, 12, 22, 03, 13, 23, 33. Integers are preserved without clipping
or quantization. Missing values and unused slots are -1; allele count
distinguishes padding. Evidence arrays use gzip level 4, shuffle, and bounded
SNP chunks.

| Code | Meaning |
|---|---|
| 0 PASS | Complete native PL, positive depth, nonmissing final GT |
| 1 NoHardCall | Complete native PL and positive depth, but final GT missing |
| 2 NoRead | Native DP zero or unavailable; AD/GQ/PL missing |
| 3 HardMask | Explicit assembly exclusion; all evidence missing |
| 4 NoMatch | No unique exact native position/ordered-allele match |
| 5 NoLikelihood | Positive depth, but native PL absent or incomplete |

Covered native AD/GQ remain available for NoLikelihood cells. Partial fields
retain explicit missing entries. A zero PL is not a missing-value sentinel.
Genotype posterior probabilities are not stored or inferred by this writer.
Hard masks include paralog, overlapping-indel, consensus sample-row and
restriction-site masks. NoRead and HardMask cells cannot become usable PL.

`evidence/native_header`, `evidence/final_header`, provenance, record-match
counts and status counts document the join and original caller metadata.
`evidence.attrs['complete']` is published only after every HDF5 SNP has been
reconciled against the canonical VCF. An unmatchable native SNP remains in
the standard database with NoMatch evidence rather than being silently dropped.

`ipyrad2.analysis.evidence.GenotypeEvidenceReader` provides bounded selected
reads, allele counts, final GT and a `usable_PL` mask. Consumers must retain
the mask and distinguish likelihoods from posteriors. Absent records/invariant likelihoods cannot be reconstructed from the SNP
group. The opt-in [HDF5 v4 RAD-site group](rad-site-evidence.md) captures native
evidence upstream; its symbolic alleles need an appropriate separate consumer.

`snpex -O bcf` writes filtered calls plus native evidence. It does not infer a
new genotype, recompute GQ, infer phase, or reorder the native allele catalogue.
Additional snpex call masking is labeled ExportMask. The usual biallelic SNP
filter still applies to snpex; the HDF5 reader itself supports multiallelic SNPs.
Imputation and synthetic reference samples are rejected for evidence export.
