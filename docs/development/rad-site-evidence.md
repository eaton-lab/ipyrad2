# HDF5 format 4.0: optional evidence across locus intervals

`ipyrad2 assemble --evidence loci` writes HDF5 format 4.0. The default is
lean HDF5 3.0 (`--evidence none`); `--evidence snps` writes HDF5 3.0 with
SNP evidence. These are three values of one option, not separate flags.
No permanent companion BCF is required. Temporary compressed pileups are removed with the assembly tmpdir;
`--keep-tmpdir` retains them for audit.

The scope is **every reference coordinate in the final exported RAD BED**,
including invariant and uncovered bases. Globally excluded loci are outside this
scope. Non-ACGT reference positions are represented but marked ineligible. The
scope is not the entire reference genome and does not remove RAD/locus-selection
ascertainment.

The `loci` mode applies to final retained intervals whether originally delimited
from RAD reads or supplied with `--loci-bed`. It includes SNP evidence as well.
The existing format-4.0 group name `rad_evidence` and stored scope `rad_sites`
remain unchanged; the public option name does not change the storage schema.

Native pileup is captured before `bcftools call` and variant-only filtering,
which can discard invariant-site likelihoods. The usual caller consumes exactly
that captured pileup, so no second read-evidence calculation occurs. Capture
BEDs are retained in the database and every eligible final base must lie within
the captured domains. Missing capture chunks cause failure, not apparent low
coverage. Native sample order, allele cardinality, coordinates and reference
base are validated; duplicate identical allele records are rejected.

The original canonical datasets and SNP evidence group remain unchanged. Root
`version=4.0`, `format_version="4.0"`, and `evidence_scope="rad_sites"` identify
the additional `rad_evidence` group. Its `complete` marker is published only
once the entire stream validates. Failed writes do not publish the group.

| Dataset in `rad_evidence` | Meaning |
|---|---|
| `intervals`, `interval_offsets` | Final half-open RAD intervals and their cumulative base counts |
| `site_map` | One row per base: scaffold index, zero-based position, interval index |
| `reference`, `site_status` | Reference ASCII base; NoRecord, HasRecord, or NonACGT |
| `hard_mask` | Sample × base exclusion indicator, including bases with no record |
| `record_span` | Start and count of native records at each base |
| `record_site`, `QUAL` | Base index and native quality per record |
| `allele_tokens`, `allele_ids`, `allele_offsets` | Dictionary and ragged ordered native REF/ALT tokens |
| `DP`, `GQ`, `status` | Sample × record native depth, quality and evidence status |
| `AD` | Sample × concatenated allele entries, indexed by `allele_offsets` |
| `PL`, `PL_offsets` | Sample × concatenated Number=G vectors and record offsets |
| `samples`, `scaffolds` | Explicit axis labels |
| `native_headers`, `capture_regions` | Original headers and capture BEDs, compressed as bytes |

Numeric evidence uses lossless int32 with -1 for missing native values. Ragged
axes avoid a fixed allele limit and preserve SNP, invariant, symbolic and indel
records, including multiple records anchored at the same position. Allele order
and VCF Number=G ordering are never changed. Literal `<*>` or other symbolic
alleles retain their native meaning; they are not relabelled as a particular
alternative nucleotide. Missing GQ stays missing; no GT, GQ, posterior or
reference-confidence score is invented.

Record statuses are PASS, NoRead, NoLikelihood, HardMask and NonACGT. Raw values
are retained even in excluded or zero-depth cells for exact auditing. The reader
must enforce `usable_PL`, which requires positive depth, complete native PL,
PASS status, and no hard mask. A zero-filled native PL vector at DP=0 is **not**
usable. Bases with no record have a zero-length record span, not a synthetic
reference genotype. Hard masks include sample paralog, overlapping-indel,
consensus-sample and shared restriction-site masks. Low-depth quality thresholds
remain separate analysis choices.

`ipyrad2.analysis.rad_evidence.RadSiteEvidenceReader.read(site_indices, samples)`
returns selected base metadata and all native records at those positions,
including literal allele tuples, meaningful vector lengths and usability.
Selections are bounded; returned vectors are padded with -1 only for convenience.
`GenotypeEvidenceReader` continues reading the unchanged SNP group in v3/v4.
`seqex` and `snpex` retain their existing scope. This release does not add an
all-site BCF exporter or a per-base likelihood-diversity estimator.

Storage should be reported as measured whole-database growth relative to SNP
and lean databases, with component sizes and temporary pileup costs separately.
Pilot intervals are not a full-dataset measurement. The oak benchmark preserves
all 375 samples and covers the full final BEDs on both assembly routes; it
recomputes pileup evidence from existing filtered BAMs without new genotype
calls or assemblies. All-sites evidence remains opt-in while these costs and
appropriate downstream models are evaluated.
