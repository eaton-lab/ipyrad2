# Assembly statistics schema 2.0

`NAME.stats.json` is the machine-readable report; `NAME.stats.txt` renders the
same statistics for review. Counts are integers, fractions/times are numbers,
and unavailable values are null. Consumers should check `run.status` before
using biological summaries. A failure report intentionally lacks those summaries.

| Key | Contents |
| --- | --- |
| `stats_schema_version` | `"2.0"`; independent of HDF5 versions 3.0/4.0 |
| `run` | Unique run ID, status, UTC start/end, elapsed seconds, package version, working directory |
| `run.invocation` | CLI `argv`, shell-quoted `command`, canonical `display_argv`; Python API runs have null CLI tokens |
| `run.parameters` | Resolved function parameters, including defaults; paths expanded to absolute paths |
| `run.inputs` | Final sample names mapped to selected original RAD/WGS BAM paths |
| `run.calling` | Independent/grouped mode and normalized sample/population assignments |
| `run.configuration_files` | Supplied populations, rename, subsample and BED paths and SHA256 hashes |
| `command` | Full reconstructed CLI command; null for API calls unless a descriptive command was explicitly supplied |
| `summary`, `alignment_summary` | Locus/site/SNP counts, occupancy and coverage summaries |
| `locus_filtering`, `sample_masking` | Filter and mask counts |
| `sample_summary`, `locus_occupancy` | Complete sample statistics and occupancy tables |
| `read_filter_diagnostics` | Alignment-record failure accounting by sample; independent failures can overlap |
| `mixed_rad_wgs_diagnostics` | Present when mixed RAD/WGS diagnostics were computed |
| `definitions` | Units, numerator/denominator definitions, depth aggregation and evidence semantics |
| `evidence` | Mode, storage status, SNP/all-locus summaries, writing times, temporary BCF bytes/disposition |
| `storage` | Whole HDF5 file bytes and physical evidence dataset allocation |
| `hdf5_format_version` | Actual database format, separate from stats schema |

Each non-null `evidence.snps` / `evidence.loci` summary contains `samples`, `cells`,
`status_counts`, `usable_pl_cells`, `complete_field_cells`, `per_sample`, and
`by_sample_type`. Every per-sample row has its own cell denominator. SNP summaries
also record `sites` and `matched_native_records`. Locus summaries additionally
record `records`, `site_status_counts`, and `symbolic_allele_records`. Status
counts partition the cell denominator; field availability counts overlap them.
SNP array padding is excluded when checking complete vectors.

Usable SNP PL requires `PASS` or `NoHardCall` (positive DP and complete PL).
Usable locus PL requires `PASS`. Locus raw fields retained behind hard masks
can contribute to complete-field counts, but not usable PL. A usable native
likelihood need not be suitable for every downstream model: symbolic alleles,
indels and locus ascertainment still require analysis-specific treatment.

Evidence writers accumulate these counters during streaming; report creation
does not reread large genotype arrays. Each evidence group stores its compact
writer summary in a compressed `summary_json` byte dataset. The database's
`assembly_report/summary_json` stores the stats schema, run and evidence metadata
as UTF-8 JSON bytes, excluding filesystem storage counts. Standard genotype,
sequence, SNP-map and evidence arrays retain their original values and ordering.

`run.invocation.argv` is the lossless argument representation, including spaces
or control characters inside arguments. Shell reconstruction uses POSIX quoting.
The bounded text command escapes control characters for display; do not use it
as provenance or as a script. CLI tokens cannot recover unexpanded shell globs,
redirections or shell variables. Configuration hashes identify small control
files; BAMs and the reference FASTA are not rehashed for this report.
