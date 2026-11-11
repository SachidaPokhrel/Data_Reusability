# Data dictionary and parser contracts

These are implemented contracts, **not claims about the unseen deposited schemas**.

## Metadata and input mapping

The API response is preserved unchanged. download_manifest.tsv contains source, record_id,
version, file_name, file_size, checksum (algorithm:digest), download_url, download_date,
local_path and validation_status. Receipt and verification tables retain transfer outcomes.
archive_inventory.tsv lists members and compressed/uncompressed lengths before safe extraction.

config/isolate_files.tsv columns: isolate, haplotype, assembly_role, fasta, gff, repeats.
All 19 isolates are required. Representations cannot be mixed silently. Supported roles are
haploid_primary, haploid_consensus, phased_haplotype and synthetic (TEST only). The haplotype
label must be source-supported, including an explicit unphased/consensus label when applicable.
Paths may be absolute or relative to PROJECT_ROOT.

## Sequence and annotation contracts

FASTA first tokens are unique contig IDs. GFF3 is nine-column, one-based inclusive, with
unique gene IDs, transcript Parent links and phased CDS segments. Internal/output intervals
are zero-based half-open. Duplicate gene/transcript IDs and contig/coordinate mismatches fail.
The longest CDS transcript is selected deterministically, with transcript ID breaking ties.
Initial partial CDS, inconsistent phase and internal standard-code stops are flagged/rejected.
CDS-specific alternative genetic codes are not supported.

GTF must first be converted with scripts/qc/gtf_to_gff3.py. Source gene boundaries are retained;
otherwise boundary_origin is marked derived_from_GTF_feature_extents and needs review.
Unknown repeat formats need explicit conversion to GFF3. Broad container features must be
removed before calculating repeat coverage if the source includes them.

## Orthology

Normalized table: orthogroup, isolate, gene_id, one gene membership per row. Long input can
be configured with explicit column names. Automatic recognition only accepts these exact
long columns or a wide table whose non-group headers exactly equal manifest isolate IDs.
Numeric copy counts cannot be interpreted as gene IDs. Explicit mapping tables contain
isolate, source_id, gene_id and must be complete/unique. All memberships are checked against GFF3.

## Output metrics

CDS/upstream pairwise tables contain alignment_length, callable_sites, callable_fraction,
identity, snp_density_per_kb, indel_events, indel_bases, gap_fraction, large_indel_events,
combined_divergence, length_difference and unaligned_fraction. Identity uses only callable
aligned bases; callable fraction uses the larger unambiguous ungapped length. Indels count
contiguous pairwise gap runs and affected bases separately; ancestral direction is unknown.

Upstream coordinates retain requested/actual length, strand, annotation anchor, genomic bounds,
upstream-gene/intergenic distance, overlap, ambiguity, gap-run and contig-edge flags. Clipped
sequences remain recorded; failed comparability is unresolved, not inferred variation.

Master locus rows contain independent-CDS and matched-pair measures, protein identity, all-19
status, repeat/TE features, synteny status, assembly flags, candidate class and unresolved reasons.
Missing numeric fields are empty cells (R reads NA); they are never replaced with biological zero.
Repeat class is unknown unless supported by Classification/Class attributes or an explicit
Name suffix after '#'. TE-only features use configured recognized class roots.

upstream_length_variability is variability in the actual extracted window lengths, primarily a QC
quantity under fixed-window extraction. intergenic_length_variability separately records differences
in the annotation-derived upstream gene distance. Neither is automatically an insertion event.
