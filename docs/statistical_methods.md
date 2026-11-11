# Statistical and sequence methods

CDS nucleotides and translated proteins are aligned separately. Longest-CDS representative
transcripts are selected without inspecting variability. Nucleotide alignments are not codon-aware;
protein-guided codon alignments can validate borderline CDS candidates. dN/dS is not estimated.

Callable columns contain A/C/G/T in both sequences (standard amino acids for protein identity).
Substitutions/callable columns yield SNP density and identity. Pairwise gap runs count indel
events; gap bases/(callable+gap bases) measures burden. Combined divergence is
(substitutions+gap bases)/(callable+gap bases), not a model-corrected evolutionary distance.
Callability uses callable columns/max(unambiguous ungapped lengths). Double gaps and ambiguous
bases do not contribute to identity. Missingness and length ranges are retained separately.

Matched valid isolate pairs are used for CDS/upstream comparisons. Independent CDS conservation
also remains available for loci whose upstream sequence fails QC, preserving unresolved denominators.
Upstream identity below 50% is unresolved homology, not confidently biological variation.
Equal-window extraction is anchored to the gene model and does not guarantee distal homology;
large indels require conserved-flank review or graph/local-alignment validation.

R reports A–E counts, conditional prevalence among high-quality coding-conserved loci, and
bounds including unresolved coding-conserved loci. Chromosome-block percentile bootstrap
(2,000 replicates, seed 2026) requires at least five reference chromosomes and conditions on
sampled isolates. No species-wide prevalence or independence of related isolates is asserted.

Exploratory gene-level OLS uses CDS divergence, upstream repeat fraction, log intergenic length,
local repeat density, contig-end distance, missingness and chromosome fixed effects, with
chromosome-clustered sandwich uncertainty. It requires ≥50 loci, ≥8 chromosomes, full-rank
predictors and ≥5 observations per coefficient. Failure of these checks is reported, not bypassed.
BH adjustment covers reported coefficient tests across windows; window analyses remain dependent.
Bounded responses, correlated predictors and small cluster counts limit inference. Pairwise
observations are never independent replicates.

TE distance is retained as a distinct feature, but the current formal regression tests general
repeat overlap/density. A TE-specific model and phylogenetic/relatedness-aware uncertainty are
not claimed implemented. The common isolate set conditions gene-level contrasts; it does not
eliminate lineage composition effects. Repeat associations are descriptive and partly share the
measured sequence substrate. No causal mechanism or positive selection follows from them.
