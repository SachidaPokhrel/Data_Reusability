# Scientific and execution limitations

The biological pilot ran on HPC and completed its report. Full mode was not executed, and
its gate did not pass at either window size. The original fixed-window candidate did not
validate when projected into the larger alignment. This does not resolve the unclassified loci.

Scientific limitations:

- Haploid assembly paths omit within-isolate heterozygosity.
- Core orthogroups may still contain orthology or gene-model errors.
- Coding/synteny QC selects a nonrandom subset of genomic contexts.
- Annotation boundaries are provisional anchors, especially after GTF-derived extents.
- Fixed windows may compare different distal ancestral intervals after an insertion.
- Low-callability or ambiguous regions are unresolved, not gene absences or variable promoters.
- Repeat annotation alone does not demonstrate a TE insertion polymorphism or its direction.
- Class B is a computational screen; automatic evidence is moderate/uncertain until validation.
- Related race isolates are not independent random samples from the species.
- Noncoding divergence can reflect weaker constraint; it does not show selection or expression effects.
- The exploratory repeat model cannot remove all local mutation or architecture confounding.

Optional graph validation requires supplied path mappings and exact full-path/assembly sequence
agreement. Partial walks, missing segment sequences and unsupported coordinates fail. Expression,
motif activity, dN/dS and distributed multi-node scheduling are future extensions, not completed analyses.
