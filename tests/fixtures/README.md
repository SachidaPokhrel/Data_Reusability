# Synthetic fixtures

src/pe_upstream/fixtures.py deterministically creates six synthetic assemblies, GFF3 gene/repeat
annotations, orthology memberships and genotype distances. Truth includes multi-exon genes,
negative-strand genes, upstream differences and ambiguous sequence. Unit tests separately construct
known insertion/gap alignments and unsafe/incomplete downloads. Generated files stay in ignored
tests/work; this tracked description documents their purpose, not biological observations.
