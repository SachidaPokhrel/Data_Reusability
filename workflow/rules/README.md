# Rule organization

Rules currently reside together in the root Snakefile so the full dependency structure is
reviewable in one file. This directory is reserved for a future refactor only; it does not
contain placeholder rule implementations. Source stage functions are in src/pe_upstream.
