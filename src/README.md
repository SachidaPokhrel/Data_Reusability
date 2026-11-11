# Reusable implementation

pe_upstream contains the executable stage library used by scripts/pe.py and Snakemake.
common handles configuration, atomic tables and provenance; acquisition handles API/downloads;
quality and selection handle mapping/core sampling; sequences handles extraction/alignment;
variation, repeats and validation calculate evidence; graph is optional; reporting produces outputs.

Code is tracked and tested with synthetic fixtures. Run through run.sh or the stage CLI so that
configuration, paths and receipts are retained. These parsers implement explicit contracts;
compatibility with the unpublished-to-this-workspace source schemas remains unverified.
