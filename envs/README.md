# Environments

workflow.yaml pins Python, Snakemake and QA tools. analysis.yaml pins sequence/scientific tools.
r-analysis.yaml pins R, plotting/statistics packages and the Python standard-library CLI runtime.
Bootstrap creates these only in project-local CONDA_ROOT, then runs version/import smoke checks.

The linux-64 exact-build exports in locks/ and adjacent *.pin.txt files were generated on
the source HPC. Environments themselves are ignored. Different platforms need their own solve.

For offline transfer, copy a verified installation/cache to shared approved storage and retain
prefix compatibility. Environments and package caches are ignored; specs and appropriate exact
exports are tracked. See ../docs/installation.md for prerequisites and platform limits.
