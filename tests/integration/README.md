# Integration tests

The offline integration invokes the complete synthetic stage harness, repeats it, verifies stable
master tables and tests audit rejection after tampering. A real MAFFT test is skipped only when
the executable is unavailable. Snakemake/R scientific integration requires installed environments;
CI provides a separate job for it. This environment has only executed the offline harness.
