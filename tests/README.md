# Tests

Run ./run.sh test --offline-test or `python3 -m unittest discover -s tests -v`. pytest can run the
same tests when installed. Unit tests cover extraction, coordinate conversion, CDS phase,
strand, repeat intervals, indel evidence, coverage, orthology and transfer integrity. Integration
runs all synthetic stage functions twice, checks reproducibility, then verifies audit failure on
output corruption. A real MAFFT test runs only when MAFFT is installed.

No external biological data are required. Synthetic work/results/reports go under ignored
 tests/work. Generated plots are never written into biological result folders. The optional
scientific CI executes Snakemake/R integration after environment installation; that integration
has not run in the current environment.
