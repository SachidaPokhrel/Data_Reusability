# Configuration

The .yaml files use JSON syntax, a valid YAML subset that standard-library tooling can also
read. test.yaml confines generated synthetic data/results to tests/work. pilot.yaml selects six
isolates and up to 100 loci; full.yaml uses all 19 and every eligible locus. Seed is 2026.

sources.yaml fixes source identity and requested filenames, not guessed download URLs. API
metadata supplies links/checksums. isolate_files.tsv is an explicit user-reviewed mapping and
intentionally contains only headers until source files are inspected. Set orthology column names
or an ID map only after inspection. Configuration is tracked; secrets and HPC credentials are not.

Storage can be overridden with PROJECT_ROOT, DATA_ROOT, RESULTS_ROOT, REPORTS_ROOT,
FIGURES_ROOT and CONDA_ROOT. Keep overrides consistent across reruns. The CLI validates mode,
thresholds, windows, strict core size and forbids synthetic alignment in biological modes.

The completed biological pilot uses [pilot.zenodo_v7.yaml](pilot.zenodo_v7.yaml) and
[isolate_files.zenodo_v7.tsv](isolate_files.zenodo_v7.tsv). Recreate the gene memberships with
`python3 scripts/qc/prepare_zenodo_v7.py` after the source files are retrieved. Pe6's blank
repeat path is intentional and marked missing; it was not among the six selected isolates.
