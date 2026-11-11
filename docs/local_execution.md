# Local execution

Run from the repository root; run.sh also sets its working directory reliably.

```bash
./run.sh test --offline-test
./bootstrap.sh
./run.sh download --cores 2
# Inspect data and populate config/isolate_files.tsv / parser configuration.
./run.sh pilot --cores 8
./run.sh report --config config/pilot.yaml
./run.sh audit --config config/pilot.yaml
./run.sh full --cores 16
./run.sh report --config config/full.yaml
```

For a dry run append --dry-run after installation. Interrupted workflows resume using the
same command and Snakemake --rerun-incomplete. One writer per mode is required. Sources,
config/code changes and checksums invalidate downstream results; run the audit after generating reports.

PROJECT_ROOT overrides the source root. DATA_ROOT, RESULTS_ROOT, REPORTS_ROOT and FIGURES_ROOT
replace storage bases; each results/report/figure base gets a test/pilot/full subdirectory.
CONDA_ROOT sets environment storage; SCRATCH_ROOT is used only by explicit HPC scratch mode.
Use the same settings on reruns. External paths in sample and genotype manifests must remain
accessible. Immutable source data can be staged once for both pilot and full mode.
