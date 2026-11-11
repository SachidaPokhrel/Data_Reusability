# Reproduce the biological pilot

The completed HPC pilot used 19 source assemblies for gene mapping and quality control, then Pe2, Pe5, Pe14, Pe18, Pe9 and Pe8 for the sequence comparison. The original run was staged interactively on HPC. Its curated results and software versions are supplied; the original machine-specific receipts are not part of Git. A complete historical biological-run audit is not claimed.

## 1. Install and download

Use a supported Linux environment with network access and writable storage. On HPC, install and download on a site-approved login or transfer node. Do not launch the biological analysis there.

```bash
bash bootstrap.sh
bash run.sh test --cores 2
export RESULTS_ROOT="$PWD/.replay/results"
export REPORTS_ROOT="$PWD/.replay/reports"
export FIGURES_ROOT="$PWD/.replay/figures"
bash run.sh download --cores 2
python3 scripts/qc/prepare_zenodo_v7.py
```

The verified source list is [here](../publication/verified_downloads.tsv). The preparation step resolves haplotype transcript IDs through GFF3 parent genes, then rejects cross-group assignments. Expect 19 isolates and 8,200 retained groups. Pe6 has no repeat annotation in this release. It remains in the core set; the original six-isolate pilot does not include Pe6. Do not replace the absent annotation with a zero-coverage file.

## 2. Run the same pilot stages on a compute allocation

Run from the repository root. Repeat the exports if this is a fresh shell. The analysis stages read the downloaded sources under `data/`. They write into the replay directory rather than replacing published tables.

```bash
export RESULTS_ROOT="$PWD/.replay/results"
export REPORTS_ROOT="$PWD/.replay/reports"
export FIGURES_ROOT="$PWD/.replay/figures"

bash scripts/replay_pilot.sh
```

Compare the new prevalence and classification tables against [the original summary](../results/pilot/statistics/prevalence.tsv). Review alignment and environment differences if results change. Bootstrap and synthetic Snakemake execution were observed on the source HPC; the packaged replay wrapper is a convenience around the same stage commands, not a claim of a second real run.

## 3. Validate the hit

The replay wrapper reconstructs the larger Chr5_38 alignment and writes the window projection and between-group identity tables. The original manual review and figure are [documented here](../reports/pilot/Chr5_38_validation.md). These are diagnostic larger intervals, not putative promoters. The automated screen’s class-B value is deliberately preserved in the published master table.

Inspect `.replay/reports/pilot/analysis_gate.json` after a replay. The historical gate did not pass. Do not bypass it or treat a successful report build as a scientific validation.

## 4. Verify the publication snapshot

```bash
sha256sum -c publication/SHA256SUMS
python3 -m unittest discover -s tests -v
```

The snapshot manifest covers the curated evidence files. It is distinct from an end-to-end audit of the biological workflow. Assembly/annotation hashes in `publication/assembly_annotation_sha256.tsv` were recorded on HPC; those files must first be retrieved and extracted to verify them on a new machine.
