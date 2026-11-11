# Reproducibility and publication

Each stage records timestamp, command, Git commit, config and seed, and hashes for inputs,
outputs and implementation. The source API record preserves version, checksums and licensing.
Stage receipts live in reports/{mode}/provenance. software_versions.tsv records actual availability
and R statistics additionally saves sessionInfo. Final run_manifest.json and input_checksums.tsv
aggregate the run. An incomplete or failed stage is never labeled complete by an audit.

```bash
./run.sh audit --config config/pilot.yaml
./run.sh validate
```

The audit recomputes hashes and rejects modified files/config/code. Synthetic integration tests
intentionally corrupt an output to verify failure. Deterministic seeds control sampling/bootstrap;
MAFFT/tool/platform effects can differ and are not guaranteed bitwise identical across platforms.
Pinned YAML specs are not solved lockfiles. Exact-build exports are produced only after real
installation and support reuse; never invent an environment lock.

Raw archives, environments, caches, secrets and large alignments are ignored. Track source,
configs, documentation, small QC/summary tables and suitable final reports/figures. Runtime
external storage paths may need review before publishing provenance. No credentials are written.

```bash
git status
git diff --check
git log --oneline
# Only after authentication and inspection:
gh auth login
gh repo create Pe-effusa-conserved-genes-variable-neighborhoods --private --source . --remote origin --push
```

No force-push or remote history replacement is supported. If a remote exists, inspect it first.
The final commit hash and publication status are verified separately from generated scientific
run commits; a later documentation commit does not rewrite a previous run's provenance.

The completed Week 12 pilot and its curated snapshot are described in [the replay guide](reproduce_pilot.md). Original HPC receipts, including incomplete historical audit context, remain local. See the published source and evidence manifests under publication/ for the portable record.
