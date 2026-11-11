# Stage entry point

Workflow jobs invoke scripts/pe.py with `--config CONFIG stage STAGE_NAME`. Implementations live
in src/pe_upstream; the same entry point is usable for diagnosis and resuming individual stages.
The complete offline integration harness invokes those functions in dependency order, rather than
substituting pseudocode. Snakemake remains the main manager for real runs.
