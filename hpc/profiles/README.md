# Execution profiles

The implementation uses scheduler-native single-allocation submission templates, not distributed
Snakemake executor plugins. Therefore no misleading per-rule PBS/SLURM profile is installed.
Inside an allocation the same local Snakemake DAG uses configured resources and existing Conda
environments. A site-specific distributed profile must be separately installed and validated if needed.
