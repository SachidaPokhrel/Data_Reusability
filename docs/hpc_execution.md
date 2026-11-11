# HPC execution

See [HPC templates](../hpc/README.md). Scheduler detection is based on commands/version and
queue inspection, never institution name. This environment has neither PBS Pro nor SLURM;
no real cluster execution has been tested or submitted.

The implementation submits one allocation and uses Snakemake's local executor inside it.
This provides per-rule dependency/resource scheduling within an explicit allocation without
claiming distributed per-rule execution. For multi-node distribution, adopt and validate a
site-supported Snakemake executor plugin separately.

Use approved internet-enabled nodes for installation and transfer. Compute jobs use prebuilt
shared environments and staged inputs. HPC account/queue must be provided explicitly. Walltime
and memory use configurable defaults or HPC_RUNTIME_MIN/HPC_MEMORY_MB. No email is enabled.
Native afterok dependencies and explicit resubmission support workflow restart. Logs include
scheduler output and /usr/bin/time resource measurements where available.
