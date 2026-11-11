# HPC support: explicit allocation templates

The supported model submits one PBS Pro or SLURM allocation; Snakemake handles dependency and
resource scheduling locally within it. These templates are implemented but not cluster-tested.
No distributed multi-node executor is claimed. Scheduler availability/version, queue and shared
environment accessibility are checked before submission. No institution-based assumptions occur.

```bash
# On an approved internet-enabled node:
./bootstrap.sh
./run.sh download
# Inspect schemas, populate mappings, run and review a real pilot before full mode.
./run.sh pilot --executor slurm --cores 8 --account YOUR_ACCOUNT --queue YOUR_PARTITION
./run.sh full --executor slurm --cores 16 --account YOUR_ACCOUNT --queue YOUR_PARTITION
./run.sh full --executor pbs --cores 16 --account YOUR_ACCOUNT --queue YOUR_QUEUE
```

Append --dry-run to prepare without submission (scheduler/queue validation still applies).
Pass --dependency JOB_ID for native afterok dependencies. HPC_ACCOUNT/HPC_QUEUE may replace
CLI settings; HPC_MEMORY_MB/HPC_RUNTIME_MIN replace resource defaults. No email notification is
configured. Inspect scheduler policy before choosing values.

--scratch requires site-approved SCRATCH_ROOT and rsync. Inputs remain on shared storage;
mode-specific output directories are staged with a lock, synchronized on success or failure,
and provenance paths relocated only after checksum verification. Original receipts are retained.
Never clear a scratch lock without verifying no active writer exists. Retry the same command to
resume incomplete jobs; no automatic resubmission loop or unrequested cluster submission occurs.

Scheduler logs, generated job scripts and /usr/bin/time resource diagnostics are ignored under
logs/hpc. Snakemake tracks per-rule mem_mb/runtime/disk_mb/threads within the allocation.
