# Troubleshooting and resume

## DNS or download failure

If DNS resolution fails, move installation/download to an approved network-enabled machine
or HPC transfer node. Retry bootstrap and the download stage there. Resume files are
retained for Range-capable servers; a server that ignores Range causes a clean restart. Checksum
mismatch deletes the corrupt partial transfer, never renames it to a valid raw file.

## Unknown schemas or identifiers

Read results/pilot/qc/data_inventory.tsv and source headers. Populate the explicit isolate manifest;
configure orthology columns and, if required, complete source-to-gene mapping. Numeric presence
counts cannot map genes. Convert GTF or unsupported repeat formats explicitly. Never fix mapping
by silent name truncation or dropping unmatched records.

## Environment solve or platform failure

Check official package availability on the selected platform. Preserve the failure receipt. A
changed spec requires a new CONDA_ROOT or reviewed environment recreation; do not overwrite an
unrelated environment. Store solved explicit exports after success. No root installation is needed.

## Alignment or full-mode failure

Inspect per-mode logs and unresolved reasons. More divergence is not a reason to relax homology
checks blindly. Full mode needs a real passing pilot gate and unchanged source hashes. Restart
the original command after correcting the cause. Do not run two writers in the same mode.

## HPC failure

Check actual scheduler commands, queue, account, walltime, storage and environment accessibility.
Resubmit explicitly; afterok dependencies do not run downstream jobs after a failed prerequisite.
Scratch outputs are retained on synchronization failure; inspect relocation receipts before cleanup.
No cluster testing is claimed here.
