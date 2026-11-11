# Sourced inside an explicitly submitted allocation. Inputs remain on shared storage.
: "${SCRATCH_ROOT:?Set SCRATCH_ROOT to site-approved temporary storage}"
: "${PE_RUN_MODE:?Set PE_RUN_MODE}"
PE_SHARED_RESULTS="${RESULTS_ROOT:-$PROJECT_ROOT/results}"
PE_SHARED_REPORTS="${REPORTS_ROOT:-$PROJECT_ROOT/reports}"
PE_SHARED_FIGURES="${FIGURES_ROOT:-$PROJECT_ROOT/figures}"
PE_SCRATCH_RUN="$SCRATCH_ROOT/pe-effusa-$PE_RUN_MODE"
mkdir -p "$PE_SCRATCH_RUN"
if ! mkdir "$PE_SCRATCH_RUN/.lock"; then
  echo "Scratch analysis already locked: $PE_SCRATCH_RUN. Check for active jobs before removing a stale lock." >&2
  exit 2
fi
export RESULTS_ROOT="$PE_SCRATCH_RUN/results"
export REPORTS_ROOT="$PE_SCRATCH_RUN/reports"
export FIGURES_ROOT="$PE_SCRATCH_RUN/figures"
mkdir -p "$RESULTS_ROOT" "$REPORTS_ROOT" "$FIGURES_ROOT"
# Resume prior outputs and retain the successful real pilot gate for full mode.
for PE_MODE_COPY in pilot "$PE_RUN_MODE"; do
  if [[ -d "$PE_SHARED_REPORTS/$PE_MODE_COPY" ]]; then
    rsync -a "$PE_SHARED_REPORTS/$PE_MODE_COPY/" "$REPORTS_ROOT/$PE_MODE_COPY/"
  fi
done
if [[ -d "$PE_SHARED_RESULTS/$PE_RUN_MODE" ]]; then
  rsync -a "$PE_SHARED_RESULTS/$PE_RUN_MODE/" "$RESULTS_ROOT/$PE_RUN_MODE/"
fi
pe_sync_outputs() {
  PE_JOB_EXIT=$?
  trap - EXIT
  mkdir -p "$PE_SHARED_RESULTS/$PE_RUN_MODE" "$PE_SHARED_REPORTS/$PE_RUN_MODE" "$PE_SHARED_FIGURES/$PE_RUN_MODE"
  PE_SYNC_FAILED=0
  rsync -a "$RESULTS_ROOT/$PE_RUN_MODE/" "$PE_SHARED_RESULTS/$PE_RUN_MODE/" || PE_SYNC_FAILED=1
  rsync -a "$REPORTS_ROOT/$PE_RUN_MODE/" "$PE_SHARED_REPORTS/$PE_RUN_MODE/" || PE_SYNC_FAILED=1
  rsync -a "$FIGURES_ROOT/$PE_RUN_MODE/" "$PE_SHARED_FIGURES/$PE_RUN_MODE/" || PE_SYNC_FAILED=1
  if [[ "$PE_SYNC_FAILED" -eq 0 ]]; then
    python3 hpc/rebase_provenance.py --reports "$PE_SHARED_REPORTS" --mode "$PE_RUN_MODE" \
      --mapping "$RESULTS_ROOT" "$PE_SHARED_RESULTS" \
      --mapping "$REPORTS_ROOT" "$PE_SHARED_REPORTS" \
      --mapping "$FIGURES_ROOT" "$PE_SHARED_FIGURES" || PE_SYNC_FAILED=1
  fi
  rmdir "$PE_SCRATCH_RUN/.lock"
  if [[ "$PE_SYNC_FAILED" -ne 0 ]]; then
    echo "Output synchronization failed; retained scratch: $PE_SCRATCH_RUN" >&2
    exit 3
  fi
  exit "$PE_JOB_EXIT"
}
trap pe_sync_outputs EXIT
