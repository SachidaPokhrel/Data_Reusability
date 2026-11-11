#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export RESULTS_ROOT="${RESULTS_ROOT:-$PWD/.replay/results}"
export REPORTS_ROOT="${REPORTS_ROOT:-$PWD/.replay/reports}"
export FIGURES_ROOT="${FIGURES_ROOT:-$PWD/.replay/figures}"
# Prevent the convenience wrapper from replacing the published evidence.
python3 - <<'PY'
import os
from pathlib import Path
for env, original in [('RESULTS_ROOT','results'), ('REPORTS_ROOT','reports'), ('FIGURES_ROOT','figures')]:
    if (Path(os.environ[env]) / 'pilot').resolve() == Path(original, 'pilot').resolve():
        raise SystemExit(f'{env} must point to a replay location')
PY
cfg=config/pilot.zenodo_v7.yaml
conda_root="${CONDA_ROOT:-$PWD/.conda}"
analysis="$conda_root/envs/analysis"
r_analysis="$conda_root/envs/r-analysis"
workflow="$conda_root/envs/workflow"
if [[ ! -x "$analysis/bin/python" || ! -x "$workflow/bin/conda" ]]; then
    echo 'Run bootstrap and retrieve/prepare the inputs before the replay.' >&2
    exit 2
fi
for stage in build_isolate_manifest validate_annotations select_single_copy_core_genes select_pilot_subset extract_CDS extract_upstream align_CDS align_upstream calculate_divergence analyze_repeats validate_synteny classify_loci; do
    PATH="$analysis/bin:$PATH" "$analysis/bin/python" scripts/pe.py --config "$cfg" stage "$stage"
done
for stage in statistical_analysis generate_figures; do
    "$workflow/bin/conda" run -p "$r_analysis" python scripts/pe.py --config "$cfg" stage "$stage"
done
"$analysis/bin/python" scripts/pe.py --config "$cfg" stage candidate_visualization
"$workflow/bin/python" scripts/pe.py --config "$cfg" stage generate_report
"$workflow/bin/python" scripts/pe.py --config "$cfg" provenance
"$analysis/bin/python" scripts/qc/review_chr5_38.py --results "$RESULTS_ROOT/pilot" \
    --rebuild --mafft "$analysis/bin/mafft"
