# Visualization entry points

Implemented in src/pe_upstream/reporting.py and figures.R, invoked from the root with
`python3 scripts/pe.py --config config/pilot.yaml stage STAGE`.
Supported stages for this category: generate_figures, candidate_visualization. See the main workflow for dependencies.

Inputs come from the immediately upstream stage/config; outputs go to data/, results/{mode},
figures/{mode} or reports/{mode} as appropriate. Source code and this documentation are tracked;
large intermediates and test work are ignored. Reproduce through run.sh or the explicit stage CLI;
this directory does not duplicate the shared Python implementation.
