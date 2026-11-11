# Workflow implementation

../Snakefile is the executable main DAG. Stages call ../scripts/pe.py, backed by
../src/pe_upstream modules. Rules declare actual tables/directories and stage receipts;
all modes use the same analysis functions. Config and implementation files are dependencies.
Optional graph validation is separate and RNA-seq follow-up is not implemented.

The source DAG design appears in ../README.md; render the actual Snakemake DAG after bootstrap
using the command in ../docs/workflow_description.md. Engine execution could not run here.
