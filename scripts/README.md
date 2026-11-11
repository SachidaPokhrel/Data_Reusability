# Executable scripts

pe.py exposes staged acquisition, QC, extraction, comparison, reporting and audit commands.
bootstrap.py manages project-local environments; runner.py dispatches run.sh. Statistical/plotting
R scripts are under statistics/ and visualization/. qc/gtf_to_gff3.py is an explicit source converter.
Reusable Python implementations are in ../src/pe_upstream to avoid duplicated parsers.

Directory READMEs map each task to its real entry point. Scripts/config/docs are tracked; large
outputs and local logs are not. `python3 scripts/pe.py --help` lists actions. Run all commands
from the project root, preferably via run.sh for reproducible environment selection.
