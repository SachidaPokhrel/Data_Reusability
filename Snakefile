# Snakemake 8.30.0. Configuration is shared with the independently runnable CLI.
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd() / "src"))
from pe_upstream.common import load_config, path, ROOT
CFGFILE = str(Path(config.get("project_config", "config/pilot.yaml")).resolve())
C = load_config(CFGFILE)
MODE = C["mode"]
def io_path(p):
    p = Path(p).resolve()
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
D = io_path(path(C,"data"))
O = io_path(path(C,"results"))
P = io_path(path(C,"reports"))
F = io_path(path(C,"figures"))
CODE = [io_path(p) for p in sorted((ROOT/"src").rglob("*.py"))]
CODE += [io_path(p) for p in sorted((ROOT/"scripts").rglob("*.py"))]
CODE += [io_path(p) for p in sorted((ROOT/"scripts").rglob("*.R"))]
CODE += [io_path(CFGFILE), "config/sources.yaml"]
# Existing source files are explicit dependencies after acquisition/schema mapping.
# The download-only pass intentionally precedes source mapping.
from pe_upstream.common import rows
SOURCE_INPUTS = []
if (ROOT/C['samples']).exists():
    for sample in rows(ROOT/C['samples']):
        for field in ['fasta','gff','repeats']:
            source = Path(sample[field])
            source = source if source.is_absolute() else ROOT/source
            if source.is_file(): SOURCE_INPUTS.append(io_path(source))
for source in [path(C,'data')/'raw'/C['orthology']['file'],
               ROOT/C['distances'] if C.get('distances') else path(C,'data')/'raw/Pe_variant_matrix.phy']:
    if source.is_file(): SOURCE_INPUTS.append(io_path(source))
if C['orthology'].get('id_map'):
    SOURCE_INPUTS.append(io_path(ROOT/C['orthology']['id_map']))
SOURCE_INPUTS = sorted(set(SOURCE_INPUTS))
ANALYSIS = "envs/analysis.yaml"
RENV = "envs/r-analysis.yaml"
WORKFLOW = "envs/workflow.yaml"

def marker(stage):
    return P + "/provenance/" + stage + ".json"

rule all:
    input:
        P + "/" + MODE + "_report.html",
        P + "/provenance/audit.json",
        P + "/provenance/run_manifest.json"

if MODE == "test":
    rule test_fixtures:
        input:
            CODE + []
        output:
            [io_path(ROOT/C["samples"]), D + "/raw/Pe_orthogroups.tsv", marker('test_fixtures')]
        threads: 1
        resources:
            mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
        params: configfile=CFGFILE
        conda: WORKFLOW
        log: "logs/test_fixtures." + MODE + ".log"
        shell:
            "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage test_fixtures > {log:q} 2>&1"

else:
    rule data_manifest:
        input:
            CODE + []
        output:
            [D + "/metadata/zenodo_record.json", D + "/metadata/download_manifest.tsv", marker('data_manifest')]
        threads: 1
        resources:
            mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
        params: configfile=CFGFILE
        conda: WORKFLOW
        log: "logs/data_manifest." + MODE + ".log"
        shell:
            "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage data_manifest > {log:q} 2>&1"

    rule download:
        input:
            CODE + [marker("data_manifest")]
        output:
            [directory(D + "/raw"), D + "/metadata/download_receipts.tsv", marker('download')]
        threads: 1
        resources:
            mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
        params: configfile=CFGFILE
        conda: WORKFLOW
        log: "logs/download." + MODE + ".log"
        shell:
            "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage download > {log:q} 2>&1"

    rule verify_checksums:
        input:
            CODE + [marker("download"), D + "/raw"]
        output:
            [D + "/metadata/verified_files.tsv", marker('verify_checksums')]
        threads: 1
        resources:
            mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
        params: configfile=CFGFILE
        conda: WORKFLOW
        log: "logs/verify_checksums." + MODE + ".log"
        shell:
            "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage verify_checksums > {log:q} 2>&1"

    rule extract_archives:
        input:
            CODE + [marker("verify_checksums")]
        output:
            [directory(D + "/external"), D + "/metadata/archive_inventory.tsv", marker('extract_archives')]
        threads: 1
        resources:
            mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
        params: configfile=CFGFILE
        conda: WORKFLOW
        log: "logs/extract_archives." + MODE + ".log"
        shell:
            "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage extract_archives > {log:q} 2>&1"

rule inspect_schemas:
    input:
        CODE + [marker("test_fixtures") if MODE=="test" else marker("extract_archives")]
    output:
        [O + "/qc/data_inventory.tsv", marker('inspect_schemas')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: WORKFLOW
    log: "logs/inspect_schemas." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage inspect_schemas > {log:q} 2>&1"

rule build_isolate_manifest:
    input:
        CODE + SOURCE_INPUTS + [marker("inspect_schemas"), io_path(ROOT/C["samples"])]
    output:
        [O + "/qc/isolate_manifest.tsv", marker('build_isolate_manifest')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: WORKFLOW
    log: "logs/build_isolate_manifest." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage build_isolate_manifest > {log:q} 2>&1"

rule validate_annotations:
    input:
        CODE + SOURCE_INPUTS + [marker("build_isolate_manifest")]
    output:
        [O + "/qc/annotation_validation.tsv", O + "/qc/assembly_summary.tsv", O + "/qc/qc_report.md", marker('validate_annotations')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/validate_annotations." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage validate_annotations > {log:q} 2>&1"

rule select_single_copy_core_genes:
    input:
        CODE + SOURCE_INPUTS + [marker("validate_annotations")]
    output:
        [O + "/orthogroups/core.json", O + "/orthogroups/single_copy_core_genes.tsv", O + "/qc/orthogroup_summary.tsv", marker('select_single_copy_core_genes')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/select_single_copy_core_genes." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage select_single_copy_core_genes > {log:q} 2>&1"

rule select_pilot_subset:
    input:
        CODE + SOURCE_INPUTS + [marker("select_single_copy_core_genes"), *([io_path(path(dict(C, mode="pilot"),"reports")/"analysis_gate.json")] if MODE=="full" else [])]
    output:
        [O + "/orthogroups/selection.json", O + "/orthogroups/selected_genes.tsv", marker('select_pilot_subset')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/select_pilot_subset." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage select_pilot_subset > {log:q} 2>&1"

rule extract_CDS:
    input:
        CODE + SOURCE_INPUTS + [marker("select_pilot_subset")]
    output:
        [directory(O + "/cds/sequences"), O + "/cds/cds_qc.tsv", marker('extract_CDS')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/extract_CDS." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage extract_CDS > {log:q} 2>&1"

rule extract_upstream:
    input:
        CODE + SOURCE_INPUTS + [marker("select_pilot_subset")]
    output:
        [directory(O + "/upstream/sequences"), O + "/upstream/upstream_coordinates.tsv", O + "/upstream/upstream_qc.tsv", marker('extract_upstream')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/extract_upstream." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage extract_upstream > {log:q} 2>&1"

rule align_CDS:
    input:
        CODE + SOURCE_INPUTS + [marker("extract_CDS"), O + "/cds/sequences"]
    output:
        [directory(O + "/cds/alignments"), O + "/cds/alignment_status.tsv", marker('align_CDS')]
    threads: C["threads"]
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/align_CDS." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage align_CDS > {log:q} 2>&1"

rule align_upstream:
    input:
        CODE + SOURCE_INPUTS + [marker("extract_upstream"), O + "/upstream/sequences"]
    output:
        [directory(O + "/upstream/alignments"), O + "/upstream/alignment_status.tsv", marker('align_upstream')]
    threads: C["threads"]
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/align_upstream." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage align_upstream > {log:q} 2>&1"

rule calculate_divergence:
    input:
        CODE + SOURCE_INPUTS + [marker("align_CDS"), marker("align_upstream"), O + "/cds/alignments", O + "/upstream/alignments"]
    output:
        [directory(O + "/variation/intermediate"), O + "/cds/cds_pairwise.tsv", O + "/upstream/upstream_pairwise.tsv", marker('calculate_divergence')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/calculate_divergence." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage calculate_divergence > {log:q} 2>&1"

rule analyze_repeats:
    input:
        CODE + SOURCE_INPUTS + [marker("extract_upstream")]
    output:
        [O + "/repeats/upstream_repeats.tsv", O + "/repeats/repeat_class_inventory.tsv", marker('analyze_repeats')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/analyze_repeats." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage analyze_repeats > {log:q} 2>&1"

rule validate_synteny:
    input:
        CODE + SOURCE_INPUTS + [marker("align_upstream"), marker("analyze_repeats")]
    output:
        [O + "/synteny/locus_evidence.tsv", O + "/synteny/structural_events.tsv", marker('validate_synteny')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/validate_synteny." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage validate_synteny > {log:q} 2>&1"

rule classify_loci:
    input:
        CODE + SOURCE_INPUTS + [marker("calculate_divergence"), marker("analyze_repeats"), marker("validate_synteny")]
    output:
        [O + "/candidates/master_locus_table.tsv", O + "/statistics/sensitivity_locus.tsv", marker('classify_loci')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/classify_loci." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage classify_loci > {log:q} 2>&1"

rule statistical_analysis:
    input:
        CODE + SOURCE_INPUTS + [marker("classify_loci")]
    output:
        [O + "/statistics/prevalence.tsv", O + "/statistics/classification_counts.tsv", O + "/statistics/statistics_status.json", marker('statistical_analysis')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: RENV
    log: "logs/statistical_analysis." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage statistical_analysis > {log:q} 2>&1"

rule generate_figures:
    input:
        CODE + SOURCE_INPUTS + [marker("statistical_analysis")]
    output:
        [F + "/figure_status.json", marker('generate_figures')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: RENV
    log: "logs/generate_figures." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage generate_figures > {log:q} 2>&1"

rule candidate_visualization:
    input:
        CODE + SOURCE_INPUTS + [marker("classify_loci")]
    output:
        [F + "/candidate_figure_status.json", marker('candidate_visualization')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/candidate_visualization." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage candidate_visualization > {log:q} 2>&1"

rule generate_report:
    input:
        CODE + SOURCE_INPUTS + [marker("statistical_analysis"), marker("generate_figures"), marker("candidate_visualization")]
    output:
        [P + "/" + MODE + "_report.html", P + "/" + MODE + "_report.md", P + "/analysis_gate.json", marker('generate_report')]
    threads: 1
    resources:
        mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: WORKFLOW
    log: "logs/generate_report." + MODE + ".log"
    shell:
        "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage generate_report > {log:q} 2>&1"

rule provenance:
    input: marker("generate_report")
    output: P + "/provenance/run_manifest.json"
    resources: mem_mb=1000, runtime=10, disk_mb=100
    params: configfile=CFGFILE
    conda: WORKFLOW
    log: "logs/provenance." + MODE + ".log"
    shell: "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} provenance > {log:q} 2>&1"

rule audit:
    input: P + "/provenance/run_manifest.json"
    output: P + "/provenance/audit.json"
    resources: mem_mb=1000, runtime=10, disk_mb=100
    params: configfile=CFGFILE
    conda: WORKFLOW
    log: "logs/audit." + MODE + ".log"
    shell: "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} audit > {log:q} 2>&1"

# Optional target only. The main DAG does not require graph validation.
rule graph_validation:
    input: marker("classify_loci"), C["graph"]["file"], C["graph"]["path_mapping"]
    output: O + "/candidates/graph_validation.tsv", marker("graph_validation")
    threads: 1
    resources: mem_mb=C["memory_mb"], runtime=C["runtime_min"], disk_mb=C["temporary_disk_mb"]
    params: configfile=CFGFILE
    conda: ANALYSIS
    log: "logs/graph_validation." + MODE + ".log"
    shell: "PE_THREADS={threads} python3 scripts/pe.py --config {params.configfile:q} stage graph_validation > {log:q} 2>&1"
