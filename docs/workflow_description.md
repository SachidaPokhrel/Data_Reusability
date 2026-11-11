# Workflow description and DAG

The main Snakefile defines data_manifest, download, verify_checksums, extract_archives,
inspect_schemas, build_isolate_manifest, validate_annotations, select_single_copy_core_genes,
select_pilot_subset, extract_CDS, extract_upstream, align_CDS, align_upstream,
calculate_divergence, analyze_repeats, validate_synteny, classify_loci, statistical_analysis,
generate_figures, candidate_visualization, generate_report, provenance and audit.

TEST substitutes synthetic fixture generation for acquisition. Graph validation is a separate
optional target. RNA-seq is a proposed future extension, not an implemented dependency.
Rules declare threads, mem_mb, runtime and disk_mb. Stage outputs include real tables/directories
and checksummed stage receipts. No marker alone is scientific evidence. Input code/config are
part of the DAG, so changes trigger reruns; broad code dependency tracking favors safety over
minimal recomputation. CLI checksums provide a second content-based audit.

```bash
# After bootstrap, render the actual engine DAG:
.conda/envs/workflow/bin/snakemake --snakefile Snakefile --configfile config/pilot.yaml \
  --config project_config=config/pilot.yaml --dag > docs/workflow.dag.dot
```

The tracked Mermaid diagram in ../README.md describes dependencies; an actual engine-generated
DAG is not supplied as if it had run here. Missing source schema or ID mappings fail affected
stages. Complete rule environments must be available before offline cluster execution.
