"""Command-line entry point for restartable, independently logged stage operations."""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from .common import ROOT,load_config,path,record_stage,read_json,write_json,rows,table,audit,now,software,digest,rel,git_hash
from . import acquisition,quality,selection,sequences,variation,repeats,validation,reporting,fixtures,graph

ACTIONS={
 'data_manifest':acquisition.manifest,'download':acquisition.download,
 'verify_checksums':acquisition.verify,'extract_archives':acquisition.extract,
 'inspect_schemas':quality.inspect,'build_isolate_manifest':quality.build_manifest,
 'validate_annotations':quality.validate,'select_single_copy_core_genes':quality.core_groups,
 'select_pilot_subset':selection.select,'extract_CDS':sequences.cds,'extract_upstream':sequences.upstream,
 'align_CDS':lambda c:sequences.align(c,'cds'),'align_upstream':lambda c:sequences.align(c,'upstream'),
 'calculate_divergence':variation.calculate,'analyze_repeats':repeats.analyze,
 'validate_synteny':validation.validate,'classify_loci':variation.classify,
 'statistical_analysis':reporting.statistics,'generate_figures':reporting.figures,
 'candidate_visualization':reporting.extra_figures,'generate_report':reporting.report,
 'test_fixtures':fixtures.create,'graph_validation':graph.validate}


def stage_inputs(c,stage):
    """Only actual upstream dependencies; old downstream outputs are never inputs."""
    result=[Path(c['_config_path'])]
    out=path(c,'results');data=path(c,'data')
    source_specs=ROOT/'config/sources.yaml'
    acquisition_deps={
        'data_manifest':[],
        'download':['zenodo_record.json','download_manifest.tsv'],
        'verify_checksums':['download_manifest.tsv','download_receipts.tsv'],
        'extract_archives':['verified_files.tsv']}
    if stage in acquisition_deps:
        result.append(source_specs)
        result.extend(data/'metadata'/name for name in acquisition_deps[stage])
        if stage in {'verify_checksums','extract_archives'}:
            mf=data/'metadata/download_manifest.tsv'
            if mf.exists():result.extend(Path(r['local_path']) for r in rows(mf))
        return sorted({p for p in result if p.is_file()})
    deps={
        'inspect_schemas':[],
        'build_isolate_manifest':['qc/data_inventory.tsv'],
        'validate_annotations':['qc/isolate_manifest.tsv'],
        'select_single_copy_core_genes':['qc/isolate_manifest.tsv','qc/annotation_validation.tsv'],
        'select_pilot_subset':['qc/isolate_manifest.tsv','orthogroups/core.json','orthogroups/membership_filtering_summary.tsv'],
        'extract_CDS':['qc/isolate_manifest.tsv','orthogroups/selection.json','orthogroups/selected_genes.tsv'],
        'extract_upstream':['qc/isolate_manifest.tsv','orthogroups/selection.json','orthogroups/selected_genes.tsv'],
        'align_CDS':['orthogroups/selection.json'],
        'align_upstream':['orthogroups/selection.json','upstream/upstream_qc.tsv'],
        'calculate_divergence':['orthogroups/selection.json','cds/alignment_status.tsv','upstream/alignment_status.tsv'],
        'analyze_repeats':['qc/isolate_manifest.tsv','upstream/upstream_coordinates.tsv'],
        'validate_synteny':['qc/isolate_manifest.tsv','orthogroups/selection.json','orthogroups/selected_genes.tsv',
            'upstream/upstream_coordinates.tsv','upstream/alignment_status.tsv','repeats/upstream_repeats.tsv'],
        'classify_loci':['orthogroups/selection.json','orthogroups/gene_coordinates.tsv',
            'orthogroups/pre_upstream_filtering_summary.tsv','upstream/upstream_coordinates.tsv','repeats/upstream_repeats.tsv'],
        'statistical_analysis':['candidates/master_locus_table.tsv','statistics/sensitivity_locus.tsv'],
        'generate_figures':['candidates/master_locus_table.tsv','statistics/classification_counts.tsv','statistics/sensitivity_summary.tsv'],
        'candidate_visualization':['candidates/master_locus_table.tsv','qc/isolate_manifest.tsv',
            'orthogroups/selection.json','upstream/upstream_coordinates.tsv'],
        'generate_report':['orthogroups/selection.json','candidates/master_locus_table.tsv',
            'statistics/prevalence.tsv','synteny/locus_evidence.tsv'],
        'graph_validation':['qc/isolate_manifest.tsv','upstream/upstream_coordinates.tsv','candidates/master_locus_table.tsv'],
        'test_fixtures':[]}
    result.extend(out/name for name in deps[stage])
    if stage=='inspect_schemas':
        for base in ['raw','external']:
            result.extend(p for p in (data/base).rglob('*') if p.is_file())
    if stage=='build_isolate_manifest':result.append(ROOT/c['samples'])
    if stage=='select_single_copy_core_genes':
        result.append(data/'raw'/c['orthology']['file'])
        if c['orthology'].get('id_map'):result.append(ROOT/c['orthology']['id_map'])
    if stage=='select_pilot_subset' and c['mode']!='full':
        result.append(ROOT/c['distances'] if c.get('distances') else data/'raw/Pe_variant_matrix.phy')
    if stage in {'validate_annotations','select_single_copy_core_genes','select_pilot_subset',
                 'extract_CDS','extract_upstream','analyze_repeats','candidate_visualization','graph_validation'}:
        manifest=out/'qc/isolate_manifest.tsv'
        if manifest.exists():
            for row in rows(manifest):result.extend(Path(row[k]) for k in ['fasta','gff','repeats'])
    patterns={
        'align_CDS':['cds/sequences/*.fa'], 'align_upstream':['upstream/sequences/*.fa'],
        'calculate_divergence':['cds/alignments/*.fa','upstream/alignments/*.fa'],
        'validate_synteny':['upstream/alignments/*.fa'],
        'classify_loci':['variation/intermediate/*.json'],
        'candidate_visualization':['variation/intermediate/*.json','upstream/alignments/*.fa']}
    for pattern in patterns.get(stage,[]):result.extend(out.glob(pattern))
    if stage=='graph_validation':
        result.extend([ROOT/c['graph']['file'],ROOT/c['graph']['path_mapping']])
    if stage=='generate_report':result.extend(path(c,'figures').glob('*.png'))
    return sorted({p for p in result if p.is_file()})


def execute(c,stage):
    inputs=stage_inputs(c,stage)
    record_stage(c,stage,'running','Stage started',inputs=inputs)
    try:
        outputs=ACTIONS[stage](c)
        if not outputs or any(not Path(p).is_file() for p in outputs):
            raise ValueError('Stage returned missing or no verified outputs')
        # Output hashes define the current stage boundary; later audit detects tampering.
        input_filtered=[p for p in inputs if p not in {Path(o) for o in outputs}]
        record_stage(c,stage,'complete','Outputs written and checksummed',outputs,inputs=input_filtered)
        print(f'{stage}: complete ({len(outputs)} files)',flush=True)
        return outputs
    except Exception as error:
        record_stage(c,stage,'failed',f'{type(error).__name__}: {error}',inputs=inputs)
        raise


def preflight(c,executor='local'):
    info={'timestamp_utc':now(),'os':platform.platform(),'architecture':platform.machine(),
        'cpu_count':os.cpu_count(),'disk_free_bytes':shutil.disk_usage(ROOT).free,
        'executor':executor,'network_assumption':'API/download stages need internet; compute stages use staged files',
        'tools':{x:shutil.which(x) for x in ['bash','python3','snakemake','mafft','Rscript','micromamba','conda','sbatch','qsub','gh']}}
    if executor=='slurm' and not shutil.which('sbatch'):raise ValueError('SLURM requested but sbatch unavailable')
    if executor=='pbs':
        if not shutil.which('qsub'):raise ValueError('PBS requested but qsub unavailable')
        p=subprocess.run(['qsub','--version'],capture_output=True,text=True)
        if 'pbs' not in (p.stdout+p.stderr).lower():raise ValueError('qsub is not verified as PBS Pro; inspect scheduler version')
    write_json(path(c,'reports')/'provenance/preflight.json',info)
    print(json.dumps(info,indent=2))
    return info


def run_test(c):
    if c['mode']!='test':raise ValueError('Test runner requires config/test.yaml')
    for action in ['test_fixtures','inspect_schemas','build_isolate_manifest','validate_annotations',
        'select_single_copy_core_genes','select_pilot_subset','extract_CDS','extract_upstream',
        'align_CDS','align_upstream','calculate_divergence','analyze_repeats','validate_synteny',
        'classify_loci','statistical_analysis','generate_figures','candidate_visualization','generate_report']:
        execute(c,action)
    audit(c)


def provenance(c):
    reports=path(c,'reports');records=[]
    for p in sorted((reports/'provenance').glob('*.json')):
        obj=read_json(p)
        if 'stage' in obj:records.append(obj)
    inputs={}
    for obj in records:inputs.update(obj['input_sha256'])
    table(reports/'provenance/input_checksums.tsv',[{'path':k,'sha256':v} for k,v in sorted(inputs.items())],['path','sha256'])
    software(c)
    obj={'mode':c['mode'],'timestamp_utc':now(),'git_commit':git_hash(),'config':dict(c,_config_path=rel(c['_config_path'])),
        'seed':c['random_seed'],'stages':[{k:r[k] for k in ['stage','status','timestamp_utc']} for r in records],
        'complete':all(r['status']=='complete' for r in records) and bool(records)}
    write_json(reports/'provenance/run_manifest.json',obj)
    return obj


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='config/pilot.yaml')
    sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('stage');s.add_argument('name',choices=sorted(ACTIONS))
    sub.add_parser('audit');sub.add_parser('test');sub.add_parser('provenance')
    q=sub.add_parser('preflight');q.add_argument('--executor',choices=['local','pbs','slurm'],default='local')
    a=p.parse_args();c=load_config(ROOT/a.config)
    try:
        if a.command=='stage':execute(c,a.name)
        elif a.command=='test':run_test(c)
        elif a.command=='audit':print(json.dumps(audit(c),indent=2))
        elif a.command=='provenance':print(json.dumps(provenance(c),indent=2))
        else:preflight(c,a.executor)
    except Exception as error:
        print(f'ERROR: {type(error).__name__}: {error}',file=sys.stderr)
        sys.exit(1)

if __name__=='__main__':main()
