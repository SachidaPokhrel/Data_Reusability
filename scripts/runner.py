#!/usr/bin/env python3
"""Portable dispatcher; real runs use Snakemake, offline tests have a documented fallback."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,load_config,read_json,path,digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['test','pilot','full','validate','report','audit','preflight','download'])
    p.add_argument('--cores',type=int,default=2)
    p.add_argument('--executor',choices=['local','pbs','slurm'],default='local')
    p.add_argument('--config',help='Alternate config; defaults to mode-specific config')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--offline-test',action='store_true',help='Run the same stage functions without Snakemake/R/MAFFT using TEST fixtures only')
    p.add_argument('--account');p.add_argument('--queue');p.add_argument('--dependency')
    p.add_argument('--scratch',action='store_true',help='Use configured SCRATCH_ROOT in a scheduler allocation')
    a=p.parse_args()
    if a.cores<1:p.error('--cores must be positive')
    mode=a.mode if a.mode in {'test','pilot','full'} else 'pilot'
    config=Path(a.config or f'config/{mode}.yaml').resolve();c=load_config(config)
    if a.mode in {'test','pilot','full'} and c['mode']!=a.mode:p.error('Mode and config mode disagree')
    logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    if a.mode=='preflight':
        subprocess.run([sys.executable,'scripts/pe.py','--config',str(config),'preflight','--executor',a.executor],check=True);return
    if a.executor!='local':
        if a.mode not in {'pilot','full'}:p.error('HPC submission accepts pilot or full modes only')
        command=[sys.executable,str(ROOT/'hpc/submit.py'),a.mode,'--executor',a.executor,
                 '--cores',str(a.cores),'--config',str(config)]
        for key in ['account','queue','dependency']:
            if getattr(a,key):command+=['--'+key,getattr(a,key)]
        if a.scratch:command.append('--scratch')
        if a.dry_run:command.append('--dry-run')
        subprocess.run(command,check=True);return
    py=sys.executable
    cli=[py,str(ROOT/'scripts/pe.py'),'--config',str(config)]
    conda_root=Path(os.environ.get('CONDA_ROOT',ROOT/'.conda'))
    snake=conda_root/'envs/workflow/bin/snakemake'
    if not snake.exists():snake=Path(shutil.which('snakemake') or '/nonexistent/snakemake')
    env=os.environ.copy();env['PATH']=str(snake.parent)+os.pathsep+env['PATH']
    if a.mode=='preflight':subprocess.run(cli+['preflight'],check=True);return
    if a.mode=='audit':subprocess.run(cli+['audit'],check=True);return
    if a.mode=='validate':
        subprocess.run([py,'-m','unittest','discover','-s','tests','-v'],check=True)
        subprocess.run([py,'scripts/validate_repository.py'],check=True);return
    if a.mode=='test' and (a.offline_test or not snake.exists()):
        if not a.offline_test:print('Snakemake unavailable: running explicitly SYNTHETIC offline stage harness.',flush=True)
        subprocess.run([py,'-m','unittest','discover','-s','tests','-v'],check=True)
        subprocess.run(cli+['test'],check=True);return
    if not snake.exists():
        raise ValueError('Snakemake not installed. Run ./bootstrap.sh on a network-enabled machine; see offline HPC instructions.')
    if a.mode=='full':
        pc=load_config(ROOT/'config/pilot.yaml');gate=read_json(path(pc,'reports')/'analysis_gate.json')
        if gate['status']!='pass':raise ValueError('Full mode requires a successful real pilot QC gate')
        for filename,sha in gate['input_files'].items():
            if digest(filename)!=sha:raise ValueError('Pilot source data changed; rerun the pilot before full mode')
    target=[]
    if a.mode=='report':target=['generate_report','provenance','audit']
    if a.mode=='download':target=['inspect_schemas']
    command=[str(snake),'--snakefile',str(ROOT/'Snakefile'),'--configfile',str(config),
        '--config','project_config='+str(config),'--cores',str(a.cores),'--use-conda',
        '--conda-prefix',str(conda_root/'snakemake'),'--rerun-incomplete','--printshellcmds',
        '--resources','mem_mb='+str(c['memory_mb']),'disk_mb='+str(c['temporary_disk_mb'])]
    if a.dry_run:command.append('--dry-run')
    command += ['--'] + target
    subprocess.run(command,cwd=ROOT,env=env,check=True)

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.CalledProcessError) as e:sys.exit(f'ERROR: {e}')
