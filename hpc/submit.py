#!/usr/bin/env python3
"""Submit one explicit scheduler allocation running the Snakemake DAG locally inside it."""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,load_config,write_json,now
from pe_upstream.cli import preflight


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['pilot','full'])
    p.add_argument('--executor',choices=['pbs','slurm'],required=True)
    p.add_argument('--cores',type=int,default=8);p.add_argument('--config')
    p.add_argument('--account',default=os.environ.get('HPC_ACCOUNT'))
    p.add_argument('--queue',default=os.environ.get('HPC_QUEUE'))
    p.add_argument('--dependency',help='Scheduler job ID; submit with afterok dependency')
    p.add_argument('--scratch',action='store_true');p.add_argument('--dry-run',action='store_true')
    p.add_argument('--preflight',action='store_true')
    a=p.parse_args();config=Path(a.config or ROOT/f'config/{a.mode}.yaml').resolve();c=load_config(config)
    preflight(c,a.executor)
    if not a.account or not a.queue:raise ValueError('Supply site-approved --account/--queue or HPC_ACCOUNT/HPC_QUEUE; no allocation is assumed')
    if a.cores<1:raise ValueError('Invalid CPU allocation')
    for key,value in [('account',a.account),('queue',a.queue)]:
        if any(ch.isspace() for ch in value) or value.startswith('-'):raise ValueError(f'Invalid {key}')
    # Site-native queue inspection. Fail instead of submitting to an unknown queue.
    command=['sinfo','--noheader','--partition',a.queue] if a.executor=='slurm' else ['qstat','-Qf',a.queue]
    if not shutil.which(command[0]):raise ValueError(f'Queue inspection command unavailable: {command[0]}')
    checked=subprocess.run(command,capture_output=True,text=True,check=True)
    if not checked.stdout.strip():raise ValueError('Requested queue/partition could not be validated')
    conda_root=Path(os.environ.get('CONDA_ROOT',ROOT/'.conda')).resolve()
    if not (conda_root/'envs/workflow/bin/snakemake').is_file():raise ValueError('Pre-install shared environments with ./bootstrap.sh on an approved internet-enabled node')
    if a.scratch and (not os.environ.get('SCRATCH_ROOT') or not shutil.which('rsync')):
        raise ValueError('--scratch needs SCRATCH_ROOT and rsync')
    if a.preflight:print('Scheduler, queue, allocation settings and shared environment passed preflight');return
    logs=ROOT/'logs/hpc';logs.mkdir(parents=True,exist_ok=True)
    export={k:v for k,v in os.environ.items() if k in ['DATA_ROOT','RESULTS_ROOT','REPORTS_ROOT','FIGURES_ROOT','SCRATCH_ROOT']}
    export.update(PROJECT_ROOT=str(ROOT),CONDA_ROOT=str(conda_root))
    script=['#!/usr/bin/env bash','set -euo pipefail']
    for k,v in export.items():script.append(f'export {k}={shlex.quote(v)}')
    script.append('cd "$PROJECT_ROOT"')
    if a.scratch:
        script += [f'export PE_RUN_MODE={shlex.quote(a.mode)}','source hpc/scratch.sh']
    script+=['date -u',f'/usr/bin/time -v ./run.sh {shlex.quote(a.mode)} --config {shlex.quote(str(config))} --executor local --cores {a.cores}']
    job=logs/f'{a.mode}.{a.executor}.sh';job.write_text('\n'.join(script)+'\n');job.chmod(0o755)
    memory=int(os.environ.get('HPC_MEMORY_MB',c['memory_mb']));minutes=int(os.environ.get('HPC_RUNTIME_MIN',c['runtime_min']))
    if a.executor=='slurm':
        cmd=['sbatch','--parsable','--account',a.account,'--partition',a.queue,'--job-name','pe-'+a.mode,
            '--cpus-per-task',str(a.cores),'--mem',str(memory)+'M','--time',str(minutes),
            '--output',str(logs/f'{a.mode}.%j.out'),'--error',str(logs/f'{a.mode}.%j.err')]
        if a.dependency:cmd+=['--dependency','afterok:'+a.dependency]
    else:
        wall=f'{minutes//60:02}:{minutes%60:02}:00'
        cmd=['qsub','-A',a.account,'-q',a.queue,'-N','pe-'+a.mode,'-l',f'select=1:ncpus={a.cores}:mem={memory}mb',
            '-l','walltime='+wall,'-o',str(logs/f'{a.mode}.pbs.out'),'-e',str(logs/f'{a.mode}.pbs.err')]
        if a.dependency:cmd+=['-W','depend=afterok:'+a.dependency]
    cmd.append(str(job))
    receipt={'executor':a.executor,'timestamp_utc':now(),'command':cmd,
        'status':'prepared_not_submitted','execution_model':'single allocation; local Snakemake DAG manages dependencies',
        'tested_on_cluster':False}
    if not a.dry_run:
        result=subprocess.run(cmd,text=True,capture_output=True,check=True)
        receipt.update(status='submitted',scheduler_response=result.stdout.strip())
    write_json(ROOT/'reports/provenance/hpc_submission.json',receipt)
    print(shlex.join(cmd));print(receipt['status'])

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.CalledProcessError) as e:sys.exit(f'HPC PREFLIGHT/SUBMISSION FAILED: {e}')
