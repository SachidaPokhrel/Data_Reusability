#!/usr/bin/env python3
"""Install only into project-local storage; verify official micromamba package SHA256."""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,now,write_json,digest,version
from pe_upstream.acquisition import get_json,download_one
PIN='2.0.5'


def target_platform():
    osname=platform.system();arch=platform.machine().lower()
    mapping={('Linux','x86_64'):'linux-64',('Linux','aarch64'):'linux-aarch64',
             ('Darwin','x86_64'):'osx-64',('Darwin','arm64'):'osx-arm64'}
    if (osname,arch) not in mapping:
        raise ValueError(f'Unsupported platform {osname}/{arch}; Windows requires WSL2 Linux')
    return mapping[(osname,arch)]


def install_manager(subdir,offline):
    existing=shutil.which('micromamba') or shutil.which('conda')
    if existing:return Path(existing)
    binary=ROOT/'tools/bin/micromamba'
    if binary.exists():return binary
    if offline:raise ValueError('No environment manager installed; supply a verified offline micromamba binary in tools/bin')
    metadata=get_json('https://api.anaconda.org/package/conda-forge/micromamba')
    candidates=[f for f in metadata['files'] if f.get('version')==PIN and
        f.get('attrs',{}).get('subdir')==subdir and f['basename'].endswith('.tar.bz2')]
    if not candidates:raise ValueError(f'Pinned micromamba {PIN} unavailable for {subdir}; no unverified fallback')
    f=sorted(candidates,key=lambda r:r['basename'])[-1]
    checksum=f.get('sha256') or f.get('attrs',{}).get('sha256')
    if not checksum or len(checksum)!=64:raise ValueError('Official package metadata lacks SHA256; refusing bootstrap')
    url=f['download_url'];url='https:'+url if url.startswith('//') else url
    if not url.startswith('https://'):raise ValueError('Bootstrap download must use HTTPS')
    archive=ROOT/'tools/cache'/Path(f['basename']).name
    download_one(url,archive,'sha256:'+checksum,int(f['size']))
    binary.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(archive,'r:bz2') as t:
        member=t.getmember('bin/micromamba')
        if not member.isfile():raise ValueError('Bootstrap archive executable is not a regular file')
        with t.extractfile(member) as source,open(binary.with_suffix('.part'),'wb') as dest:
            shutil.copyfileobj(source,dest)
    binary.with_suffix('.part').replace(binary);binary.chmod(0o755)
    write_json(ROOT/'reports/provenance/bootstrap_binary.json',{'version':PIN,'platform':subdir,
        'download_url':url,'package_sha256':checksum,'binary_sha256':digest(binary),
        'metadata_source':'https://api.anaconda.org/package/conda-forge/micromamba','timestamp_utc':now()})
    return binary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--offline',action='store_true',help='Use installed manager and locally cached packages only')
    p.add_argument('--check',action='store_true',help='Inspect prerequisites without installing')
    a=p.parse_args();subdir=target_platform()
    conda_root=Path(os.environ.get('CONDA_ROOT',ROOT/'.conda')).resolve()
    env=os.environ.copy();env['MAMBA_ROOT_PREFIX']=str(conda_root)
    receipt={'timestamp_utc':now(),'platform':subdir,'status':'running','steps':[]}
    output=ROOT/'reports/provenance/bootstrap_status.json'
    try:
        if a.check:
            receipt.update(status='inspection_only',python=sys.version,
                manager=shutil.which('micromamba') or shutil.which('conda'),bash=shutil.which('bash'))
            write_json(output,receipt);print(json.dumps(receipt,indent=2));return
        manager=install_manager(subdir,a.offline);receipt['manager']=str(manager)
        for name in ['workflow','analysis','r-analysis']:
            spec=ROOT/f'envs/{name}.yaml';prefix=conda_root/'envs'/name
            command=[str(manager),'create' if manager.name=='micromamba' else 'env','--help']
            lock=ROOT/f'envs/locks/{subdir}-{name}.explicit.txt'
            # Exact-build exports, once available, supersede a fresh solve.
            if lock.exists():
                command=[str(manager),'create','-y','-p',str(prefix),'-f',str(lock)]
            elif manager.name=='micromamba':
                command=[str(manager),'create','-y','-p',str(prefix),'-f',str(spec)]
            else:
                command=[str(manager),'env','create','-y','-p',str(prefix),'-f',str(spec)]
            if prefix.exists():
                stamp=prefix/'pe_spec.sha256'
                if not stamp.exists() or stamp.read_text().strip()!=digest(spec):
                    raise ValueError(f'Existing environment differs from spec: {prefix}; retain it and select a new CONDA_ROOT')
                receipt['steps'].append({'environment':name,'status':'reused_spec_match'})
            else:
                if a.offline:command.append('--offline')
                subprocess.run(command,env=env,check=True)
                (prefix/'pe_spec.sha256').write_text(digest(spec)+'\n')
                receipt['steps'].append({'environment':name,'status':'created','command':command})
            result=subprocess.run([str(manager),'list','-p',str(prefix),'--explicit'],
                env=env,capture_output=True,text=True,check=True)
            lock.parent.mkdir(parents=True,exist_ok=True);lock.write_text(result.stdout)
            # Snakemake recognizes platform-specific explicit environment pin files.
            (ROOT/f'envs/{name}.{subdir}.pin.txt').write_text(result.stdout)
        work=conda_root/'envs/workflow/bin'
        subprocess.run([str(work/'python'),'--version'],check=True)
        subprocess.run([str(work/'snakemake'),'--version'],check=True)
        subprocess.run([str(conda_root/'envs/analysis/bin/mafft'),'--version'],check=True)
        subprocess.run([str(conda_root/'envs/r-analysis/bin/Rscript'),'-e',
            'library(ggplot2);library(dplyr);library(sandwich);sessionInfo()'],check=True)
        # Resolve all rule environments while connected. No biological job executes.
        command=[str(work/'snakemake'),'--snakefile',str(ROOT/'Snakefile'),'--configfile',str(ROOT/'config/test.yaml'),
                 '--config','project_config='+str(ROOT/'config/test.yaml'),'--cores','1','--use-conda',
                 '--conda-prefix',str(conda_root/'snakemake'),'--conda-create-envs-only']
        env['PATH']=str(work)+os.pathsep+env['PATH']
        subprocess.run(command,cwd=ROOT,env=env,check=True)
        receipt['status']='complete'
    except Exception as error:
        receipt.update(status='failed',error=f'{type(error).__name__}: {error}')
        write_json(output,receipt)
        raise
    write_json(output,receipt)

if __name__=='__main__':
    try:main()
    except Exception as e:sys.exit(f'BOOTSTRAP FAILED: {e}\nSee docs/troubleshooting.md and reports/provenance/bootstrap_status.json')
