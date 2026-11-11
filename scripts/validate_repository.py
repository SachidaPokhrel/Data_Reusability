#!/usr/bin/env python3
"""Offline syntax, configuration, Markdown-link and required-file checks."""
import ast
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,load_config,write_json,now


def main():
    errors=[];counts={'python_files':0,'shell_files':0,'configs':0,'markdown_links':0}
    for base in ['src','scripts','hpc','tests']:
        for p in (ROOT/base).rglob('*.py'):
            if 'work' in p.relative_to(ROOT).parts:continue
            ast.parse(p.read_text(),filename=str(p));counts['python_files']+=1
            if any(line.rstrip()!=line for line in p.read_text().splitlines()):errors.append('Trailing whitespace: '+str(p.relative_to(ROOT)))
    for p in [ROOT/'bootstrap.sh',ROOT/'run.sh',ROOT/'hpc/submit.sh',ROOT/'hpc/scratch.sh']:
        result=subprocess.run(['bash','-n',str(p)],capture_output=True,text=True)
        counts['shell_files']+=1
        if result.returncode:errors.append(result.stderr)
    for mode in ['test','pilot','full']:
        load_config(ROOT/f'config/{mode}.yaml');counts['configs']+=1
    for base in ['docs','config','envs','hpc','workflow','scripts']:
        for p in (ROOT/base).rglob('*.md'):
            for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
                if target.startswith(('http:','https:','#','mailto:')):continue
                target=unquote(target.split('#')[0]);counts['markdown_links']+=1
                if not (p.parent/target).exists():errors.append(f'Broken local link {p.relative_to(ROOT)}: {target}')
    for name in ['README.md','LICENSE','CITATION.cff','.gitignore','.gitattributes','CHANGELOG.md',
                 'bootstrap.sh','run.sh','Snakefile','envs/workflow.yaml','envs/analysis.yaml','envs/r-analysis.yaml']:
        if not (ROOT/name).is_file():errors.append('Missing required file: '+name)
    lint='not_run_ruff_unavailable'
    if shutil.which('ruff'):
        subprocess.run(['ruff','check','src','scripts','hpc','tests'],cwd=ROOT,check=True);lint='passed'
    result={'timestamp_utc':now(),'status':'pass' if not errors else 'fail','counts':counts,
        'errors':errors,'ruff':lint,'Snakemake_parse':'not_checked_by_this_script',
        'R_syntax':'requires_installed_R_and_scientific_integration'}
    write_json(ROOT/'reports/provenance/repository_validation.json',result)
    print(json.dumps(result,indent=2))
    if errors:sys.exit(1)

if __name__=='__main__':main()
