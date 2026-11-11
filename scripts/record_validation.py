#!/usr/bin/env python3
"""Execute offline tests and repository checks, preserving their actual outcomes."""
import platform
import re
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,write_json,now,git_hash

logs=ROOT/'logs';logs.mkdir(exist_ok=True)
command=[sys.executable,'-m','unittest','discover','-s','tests','-v']
result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
text=result.stdout+result.stderr;(logs/'tests_final.log').write_text(text);print(text)
run=re.search(r'Ran (\d+) tests?',text)
skip=len(re.findall(r'\.\.\. skipped ',text))
fail=len(re.findall(r'\.\.\. (?:FAIL|ERROR)',text))
count=int(run[1]) if run else 0
obj={'timestamp_utc':now(),'git_commit':git_hash(),'python':platform.python_version(),
    'command':['python3']+command[1:],'tests_run':count,'passed':count-skip-fail,
    'skipped':skip,'failed':fail,'exit_code':result.returncode,
    'status':'pass' if result.returncode==0 and count else 'fail',
    'scope':'unit and complete offline synthetic stage integration; no biological data',
    'real_tool_tests':'MAFFT conditional; R and Snakemake execution not available here'}
write_json(ROOT/'reports/provenance/test_results.json',obj)
if result.returncode:sys.exit(result.returncode)
subprocess.run([sys.executable,'scripts/validate_repository.py'],cwd=ROOT,check=True)
