#!/usr/bin/env python3
"""Relocate scratch output references only after matching destination checksums."""
import argparse
import shutil
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pe_upstream.common import ROOT,read_json,write_json,digest,rel,now

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reports',required=True);p.add_argument('--mode',required=True)
p.add_argument('--mapping',nargs=2,action='append',required=True,metavar=('OLD_BASE','NEW_BASE'))
a=p.parse_args();directory=Path(a.reports)/a.mode/'provenance';log=[]
for file in directory.glob('*.json'):
    obj=read_json(file)
    if 'stage' not in obj:continue
    original=digest(file);changed=False
    for group in ['input_sha256','output_sha256']:
        new={}
        for name,sha in obj[group].items():
            target=ROOT/name
            for oldbase,newbase in a.mapping:
                if target.is_relative_to(Path(oldbase)):
                    dest=Path(newbase)/target.relative_to(Path(oldbase))
                    if not dest.exists() or digest(dest)!=sha:raise ValueError(f'Relocation checksum failed: {dest}')
                    target=dest;changed=True;break
            new[rel(target)]=sha
        obj[group]=new
    if changed:
        saved=directory/'relocation_originals'/file.name;saved.parent.mkdir(exist_ok=True)
        shutil.copy2(file,saved);write_json(file,obj)
        log.append({'stage':obj['stage'],'original_sha256':original,'relocated_sha256':digest(file)})
write_json(directory/'relocation.json',{'timestamp_utc':now(),'mappings':a.mapping,'records':log})
