"""Configuration, atomic I/O, checksums, provenance, and status contracts."""
from __future__ import annotations
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[2])).resolve()


def now():
    return datetime.now(timezone.utc).isoformat()


def load_config(path):
    """Load a mode configuration and reject incompatible thresholds or TEST settings."""
    path = Path(path).resolve()
    text = path.read_text()
    try:
        c = json.loads(text)
    except json.JSONDecodeError:
        import yaml
        c = yaml.safe_load(text)
    c['_config_path'] = str(path)
    if os.environ.get('PE_THREADS'):
        c['threads'] = int(os.environ['PE_THREADS'])
    if c['mode'] not in {'test', 'pilot', 'full'}:
        raise ValueError('mode must be test, pilot or full')
    for key in ['number_of_isolates', 'core_population_size', 'threads']:
        if not isinstance(c[key], int) or c[key] <= 0:
            raise ValueError(f'{key} must be a positive integer')
    for key in ['min_pair_callable_fraction', 'min_valid_pair_fraction',
                'max_cds_divergence', 'min_upstream_divergence', 'max_ambiguous_fraction']:
        if not 0 <= c[key] <= 1:
            raise ValueError(f'{key} outside [0,1]')
    if not c['upstream_windows_bp'] or any(not isinstance(w, int) or w <= 0
                                          for w in c['upstream_windows_bp']):
        raise ValueError('upstream_windows_bp requires positive integers')
    if len(set(c['upstream_windows_bp'])) != len(c['upstream_windows_bp']):
        raise ValueError('duplicate upstream windows')
    if c['test_stub_alignment'] and c['mode'] != 'test':
        raise ValueError('Synthetic alignment stub is forbidden for biological runs')
    if c['mode'] != 'test' and c['core_population_size'] != 19:
        raise ValueError('Main biological question requires a 19-isolate core definition')
    return c


def path(c, key):
    """Resolve a storage base and add the mode subdirectory for generated products."""
    env = {'data': 'DATA_ROOT', 'results': 'RESULTS_ROOT', 'reports': 'REPORTS_ROOT',
           'figures': 'FIGURES_ROOT'}[key]
    p = Path(os.environ.get(env, c['paths'][key]))
    p = p if p.is_absolute() else ROOT / p
    if key in {'results', 'reports', 'figures'}:
        p = p / c['mode']
    p.mkdir(parents=True, exist_ok=True)
    return p


def digest(p, algorithm='sha256'):
    """Hash a file in chunks so large sequence inputs need not fit in memory."""
    h = hashlib.new(algorithm)
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_text(p, text):
    """Write a temporary sibling and rename it so incomplete content is not published."""
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + '.part')
    tmp.write_text(text)
    tmp.replace(p)
    return p


def write_json(p, obj):
    return atomic_text(p, json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n')


def read_json(p):
    return json.loads(Path(p).read_text())


def table(p, rows, columns=None):
    """Write a TSV atomically, preserving blank values and explicit column order."""
    rows = list(rows)
    columns = columns or (list(rows[0]) if rows else None)
    if not columns:
        raise ValueError(f'No table columns for empty output {p}')
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + '.part')
    with open(tmp, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    tmp.replace(p)
    return p


def rows(p):
    """Read a tab-delimited table into dictionaries keyed by its header."""
    with open(p, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def rel(p):
    """Use a project-relative path when possible; preserve external paths explicitly."""
    p = Path(p).resolve()
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)


def git_hash():
    p = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'],
                       capture_output=True, text=True)
    return p.stdout.strip() if p.returncode == 0 else None


def code_checksums():
    """Hash the executable workflow files used for stage provenance."""
    files = sorted((ROOT / 'src').rglob('*.py')) + sorted((ROOT / 'scripts').rglob('*.py'))
    files += sorted((ROOT / 'scripts').rglob('*.R')) + [ROOT / 'Snakefile']
    files += sorted((ROOT / 'workflow').rglob('*.smk'))
    return {rel(p): digest(p) for p in files if p.is_file()}


def record_stage(c, stage, status, details, outputs=(), inputs=()):
    """Record stage status, input and output hashes, configuration, and software context."""
    dest = path(c, 'reports') / 'provenance' / f'{stage}.json'
    obj = {'stage': stage, 'status': status, 'timestamp_utc': now(), 'mode': c['mode'],
           'details': details, 'git_commit': git_hash(), 'config': c,
           'config_sha256': digest(c['_config_path']), 'code_sha256': code_checksums(),
           'input_sha256': {rel(p): digest(p) for p in inputs if Path(p).is_file()},
           'output_sha256': {rel(p): digest(p) for p in outputs if Path(p).is_file()},
           'command': sys.argv}
    tools={'python':platform.python_version()}
    if stage.startswith('align_'):tools['mafft']=version(['mafft','--version'])
    if stage in {'statistical_analysis','generate_figures'}:tools['R']=version(['Rscript','--version'])
    obj['stage_software_versions']=tools
    obj['config'] = dict(c, _config_path=rel(c['_config_path']))
    write_json(dest, obj)
    return obj


def version(command):
    if not shutil.which(command[0]):
        return 'UNAVAILABLE'
    try:
        p = subprocess.run(command, text=True, capture_output=True, timeout=20)
        return (p.stdout + p.stderr).strip()[:400]
    except (OSError, subprocess.TimeoutExpired) as e:
        return str(e)


def software(c):
    """Record tool versions from the current environment and successful stage receipts."""
    out = path(c, 'reports') / 'provenance' / 'software_versions.tsv'
    versions={
        'python': sys.version, 'platform': platform.platform(), 'mafft': version(['mafft','--version']),
        'R': version(['Rscript','--version']), 'snakemake': version(['snakemake','--version']),
        'minimap2': version(['minimap2','--version']), 'git': version(['git','--version'])}
    for receipt in (path(c,'reports')/'provenance').glob('*.json'):
        obj=read_json(receipt)
        if obj.get('status')!='complete':continue
        for key,value in obj.get('stage_software_versions',{}).items():
            if value!='UNAVAILABLE':versions[key]=value
    table(out,[{'software':k,'version':v} for k,v in versions.items()])
    return out


def audit(c):
    """Verify recorded inputs, outputs and implementation hashes against the current files."""
    problems = []
    checked = 0
    for p in sorted((path(c, 'reports') / 'provenance').glob('*.json')):
        obj = read_json(p)
        if 'stage' not in obj:
            continue
        if obj['status'] != 'complete':
            problems.append(f'{p.name}: {obj["status"]}')
            continue
        for group in ['input_sha256','output_sha256','code_sha256']:
            for name, expected in obj[group].items():
                target = ROOT / name
                checked += 1
                if not target.is_file() or digest(target) != expected:
                    problems.append(f'{p.name}: stale/missing {name}')
        cfg = ROOT / obj['config']['_config_path']
        if not cfg.exists() or digest(cfg) != obj['config_sha256']:
            problems.append(f'{p.name}: changed config')
    if checked == 0:
        problems.append('No completed stages with checksummed inputs/outputs')
    result = {'status': 'pass' if not problems else 'fail', 'checks': checked,
              'problems': problems, 'timestamp_utc': now()}
    write_json(path(c,'reports')/'provenance'/'audit.json',result)
    if problems:
        raise ValueError('\n'.join(problems))
    return result


def markdown_page(text,title):
    """Small escaped Markdown renderer for generated report headings/tables/images."""
    import html
    import re
    def inline(s):
        escaped=html.escape(s)
        return re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',escaped)
    result=[];in_table=False;table_rows=0;in_list=False
    for line in text.splitlines()+['']:
        if line.startswith('|'):
            if not in_table:
                result.append('<table>');in_table=True;table_rows=0
            cells=[x.strip() for x in line.strip('|').split('|')]
            if all(re.fullmatch(r'[:\- ]+',x) for x in cells):continue
            tag='th' if table_rows==0 else 'td'
            result.append('<tr>'+''.join(f'<{tag}>{inline(x)}</{tag}>' for x in cells)+'</tr>')
            table_rows+=1;continue
        if in_table:result.append('</table>');in_table=False
        if line.startswith('- '):
            if not in_list:result.append('<ul>');in_list=True
            result.append('<li>'+inline(line[2:])+'</li>');continue
        if in_list:result.append('</ul>');in_list=False
        if not line.strip():continue
        heading=re.match(r'^(#{1,6}) (.*)',line)
        picture=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if heading:
            level=len(heading[1]);result.append(f'<h{level}>'+inline(heading[2])+f'</h{level}>')
        elif picture:
            result.append('<figure><img src="'+html.escape(picture[2],quote=True)+'" alt="'+html.escape(picture[1],quote=True)+'"><figcaption>'+inline(picture[1])+'</figcaption></figure>')
        else:result.append('<p>'+inline(line)+'</p>')
    style='body{font:16px system-ui;max-width:1100px;margin:3em auto;padding:0 1em;line-height:1.6}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:.5em}img{max-width:100%}h1,h2{line-height:1.2}'
    return '<!doctype html><html lang="en"><meta charset="utf-8"><title>'+html.escape(title)+'</title><style>'+style+'</style><main>'+'\n'.join(result)+'</main></html>'
