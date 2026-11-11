"""API-first metadata retrieval and verified, resumable public-data downloads."""
import hashlib
import json
import shutil
import stat
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from .common import ROOT, path, rows, table, read_json, write_json, digest, now


def get_json(url, attempts=3, timeout=30):
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers={'User-Agent':'Pe-effusa-reuse/0.2'})
            with urllib.request.urlopen(req,timeout=timeout) as r:
                return json.load(r)
        except (OSError, ValueError):
            if attempt + 1 == attempts:
                raise
            time.sleep(min(2**attempt,8))


def manifest(c):
    sources = read_json(ROOT/'config/sources.yaml')
    record = get_json(sources['zenodo_api'])
    if str(record.get('id')) != sources['record_id']:
        raise ValueError('Zenodo API returned an unexpected record identifier')
    required = set(sources['required_files'])
    if c['mode'] in {'pilot','full'} and not c.get('distances'):
        required.add('Pe_variant_matrix.phy')
    if c.get('graph',{}).get('enabled'):
        required.add('Pe.full.gfa.gz')
    available = {f['key']:f for f in record.get('files',[])}
    if not required.issubset(available):
        raise ValueError(f'Required files absent: {sorted(required-set(available))}')
    out = path(c,'data')/'metadata'
    write_json(out/'zenodo_record.json',record)
    data = []
    for name in sorted(required):
        f = available[name]
        url = f['links'].get('content') or f['links'].get('self')
        if not url or urllib.parse.urlparse(url).scheme != 'https':
            raise ValueError(f'No verified HTTPS download link for {name}')
        checksum = f.get('checksum','')
        algorithm, expected = checksum.split(':',1)
        hashlib.new(algorithm)
        if not expected or not isinstance(f.get('size'),int):
            raise ValueError(f'Incomplete checksum/size metadata for {name}')
        if Path(name).name != name:
            raise ValueError('Unsafe remote filename')
        data.append({'source':'Zenodo','record_id':str(record['id']),
            'version':record.get('metadata',{}).get('version','not specified'),
            'file_name':name,'file_size':f['size'],'checksum':checksum,'download_url':url,
            'download_date':'','local_path':str(path(c,'data')/'raw'/name),
            'validation_status':'pending'})
    table(out/'download_manifest.tsv',data)
    return [out/'zenodo_record.json',out/'download_manifest.tsv']


def download_one(url, dest, checksum, expected_size, attempts=3):
    """Use Range when supported; verify checksum before atomic final rename."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    algorithm, expected = checksum.split(':',1)
    if dest.exists() and dest.stat().st_size == expected_size and digest(dest,algorithm)==expected:
        return 'verified_existing'
    part = dest.with_name(dest.name+'.part')
    for attempt in range(attempts):
        try:
            offset = part.stat().st_size if part.exists() else 0
            if offset > expected_size:
                part.unlink()
                offset = 0
            headers = {'User-Agent':'Pe-effusa-reuse/0.2'}
            if offset:
                headers['Range']=f'bytes={offset}-'
            req = urllib.request.Request(url,headers=headers)
            with urllib.request.urlopen(req,timeout=90) as response:
                resume = response.status == 206 and offset > 0
                if resume and not response.headers.get('Content-Range','').startswith(f'bytes {offset}-'):
                    raise ValueError('Invalid Content-Range in resumed transfer')
                with open(part,'ab' if resume else 'wb') as out:
                    shutil.copyfileobj(response,out)
            if part.stat().st_size != expected_size:
                raise ValueError(f'Incomplete download: {part.stat().st_size}/{expected_size}')
            if digest(part,algorithm) != expected:
                part.unlink()
                raise ValueError('Checksum mismatch; partial file removed')
            part.replace(dest)
            return 'verified_downloaded'
        except (OSError, ValueError):
            if attempt+1==attempts:
                raise
            time.sleep(min(2**attempt,8))
    raise RuntimeError('Unreachable download state')


def download(c):
    source=path(c,'data')/'metadata/download_manifest.tsv'
    data=rows(source)
    for row in data:
        row['validation_status']=download_one(row['download_url'],row['local_path'],
            row['checksum'],int(row['file_size']))
        row['download_date']=now()
        # Retain each successful transfer even if a later transfer fails.
        table(path(c,'data')/'metadata/download_receipts.tsv',data)
    return [path(c,'data')/'metadata/download_receipts.tsv']


def verify(c):
    data=rows(path(c,'data')/'metadata/download_manifest.tsv')
    for row in data:
        p=Path(row['local_path']);algorithm,expected=row['checksum'].split(':',1)
        if not p.exists() or p.stat().st_size!=int(row['file_size']) or digest(p,algorithm)!=expected:
            raise ValueError(f'Input failed checksum/size verification: {p}')
        row['validation_status']='verified'
    p=path(c,'data')/'metadata/verified_files.tsv'
    table(p,data)
    return [p]


def safe_extract(archive, destination):
    destination=Path(destination).resolve()
    destination.mkdir(parents=True,exist_ok=True)
    listing=[]
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            target=(destination/info.filename).resolve()
            if not target.is_relative_to(destination) or '\\' in info.filename:
                raise ValueError(f'Unsafe ZIP path: {info.filename}')
            if stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError(f'ZIP symlink forbidden: {info.filename}')
            listing.append({'archive':str(archive),'member':info.filename,
                            'uncompressed_size':info.file_size,'compressed_size':info.compress_size})
        # Inspect and validate every member before any extraction.
        z.extractall(destination)
    return listing


def extract(c):
    inventory=[]
    for row in rows(path(c,'data')/'metadata/verified_files.tsv'):
        p=Path(row['local_path'])
        if p.suffix=='.zip':
            inventory+=safe_extract(p,path(c,'data')/'external'/p.stem)
    dest=path(c,'data')/'metadata/archive_inventory.tsv'
    table(dest,inventory,['archive','member','uncompressed_size','compressed_size'])
    return [dest]
