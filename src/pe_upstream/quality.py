"""Explicit sample mapping, file/schema inventory and annotation validation."""
import csv
import re
from collections import Counter,defaultdict
from pathlib import Path
from .common import ROOT,path,table,rows,write_json,read_json
from .core import fasta,features,load_genes,coding_sequence,load_orthology,open_text


def source_path(c,name):
    return path(c,'data')/'raw'/name


def inspect(c):
    inventory=[]
    for base in ['raw','external']:
        for p in sorted((path(c,'data')/base).rglob('*')):
            if not p.is_file():
                continue
            kind='binary';header='';feature_types=''
            if p.suffix.lower() in {'.fa','.fna','.fasta','.fas','.gff','.gff3','.gtf','.tsv','.bed','.gz'}:
                try:
                    with open_text(p) as f:
                        preview=[]
                        for _ in range(50):
                            line=f.readline()
                            if not line:break
                            preview.append(line.rstrip())
                    data=[l for l in preview if l and not l.startswith('#')]
                    header=(data[0] if data else '').replace('\t',' | ')[:800]
                    if preview and preview[0].startswith('>'):kind='fasta'
                    elif data and len(data[0].split('\t'))==9:
                        kind='gff3' if '=' in data[0].split('\t')[-1] else 'gtf_or_gff'
                        feature_types=','.join(sorted({l.split('\t')[2] for l in data if len(l.split('\t'))==9}))
                    elif p.suffix=='.tsv':kind='tsv'
                    elif p.suffix=='.bed':kind='bed'
                    else:kind='text_unclassified'
                except (UnicodeError,OSError):kind='binary_or_compressed_unknown'
            inventory.append({'path':str(p),'bytes':p.stat().st_size,'detected_format':kind,
                              'first_data_row':header,'feature_types_preview':feature_types})
    out=path(c,'results')/'qc'
    p=table(out/'data_inventory.tsv',inventory,['path','bytes','detected_format','first_data_row','feature_types_preview'])
    return [p]


def build_manifest(c):
    src=ROOT/c['samples'];data=rows(src)
    fields={'isolate','haplotype','assembly_role','fasta','gff','repeats'}
    if not data or not fields.issubset(data[0]):
        raise ValueError(f'Populate {src} with columns {sorted(fields)} after inspecting the archives')
    if len(data)!=c['core_population_size'] or len({r['isolate'] for r in data})!=len(data):
        raise ValueError(f"Expected {c['core_population_size']} unique isolates in the core manifest")
    for row in data:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',row['isolate']):
            raise ValueError('Unsafe/blank isolate identifier')
        if not row['haplotype'] or row['assembly_role'] not in {'haploid_primary','haploid_consensus','phased_haplotype','synthetic'}:
            raise ValueError(f"Specify source-supported haplotype and assembly_role for {row['isolate']}")
        if c['mode']!='test' and row['assembly_role']=='synthetic':
            raise ValueError('Synthetic sample forbidden in biological mode')
        for key in ['fasta','gff','repeats']:
            if key == 'repeats' and not row[key] and row.get('repeat_status') == 'missing_from_release':
                continue
            p=Path(row[key]);p=p if p.is_absolute() else ROOT/p
            if not p.is_file():raise ValueError(f'Missing {key}: {p}')
            row[key]=str(p.resolve())
    if len({r['assembly_role'] for r in data})>1:
        raise ValueError('Mixed assembly representations require separate analyses, not pooled haplotypes')
    dest=path(c,'results')/'qc/isolate_manifest.tsv'
    table(dest,data)
    return [dest]


def n50(lengths):
    running=0;total=sum(lengths)
    for length in sorted(lengths,reverse=True):
        running+=length
        if running>=total/2:return length
    return 0


def validate(c):
    manifest=rows(path(c,'results')/'qc/isolate_manifest.tsv')
    validation=[];summaries=[];errors=[]
    for sample in manifest:
        isolate=sample['isolate'];genome=fasta(sample['fasta'])
        models=load_genes(sample['gff'])  # Duplicate gene IDs fail here, never silently overwritten.
        repeats=list(features(sample['repeats'])) if sample['repeats'] else []
        for f in repeats:
            if f.chrom not in genome or f.end>len(genome[f.chrom]):
                errors.append(f'{isolate}:repeat_contig_or_coordinate_mismatch')
        lengths=[len(seq) for seq in genome.values()]
        summaries.append({'isolate':isolate,'haplotype':sample['haplotype'],
            'assembly_role':sample['assembly_role'],'contig_count':len(genome),
            'total_bp':sum(lengths),'N50':n50(lengths),'max_contig_bp':max(lengths),
            'ambiguous_bp':sum(sum(ch not in 'ACGT' for ch in seq) for seq in genome.values()),
            'gap_runs':sum(len(re.findall('N+',seq)) for seq in genome.values()),
            'annotated_genes':len(models),'repeat_features':len(repeats) if sample['repeats'] else None})
        for g in models.values():
            flags=[]
            if g.chrom not in genome or g.end>len(genome.get(g.chrom,'')):
                flags.append('invalid_contig_or_coordinates');errors.append(f'{isolate}:{g.id}:coordinates')
            else:
                try:coding_sequence(g,genome)
                except ValueError as e:flags.append(str(e))
            validation.append({'isolate':isolate,'gene_id':g.id,'transcript_id':g.transcript,
                'chrom':g.chrom,'start_0':g.start,'end_0_exclusive':g.end,'strand':g.strand,
                'status':'pass' if not flags else 'flagged','flags':';'.join(flags)})
    out=path(c,'results')/'qc'
    outputs=[table(out/'annotation_validation.tsv',validation),table(out/'assembly_summary.tsv',summaries)]
    if errors:raise ValueError('Annotation/assembly mismatch: '+';'.join(errors[:20]))
    report=['# Input QC','',f'{len(manifest)} mapped isolates; {len(validation)} annotated genes inspected.',
            '',f'{sum(r["status"]!="pass" for r in validation)} gene models flagged; flags are retained for filtering.',
            '', 'Assembly representation is explicit in isolate_manifest.tsv. Gap runs and ambiguous bases '
            'are assembly properties, not biological deletions. Repeat classifications require source inspection.']
    from .common import atomic_text
    outputs.append(atomic_text(out/'qc_report.md','\n'.join(report)+'\n'))
    return outputs


def normalize_orthology(c,manifest):
    opts=c['orthology'];p=source_path(c,opts['file'])
    data=rows(p)
    if not data:raise ValueError('Empty orthogroup membership table')
    headers=list(data[0]);isolates={r['isolate'] for r in manifest}
    fmt=opts['format'];group=opts.get('group_column');iso=opts.get('isolate_column');gene=opts.get('gene_column')
    if fmt=='auto':
        if {'orthogroup','isolate','gene_id'}.issubset(headers):
            fmt,group,iso,gene='long','orthogroup','isolate','gene_id'
        elif len(headers)==len(isolates)+1 and set(headers[1:])==isolates:
            fmt,group='wide',headers[0]
        else:
            raise ValueError(f'Unrecognized orthology schema {headers}; inspect and explicitly configure parser')
    required={group,iso,gene} if fmt=='long' else {group}
    if not required.issubset(headers):raise ValueError(f'Missing configured orthology columns: {required-set(headers)}')
    if fmt=='wide' and set(headers)-{group}!=isolates:
        raise ValueError('Wide table isolate headers must exactly match the manifest; remove summary columns explicitly')
    id_map={}
    if opts.get('id_map'):
        for row in rows(ROOT/opts['id_map']):
            key=(row['isolate'],row['source_id'])
            if key in id_map:raise ValueError('Nonunique explicit gene-ID map')
            id_map[key]=row['gene_id']
    normalized=[];seen=set()
    for row in data:
        og=row[group].strip()
        if not og:raise ValueError('Empty orthogroup ID')
        cells=[(row[iso],row[gene])] if fmt=='long' else [(s,row[s]) for s in sorted(isolates)]
        for isolate,cell in cells:
            if isolate not in isolates:raise ValueError(f'Unknown orthology isolate: {isolate}')
            for gid in re.split(r'[,;\s]+',cell.strip()):
                if gid in {'','-','.','NA'}:continue
                if gid.isdecimal():raise ValueError('Numeric gene cell: copy-count tables cannot provide gene memberships')
                if id_map:
                    if (isolate,gid) not in id_map:raise ValueError(f'Missing explicit ID mapping: {isolate}:{gid}')
                    gid=id_map[(isolate,gid)]
                key=(isolate,gid)
                if key in seen:raise ValueError(f'Duplicate/multiply assigned orthology gene: {key}')
                seen.add(key);normalized.append({'orthogroup':og,'isolate':isolate,'gene_id':gid})
    out=path(c,'data')/'processed/orthology.tsv'
    table(out,normalized,['orthogroup','isolate','gene_id'])
    # Validate every normalized membership, not just eventual candidate genes.
    by_iso=defaultdict(set)
    for row in normalized:by_iso[row['isolate']].add(row['gene_id'])
    for sample in manifest:
        genes=load_genes(sample['gff']);missing=by_iso[sample['isolate']]-set(genes)
        if missing:raise ValueError(f"Unmapped gene IDs for {sample['isolate']}: {sorted(missing)[:10]}")
    return out,normalized


def core_groups(c):
    manifest=rows(path(c,'results')/'qc/isolate_manifest.tsv')
    out,normalized=normalize_orthology(c,manifest)
    isolates=[r['isolate'] for r in manifest]
    membership_core=load_orthology(out,isolates)
    annotation_status={(r['isolate'],r['gene_id']):r['status']=='pass'
        for r in rows(path(c,'results')/'qc/annotation_validation.tsv')}
    core={og:members for og,members in membership_core.items()
        if all(annotation_status[(isolate,gid)] for isolate,gid in members.items())}
    membership=defaultdict(Counter)
    for row in normalized:membership[row['orthogroup']][row['isolate']]+=1
    summary=[]
    for og,counts in sorted(membership.items()):
        summary.append({'orthogroup':og,'isolates_present':len(counts),'copies_total':sum(counts.values()),
            'single_copy_all_manifest':int(og in membership_core),
            'high_confidence_cds_all_manifest':int(og in core),
            'all_19_core':int(len(isolates)==19 and og in membership_core),
            'reason':('pass' if og in core else 'gene_model_QC' if og in membership_core else 'missing_isolate' if len(counts)<len(isolates) else 'multiple_copies')})
    dest=path(c,'results')
    outputs=[out,table(dest/'qc/orthogroup_summary.tsv',summary)]
    core_rows=[{'orthogroup':og,'isolate':s,'gene_id':gid,
                'all_19_core':int(len(isolates)==19),'single_copy_core_status':'all-19-core' if len(isolates)==19 else 'synthetic-core'}
               for og,members in sorted(core.items()) for s,gid in members.items()]
    outputs.append(table(dest/'orthogroups/single_copy_core_genes.tsv',core_rows,
                   ['orthogroup','isolate','gene_id','all_19_core','single_copy_core_status']))
    outputs.append(write_json(dest/'orthogroups/core.json',core))
    filtering=[{'stage':'membership_groups','remaining':len(membership),'removed':0},
               {'stage':'single_copy_all_manifest','remaining':len(membership_core),'removed':len(membership)-len(membership_core)},
               {'stage':'CDS_QC_all_manifest','remaining':len(core),'removed':len(membership_core)-len(core)}]
    outputs.append(table(dest/'orthogroups/membership_filtering_summary.tsv',filtering))
    if not core:raise ValueError('No single-copy core orthogroups after identifier validation')
    return outputs
