"""Repeat overlap and TE-specific features without assuming repeat equals TE."""
from collections import defaultdict,Counter
from .common import path,rows,table
from .core import features,IntervalIndex,fasta


def repeat_records(filename,te_classes):
    """Read repeat intervals; only recognized source classifications count as TE classes."""
    records=[]
    for f in features(filename):
        classification=next((f.attrs[k] for k in ['Classification','classification','Class','class','repeat_class'] if k in f.attrs),'unknown')
        if classification=='unknown' and '#' in f.attrs.get('Name',''):
            classification=f.attrs['Name'].split('#',1)[1]
        te=classification.split('/')[0] in te_classes
        records.append({'chrom':f.chrom,'start':f.start,'end':f.end,
            'class':classification,'family':f.attrs.get('Name',f.attrs.get('Target','unknown')),
            'is_TE':te})
    return records


def indices(records,te_only=False):
    """Merge overlapping repeat intervals to avoid double counting bases."""
    intervals=defaultdict(list)
    for r in records:
        if not te_only or r['is_TE']:intervals[r['chrom']].append((r['start'],r['end']))
    return IntervalIndex(intervals)


def analyze(c):
    """Summarize repeat overlap, proximity and density for each measured upstream interval."""
    out=path(c,'results');manifest={r['isolate']:r for r in rows(out/'qc/isolate_manifest.tsv')}
    coord=rows(out/'upstream/upstream_coordinates.tsv');by_iso=defaultdict(list)
    for r in coord:by_iso[r['isolate']].append(r)
    result=[];classes=[]
    for s,entries in by_iso.items():
        rr=repeat_records(manifest[s]['repeats'],c['repeat_te_classes']);all_ix=indices(rr);te_ix=indices(rr,True)
        genome=fasta(manifest[s]['fasta']);by_chrom=defaultdict(list)
        for r in rr:by_chrom[r['chrom']].append(r)
        for classification,n in Counter(r['class'] for r in rr).items():
            classes.append({'isolate':s,'classification':classification,'features':n,
                'is_known_TE_class':int(classification.split('/')[0] in c['repeat_te_classes'])})
        for r in entries:
            a,b=int(r['start']),int(r['end']);length=b-a;chrom=r['chrom']
            overlapping=[x for x in by_chrom[chrom] if max(a,x['start'])<min(b,x['end'])]
            radius=c['local_repeat_radius_bp'];start=max(0,int(r['gene_start'])-radius)
            end=min(len(genome[chrom]),int(r['gene_end'])+radius)
            result.append({'orthogroup':r['orthogroup'],'isolate':s,'window':r['window'],
                'repeat_overlaps':len(overlapping),'repeat_fraction':all_ix.overlap(chrom,a,b)/length if length else None,
                'TE_fraction':te_ix.overlap(chrom,a,b)/length if length else None,
                'nearest_TE_distance':te_ix.distance(chrom,int(r['gene_start']),int(r['gene_end'])),
                'nearest_repeat_distance':all_ix.distance(chrom,int(r['gene_start']),int(r['gene_end'])),
                'local_repeat_density':all_ix.overlap(chrom,start,end)/(end-start),
                'repeat_classes':';'.join(sorted({x['class'] for x in overlapping})),
                'repeat_families':';'.join(sorted({x['family'] for x in overlapping})),
                'TE_status':'classified_TE_annotations' if any(x['is_TE'] for x in rr) else 'TE_classification_unavailable'})
    return [table(out/'repeats/upstream_repeats.tsv',result),table(out/'repeats/repeat_class_inventory.tsv',classes)]
