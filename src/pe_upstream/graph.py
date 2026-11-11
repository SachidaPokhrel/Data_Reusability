"""Optional GFA assembly-path validation with explicit path/contig mapping."""
import re
from .common import ROOT,path,rows,table
from .core import open_text,fasta,revcomp


def load_gfa(filename):
    segments={};paths={}
    with open_text(filename) as f:
        for line in f:
            fields=line.rstrip().split('\t')
            if fields[0]=='S':
                if fields[2]=='*':raise ValueError('GFA segment sequences absent; cannot validate sequence coordinates')
                if fields[1] in segments:raise ValueError('Duplicate GFA segment')
                segments[fields[1]]=fields[2].upper()
            elif fields[0]=='P':
                if fields[1] in paths:raise ValueError('Duplicate GFA path name')
                paths[fields[1]]=[(x[:-1],x[-1]) for x in fields[2].split(',')]
            elif fields[0]=='W':
                name='#'.join(fields[1:4])
                if fields[4]!='0':raise ValueError('Partial GFA walks require explicit offset handling; not supported here')
                if name in paths:raise ValueError('Multiple GFA walk chunks for one path need explicit coordinate conversion')
                paths[name]=[(segment,'+' if orientation=='>' else '-') for orientation,segment in re.findall(r'([><])([^><]+)',fields[6])]
    if not segments or not paths:raise ValueError('GFA lacks usable segments and P/W assembly paths')
    return segments,paths


def validate(c):
    options=c.get('graph',{})
    if not options.get('enabled'):raise ValueError('Enable graph validation only after supplying an explicit graph-path mapping')
    mapping=rows(ROOT/options['path_mapping']);segments,paths=load_gfa(ROOT/options['file'])
    manifest={r['isolate']:r for r in rows(path(c,'results')/'qc/isolate_manifest.tsv')}
    coordinates=rows(path(c,'results')/'upstream/upstream_coordinates.tsv')
    master=rows(path(c,'results')/'candidates/master_locus_table.tsv')
    selected=[]
    for r in sorted(master,key=lambda r:float(r['upstream_divergence'] or 0),reverse=True):
        if r['candidate_status']=='B' and r['orthogroup_id'] not in selected:selected.append(r['orthogroup_id'])
        if len(selected)>=c['max_candidate_loci']:break
    result=[];seen=set()
    for row in mapping:
        key=(row['isolate'],row['contig'])
        if key in seen:raise ValueError('Nonunique graph path mapping')
        seen.add(key)
        if row['isolate'] not in manifest or row['graph_path'] not in paths:raise ValueError('Unmapped sample or graph path')
        assembly=fasta(manifest[row['isolate']]['fasta'])
        if row['contig'] not in assembly:raise ValueError('Mapped assembly contig absent')
        parts=[]
        for seg,orientation in paths[row['graph_path']]:
            if seg not in segments or orientation not in {'+','-'}:raise ValueError('Invalid graph walk')
            parts.append(segments[seg] if orientation=='+' else revcomp(segments[seg]))
        sequence=''.join(parts)
        if sequence!=assembly[row['contig']]:
            raise ValueError('Graph path differs from mapped full assembly contig; coordinate equivalence is not established')
        for r in coordinates:
            if r['orthogroup'] in selected and r['isolate']==row['isolate'] and r['chrom']==row['contig']:
                result.append({'orthogroup':r['orthogroup'],'window':r['window'],'isolate':row['isolate'],
                    'graph_path':row['graph_path'],'assembly_contig':row['contig'],
                    'graph_supported_sequence':1,'coordinate_validation':'full_path_exact_assembly_match',
                    'TE_insertion_history':'not_inferred','evidence':'strong_assembly_path_sequence_support'})
    if not result:raise ValueError('No candidate upstream regions covered by supplied graph path mapping')
    return [table(path(c,'results')/'candidates/graph_validation.tsv',result)]
