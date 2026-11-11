"""Strand-aware extraction, sequence QC, deterministic translation and alignment."""
import hashlib
import json
import shutil
import itertools
import subprocess
from collections import defaultdict
from pathlib import Path
from .common import path,rows,table,read_json,write_json,digest,version,atomic_text
from .core import fasta,write_fasta,load_genes,coding_sequence,revcomp,DNA

CODONS=[''.join(x) for x in itertools.product('TCAG',repeat=3)]
AA=dict(zip(CODONS,'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def translate(seq):
    protein=''.join(AA.get(seq[i:i+3],'X') for i in range(0,len(seq),3))
    return protein[:-1] if protein.endswith('*') else protein


def file_id(og):
    return hashlib.sha256(og.encode()).hexdigest()[:20]


def selected_data(c):
    out=path(c,'results');selection=read_json(out/'orthogroups/selection.json')
    selected=rows(out/'orthogroups/selected_genes.tsv')
    manifest={r['isolate']:r for r in rows(out/'qc/isolate_manifest.tsv')}
    return selection,selected,manifest


def cds(c):
    selection,selected,manifest=selected_data(c)
    by_iso=defaultdict(list)
    for r in selected:by_iso[r['isolate']].append(r)
    coding=defaultdict(dict);proteins=defaultdict(dict);qc=[]
    for s,entries in by_iso.items():
        genome=fasta(manifest[s]['fasta']);genes=load_genes(manifest[s]['gff'])
        for r in entries:
            gene=genes[r['gene_id']];seq=coding_sequence(gene,genome);prot=translate(seq)
            coding[r['orthogroup']][s]=seq;proteins[r['orthogroup']][s]=prot
            qc.append({'orthogroup':r['orthogroup'],'isolate':s,'gene_id':gene.id,
                'transcript_id':gene.transcript,'cds_length':len(seq),'protein_length':len(prot),
                'ambiguous_fraction':sum(ch not in DNA for ch in seq)/len(seq),
                'internal_stop_codons':prot.count('*'),'translation_status':'consistent_standard_code'})
    outputs=[];directory=path(c,'results')/'cds/sequences';directory.mkdir(parents=True,exist_ok=True)
    for og in selection['orthogroups']:
        for label,data in [('cds',coding),('protein',proteins)]:
            p=directory/f'{file_id(og)}.{label}.fa';write_fasta(p,data[og]);outputs.append(p)
    outputs+=[table(path(c,'results')/'cds/cds_qc.tsv',qc)]
    return outputs


def extract_one(g,genome,models,window):
    """Clip at assembly/neighbor boundaries, retain QC instead of hiding short regions."""
    length=len(genome[g.chrom]);desired=(g.start-window,g.start) if g.strand=='+' else (g.end,g.end+window)
    a=max(0,desired[0]);b=min(length,desired[1]);flags=[]
    if desired[0]<0 or desired[1]>length:flags.append('contig_boundary')
    overlap=any(n.id!=g.id and n.chrom==g.chrom and max(n.start,a)<min(n.end,b) for n in models.values())
    if overlap:flags.append('neighbor_overlap')
    neighbors=[n for n in models.values() if n.id!=g.id and n.chrom==g.chrom]
    if g.strand=='+':
        ends=[n.end for n in neighbors if n.end<=g.start]
        distance=g.start-max(ends) if ends else g.start
        if ends:a=max(a,max(ends))
    else:
        starts=[n.start for n in neighbors if n.start>=g.end]
        distance=min(starts)-g.end if starts else length-g.end
        if starts:b=min(b,min(starts))
    if a>b:a=b
    seq=genome[g.chrom][a:b]
    seq=revcomp(seq) if g.strand=='-' else seq
    ambiguous=sum(ch not in DNA for ch in seq)
    if ambiguous:flags.append('ambiguous_bases')
    if len(seq)<window:flags.append('short_upstream')
    return seq,{'chrom':g.chrom,'strand':g.strand,'gene_start':g.start,'gene_end':g.end,
        'requested_length':window,'actual_length':len(seq),'start':a,'end':b,
        'upstream_gene_distance':distance,'intergenic_length':distance,
        'overlaps_neighbor':int(overlap),'ambiguous_bases':ambiguous,
        'ambiguous_fraction':ambiguous/len(seq) if seq else 1,
        'gap_runs':len(__import__('re').findall('N+',seq)),
        'distance_to_contig_end':min(g.start,length-g.end),'flags':';'.join(flags)}


def upstream(c):
    selection,selected,manifest=selected_data(c)
    by_iso=defaultdict(list)
    for r in selected:by_iso[r['isolate']].append(r)
    seqs=defaultdict(dict);coords=[];all_sequences={}
    for s,entries in by_iso.items():
        genome=fasta(manifest[s]['fasta']);genes=load_genes(manifest[s]['gff'])
        # Limit neighbor scans to a contig to avoid cross-chromosome comparisons.
        per_chrom=defaultdict(dict)
        for g in genes.values():per_chrom[g.chrom][g.id]=g
        for r in entries:
            g=genes[r['gene_id']]
            for w in c['upstream_windows_bp']:
                seq,info=extract_one(g,genome,per_chrom[g.chrom],w)
                valid=(info['actual_length']==w and not info['overlaps_neighbor'] and
                       info['ambiguous_fraction']<=c['max_ambiguous_fraction'])
                info.update(orthogroup=r['orthogroup'],isolate=s,gene_id=g.id,window=w,
                            annotation_anchor='gene_feature_boundary',extractable=int(valid))
                coords.append(info);seqs[(r['orthogroup'],w)][s]=seq
                all_sequences[f'{file_id(r["orthogroup"])}|{s}|{w}']=seq
    directory=path(c,'results')/'upstream/sequences';directory.mkdir(parents=True,exist_ok=True)
    outputs=[]
    for (og,w),seq in seqs.items():
        p=directory/f'{file_id(og)}.upstream_{w}.fa';write_fasta(p,seq);outputs.append(p)
    out=path(c,'results')/'upstream'
    outputs.append(table(out/'upstream_coordinates.tsv',coords))
    outputs.append(table(out/'upstream_qc.tsv',[{k:r[k] for k in
        ['orthogroup','isolate','window','extractable','actual_length','flags']} for r in coords]))
    write_fasta(out/'upstream_sequences.fasta',all_sequences);outputs.append(out/'upstream_sequences.fasta')
    return outputs


def run_alignment(source,dest,c,protein=False):
    original=fasta(source)
    if c['test_stub_alignment']:
        # This deliberately simple TEST-only padder is not a biological aligner.
        maxlen=max(map(len,original.values()))
        write_fasta(dest,{s:seq.ljust(maxlen,'-') for s,seq in original.items()})
    else:
        temp=dest.with_suffix('.part')
        args=['mafft','--thread',str(c['threads']),'--auto']
        if protein:args.append('--amino')
        args.append(str(source))
        with open(temp,'w') as out,open(dest.with_suffix('.log'),'w') as err:
            subprocess.run(args,stdout=out,stderr=err,check=True)
        temp.replace(dest)
    aligned=fasta(dest)
    if set(aligned)!=set(original) or len(set(map(len,aligned.values())))!=1:
        raise ValueError('Invalid alignment rows/identifiers')
    for s,seq in aligned.items():
        if seq.replace('-','')!=original[s]:raise ValueError('Aligner modified ungapped sequence')


def align(c,kind):
    selection=read_json(path(c,'results')/'orthogroups/selection.json');out=path(c,'results')
    directory=out/kind/'alignments';directory.mkdir(parents=True,exist_ok=True)
    cache=out/kind/'intermediate/alignment_cache';cache.mkdir(parents=True,exist_ok=True)
    tool_version='TEST_STUB' if c['test_stub_alignment'] else version(['mafft','--version'])
    qc=rows(out/'upstream/upstream_qc.tsv') if kind=='upstream' else []
    valid={(r['orthogroup'],int(r['window']),r['isolate']):r['extractable']=='1' for r in qc}
    outputs=[];statuses=[]
    for og in selection['orthogroups']:
        regions=['cds','protein'] if kind=='cds' else [f'upstream_{w}' for w in c['upstream_windows_bp']]
        for region in regions:
            ok=kind=='cds' or all(valid[(og,int(region.split('_')[-1]),s)] for s in selection['isolates'])
            if not ok:
                statuses.append({'orthogroup':og,'region':region,'alignment_status':'unresolved_extraction'})
                continue
            source=out/kind/'sequences'/f'{file_id(og)}.{region}.fa'
            dest=directory/f'{file_id(og)}.{region}.fa'
            key=hashlib.sha256((digest(source)+tool_version+str(c['threads'])+digest(__file__)).encode()).hexdigest()
            cached=cache/(key+'.fa');receipt=cache/(key+'.json')
            valid_cache=cached.exists() and receipt.exists() and digest(cached)==read_json(receipt)['output_sha256']
            if not valid_cache:
                run_alignment(source,cached,c,region=='protein')
                write_json(receipt,{'input_sha256':digest(source),'output_sha256':digest(cached),
                    'software':tool_version,'threads':c['threads']})
            tmp=dest.with_suffix('.part');shutil.copyfile(cached,tmp);tmp.replace(dest);outputs.append(dest)
            statuses.append({'orthogroup':og,'region':region,'alignment_status':'synthetic_stub' if c['test_stub_alignment'] else 'mafft'})
    outputs.append(table(out/kind/'alignment_status.tsv',statuses))
    return outputs
