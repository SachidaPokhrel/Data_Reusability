"""Diversity-based isolate selection and outcome-independent locus sampling."""
import itertools
import random
from collections import defaultdict
from pathlib import Path
from .common import ROOT,path,rows,table,read_json,write_json,digest
from .core import fasta,load_genes,coding_sequence,gene_index,synteny_signatures,upstream,DNA


def distance_matrix(c,isolates):
    p=Path(c['distances']) if c.get('distances') else path(c,'data')/'raw/Pe_variant_matrix.phy'
    p=p if p.is_absolute() else ROOT/p
    distances={}
    if c['distances_format']=='tsv':
        for r in rows(p):
            a,b=r['isolate_a'],r['isolate_b'];d=float(r['distance'])
            if a not in isolates or b not in isolates or not 0<=d<=1:
                raise ValueError('Invalid isolate or distance in genotype distance table')
            k=tuple(sorted((a,b)))
            if k in distances:raise ValueError('Duplicate unordered isolate pair in distance table')
            distances[k]=d
    else:
        from Bio import AlignIO
        errors=[];alignment=None
        for fmt in ['phylip-relaxed','phylip']:
            try:
                alignment=AlignIO.read(str(p),fmt);break
            except ValueError as e:errors.append(str(e))
        if alignment is None:raise ValueError('Invalid PHYLIP matrix: '+';'.join(errors))
        seqs={r.id:str(r.seq).upper() for r in alignment}
        if set(seqs)!=set(isolates):raise ValueError('PHYLIP isolate names must exactly match the sample manifest')
        for a,b in itertools.combinations(isolates,2):
            sites=[(x,y) for x,y in zip(seqs[a],seqs[b]) if x in DNA and y in DNA]
            if not sites:raise ValueError('No callable genotype sites for isolate pair')
            distances[tuple(sorted((a,b)))]=sum(x!=y for x,y in sites)/len(sites)
    expected={tuple(sorted(k)) for k in itertools.combinations(isolates,2)}
    if not expected.issubset(distances):raise ValueError('Incomplete genotype distance matrix')
    return distances,p


def farthest_isolates(isolates,distances,n):
    if n>len(isolates):raise ValueError('Requested too many isolates')
    if n==len(isolates):return sorted(isolates)
    pairs=sorted(itertools.combinations(sorted(isolates),2))
    first=max(pairs,key=lambda pair:distances[tuple(sorted(pair))])
    selected=list(first[:n])
    while len(selected)<n:
        remaining=sorted(set(isolates)-set(selected))
        best=max(remaining,key=lambda s:min(distances[tuple(sorted((s,t)))] for t in selected))
        selected.append(best)
    return selected


def select(c):
    out=path(c,'results');manifest=rows(out/'qc/isolate_manifest.tsv')
    isolates=[r['isolate'] for r in manifest]
    if c['mode']=='full':
        gate=read_json(path(dict(c,mode='pilot'),'reports')/'analysis_gate.json')
        if gate['status']!='pass':raise ValueError('Full mode blocked: real pilot QC gate has not passed')
        for filename,sha in gate['input_files'].items():
            if digest(ROOT/filename)!=sha:raise ValueError('Pilot inputs changed; rerun pilot before full mode')
        selected=sorted(isolates);basis='all manifest isolates'
    else:
        distances,p=distance_matrix(c,isolates)
        selected=farthest_isolates(isolates,distances,c['number_of_isolates'])
        basis=f'farthest-point sampling from validated genotype distances: {p}'
    core=read_json(out/'orthogroups/core.json')
    models={};genomes={};indices={};sigs={}
    for sample in manifest:
        s=sample['isolate']
        if s not in selected:continue
        models[s]=load_genes(sample['gff']);genomes[s]=fasta(sample['fasta'])
        indices[s]=gene_index(models[s]);mapping={m[s]:og for og,m in core.items()}
        sigs[s]=synteny_signatures(models[s],mapping)
    qc=[];eligible=[]
    for og,members in sorted(core.items()):
        flags=[];signatures=[sigs[s].get(members[s]) for s in selected]
        if None in signatures or len(set(signatures))!=1:flags.append('unresolved_synteny')
        for s in selected:
            g=models[s][members[s]]
            try:
                seq=coding_sequence(g,genomes[s])
                if len(seq)<c['min_cds_bases']:flags.append(f'{s}:short_cds')
                if sum(ch not in DNA for ch in seq)/len(seq)>c['max_ambiguous_fraction']:
                    flags.append(f'{s}:ambiguous_cds')
                # Extraction QC is recorded separately; upstream failure stays unresolved,
                # rather than disappearing from the locus denominator.
            except ValueError as e:flags.append(f'{s}:{e}')
        qc.append({'orthogroup':og,'eligible_cds_synteny':int(not flags),'flags':';'.join(flags)})
        if not flags:eligible.append(og)
    chosen=sorted(eligible)
    if c['number_of_genes'] and len(chosen)>c['number_of_genes']:
        chosen=sorted(random.Random(c['random_seed']).sample(chosen,c['number_of_genes']))
    if not chosen:raise ValueError('No valid coding/synteny loci; inspect QC before changing filters')
    isolate_rows=[]
    for rank,isolate in enumerate(selected):
        marginal=None
        if c['mode']!='full' and rank and len(selected)<len(isolates):
            marginal=min(distances[tuple(sorted((isolate,t)))] for t in selected[:rank])
        isolate_rows.append({'isolate':isolate,'selection_rank':rank+1,'method':basis,
            'min_distance_to_previously_selected':marginal,
            'reason':'full cohort' if len(selected)==len(isolates) else 'member of maximum-distance initial pair' if rank<2 else 'maximized minimum genotype distance to prior choices'})
    selected_rows=[];coords=[]
    for og in chosen:
        for s in selected:
            g=models[s][core[og][s]]
            selected_rows.append({'orthogroup':og,'isolate':s,'gene_id':g.id,
                'transcript_id':g.transcript,'all_19_core':int(len(isolates)==19),
                'pilot_six_core':int(len(selected)==6),'selection_basis':basis})
            coords.append({'orthogroup':og,'isolate':s,'gene_id':g.id,'transcript_id':g.transcript,
                'chrom':g.chrom,'start':g.start,'end':g.end,'strand':g.strand,
                'contig_length':len(genomes[s][g.chrom]),
                'distance_to_contig_end':min(g.start,len(genomes[s][g.chrom])-g.end)})
    filtering=rows(out/'orthogroups/membership_filtering_summary.tsv')
    filtering += [{'stage':'CDS_synteny_selected_isolates','remaining':len(eligible),'removed':len(core)-len(eligible)},
                  {'stage':'random_sample','remaining':len(chosen),'removed':len(eligible)-len(chosen)}]
    files=[table(out/'orthogroups/selected_isolates.tsv',isolate_rows),
        table(out/'orthogroups/pre_upstream_filtering_summary.tsv',filtering),
        table(out/'orthogroups/selected_genes.tsv',selected_rows),
        table(out/'orthogroups/gene_coordinates.tsv',coords),
        table(out/'orthogroups/locus_selection_qc.tsv',qc),
        write_json(out/'orthogroups/selection.json',{'isolates':selected,'orthogroups':chosen,
            'isolate_selection_basis':basis,'seed':c['random_seed'],
            'all_manifest_core_genes':len(core),'eligible_cds_synteny_genes':len(eligible),
            'sampled_before_upstream_qc':True})]
    return files
