"""Callable sequence distances, gap-event summaries and paired locus classification."""
import csv
import itertools
import statistics
from collections import defaultdict
from .common import path,rows,table,read_json,write_json
from .core import fasta,DNA
from .sequences import file_id

PROTEIN=frozenset('ACDEFGHIKLMNPQRSTVWY')


def metrics(a,b,large_indel=50,alphabet=DNA):
    if len(a)!=len(b):raise ValueError('Unequal aligned lengths')
    callable_sites=mismatch=gap_bp=events=large=0;run=0;direction=None
    for x,y in list(zip(a,b))+[('?', '?')]:
        d='a' if x=='-' and y in alphabet else 'b' if y=='-' and x in alphabet else None
        if d is not None:
            gap_bp+=1
            if d==direction:run+=1
            else:
                if run>=large_indel:large+=1
                events+=1;run=1
        else:
            if run>=large_indel:large+=1
            run=0
        direction=d
        if x in alphabet and y in alphabet:
            callable_sites+=1;mismatch+=x!=y
    longest=max(sum(ch in alphabet for ch in a),sum(ch in alphabet for ch in b))
    denom=callable_sites+gap_bp
    return {'alignment_length':len(a),'callable_sites':callable_sites,
        'callable_fraction':callable_sites/longest if longest else 0,
        'identity':1-mismatch/callable_sites if callable_sites else None,
        'snp_density_per_kb':1000*mismatch/callable_sites if callable_sites else None,
        'substitutions':mismatch,'indel_events':events,'indel_bases':gap_bp,
        'indel_events_per_kb':1000*events/denom if denom else None,
        'gap_fraction':gap_bp/denom if denom else None,'large_indel_events':large,
        'combined_divergence':(mismatch+gap_bp)/denom if denom else None,
        'length_difference':abs(sum(ch!='-' for ch in a)-sum(ch!='-' for ch in b)),
        'unaligned_fraction':1-callable_sites/longest if longest else 1}


def alignment_summary(aligned,alphabet=DNA):
    seqs=list(aligned.values());length=len(seqs[0])
    return {'alignment_length':length,'variable_sites':sum(len({s[i] for s in seqs if s[i] in alphabet})>1
            for i in range(length)),
        'length_min':min(sum(ch!='-' for ch in s) for s in seqs),
        'length_max':max(sum(ch!='-' for ch in s) for s in seqs),
        'missing_fraction':sum(ch not in alphabet for s in seqs for ch in s)/(len(seqs)*length) if length else 1,
        'ambiguous_fraction':sum(ch not in alphabet and ch!='-' for s in seqs for ch in s)/(len(seqs)*length) if length else 1}


def calculate(c):
    out=path(c,'results');selection=read_json(out/'orthogroups/selection.json')
    statuses={}
    for kind in ['cds','upstream']:
        statuses.update({(r['orthogroup'],r['region']):r['alignment_status'] for r in rows(out/kind/'alignment_status.tsv')})
    work=out/'variation/intermediate';work.mkdir(parents=True,exist_ok=True)
    outputs=[];summaries=[]
    for og in selection['orthogroups']:
        result={}
        for region in ['cds','protein']+[f'upstream_{w}' for w in c['upstream_windows_bp']]:
            kind='upstream' if region.startswith('upstream') else 'cds'
            status=statuses[(og,region)]
            if status.startswith('unresolved'):
                result[region]={'status':status,'pairs':[]};continue
            aligned=fasta(out/kind/'alignments'/f'{file_id(og)}.{region}.fa')
            alphabet=PROTEIN if region=='protein' else DNA
            pp=[]
            for a,b in itertools.combinations(selection['isolates'],2):
                m=metrics(aligned[a],aligned[b],c['large_indel_bp'],alphabet)
                pp.append(dict(isolate_a=a,isolate_b=b,valid=int(m['callable_fraction']>=c['min_pair_callable_fraction']),**m))
            summary=alignment_summary(aligned,alphabet)
            result[region]={'status':status,'pairs':pp,'summary':summary}
            summaries.append(dict(orthogroup=og,region=region,**summary))
        p=work/f'{file_id(og)}.json';write_json(p,result);outputs.append(p)
    outputs.append(table(out/'variation/alignment_summary.tsv',summaries))
    # Stream pairwise output from per-locus records; full mode does not retain millions of rows in RAM.
    for kind in ['cds','upstream']:
        dest=out/kind/f'{kind}_pairwise.tsv';tmp=dest.with_suffix('.part');writer=None
        with open(tmp,'w',newline='') as handle:
            for og in selection['orthogroups']:
                data=read_json(work/f'{file_id(og)}.json')
                for region,obj in data.items():
                    if (region.startswith('upstream'))!=(kind=='upstream'):continue
                    for pair in obj['pairs']:
                        row=dict(orthogroup=og,region=region,**pair)
                        if writer is None:
                            writer=csv.DictWriter(handle,fieldnames=list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                        writer.writerow(row)
            if writer is None:handle.write('orthogroup\tregion\tisolate_a\tisolate_b\tvalid\n')
        tmp.replace(dest);outputs.append(dest)
    return outputs


def classify(c):
    out=path(c,'results');selection=read_json(out/'orthogroups/selection.json')
    coords=rows(out/'upstream/upstream_coordinates.tsv');repeat=rows(out/'repeats/upstream_repeats.tsv')
    by_coord=defaultdict(list);by_rep=defaultdict(list)
    for r in coords:by_coord[(r['orthogroup'],int(r['window']))].append(r)
    for r in repeat:by_rep[(r['orthogroup'],int(r['window']))].append(r)
    def mean(values):
        values=[float(v) for v in values if v is not None and v!='']
        return statistics.mean(values) if values else None
    master=[];sensitivity=[]
    coords_gene=rows(out/'orthogroups/gene_coordinates.tsv')
    genes=defaultdict(list)
    for r in coords_gene:genes[r['orthogroup']].append(r)
    for og in selection['orthogroups']:
        data=read_json(out/'variation/intermediate'/f'{file_id(og)}.json')
        cds={(r['isolate_a'],r['isolate_b']):r for r in data['cds']['pairs']}
        protein=data['protein']['pairs']
        aa=mean([r['identity'] for r in protein if r['valid']])
        for w in c['upstream_windows_bp']:
            up=data[f'upstream_{w}'];pairup={(r['isolate_a'],r['isolate_b']):r for r in up['pairs']}
            matched=[k for k in pairup if pairup[k]['valid'] and cds[k]['valid']]
            coverage=len(matched)/len(cds) if cds else 0
            independent_cd=mean([v['combined_divergence'] for v in cds.values() if v['valid']])
            independent_callable=sum(v['valid'] for v in cds.values())/len(cds)
            independently_conserved=independent_callable>=c['min_valid_pair_fraction'] and independent_cd is not None and independent_cd<=c['max_cds_divergence']
            reason=[];cd=uv=identity=None
            if up['status'].startswith('unresolved'):reason.append(up['status'])
            if coverage<c['min_valid_pair_fraction']:reason.append('insufficient_matched_callable_pairs')
            if matched:
                cd=mean([cds[k]['combined_divergence'] for k in matched])
                uv=mean([pairup[k]['combined_divergence'] for k in matched])
                identity=mean([pairup[k]['identity'] for k in matched])
                if identity<c['min_upstream_identity_for_homology']:reason.append('unresolved_upstream_homology')
            cr=by_coord[(og,w)];rr=by_rep[(og,w)]
            missing=mean([r['ambiguous_fraction'] for r in cr]);ref=genes[og][0]
            status='E'
            # The middle band is explicitly unresolved instead of forcing borderline loci to be conserved.
            if not reason:
                cc=independently_conserved and cd<=c['max_cds_divergence'];up_variable=uv>=c['min_upstream_divergence']
                up_conserved=uv<=c.get('max_conserved_upstream_divergence',0.01)
                if not up_variable and not up_conserved:reason.append('intermediate_upstream_divergence')
                else:status=('B' if cc else 'D') if up_variable else ('A' if cc else 'C')
            row={'orthogroup_id':og,'window':w,'gene_ids':';'.join(r['isolate']+':'+r['gene_id'] for r in genes[og]),
                'isolate_count':len(selection['isolates']),
                'single_copy_core_status':'all-19-core' if c['mode']!='test' else 'synthetic-core',
                'mean_CDS_identity':mean([cds[k]['identity'] for k in matched]),
                'median_CDS_identity':statistics.median(cds[k]['identity'] for k in matched) if matched else None,
                'protein_identity':aa,'CDS_SNP_density':mean([cds[k]['snp_density_per_kb'] for k in matched]),
                'CDS_divergence':cd,
                'CDS_divergence_independent':independent_cd,
                'coding_conserved_independent':int(independently_conserved),
                'mean_upstream_identity':identity,
                'upstream_divergence':uv,'upstream_SNP_density':mean([pairup[k]['snp_density_per_kb'] for k in matched]),
                'upstream_indel_burden':mean([pairup[k]['gap_fraction'] for k in matched]),
                'upstream_indel_events':mean([pairup[k]['indel_events'] for k in matched]),
                'upstream_large_indel_events':mean([pairup[k]['large_indel_events'] for k in matched]),
                'upstream_length_variability':max(int(r['actual_length']) for r in cr)-min(int(r['actual_length']) for r in cr),
                'intergenic_length_variability':max(int(r['intergenic_length']) for r in cr)-min(int(r['intergenic_length']) for r in cr),
                'upstream_alignment_coverage':mean([pairup[k]['callable_fraction'] for k in matched]),
                'matched_pair_fraction':coverage,'upstream_repeat_fraction':mean([r['repeat_fraction'] for r in rr]),
                'nearest_TE_distance':mean([r['nearest_TE_distance'] for r in rr]),
                'nearest_repeat_distance':mean([r['nearest_repeat_distance'] for r in rr]),
                'local_repeat_density':mean([r['local_repeat_density'] for r in rr]),
                'intergenic_length':mean([r['intergenic_length'] for r in cr]),
                'missing_data_fraction':missing,'chrom':ref['chrom'],
                'distance_to_contig_end':int(ref['distance_to_contig_end']),
                'synteny_status':'conserved_immediate_core_neighbors',
                'assembly_quality_flags':';'.join(sorted({r['flags'] for r in cr if r['flags']})),
                'candidate_status':status,'unresolved_reason':';'.join(reason),
                'candidate_screen':int(status=='B'),'biological_validation':'pending' if status=='B' else 'not_selected'}
            master.append(row)
            # Sensitivity classification only varies prespecified divergence thresholds;
            # extraction, callability, homology QC are held fixed.
            base_resolved=not [r for r in reason if r!='intermediate_upstream_divergence'] and cd is not None and uv is not None
            for ct in c['sensitivity_cds']:
                for ut in c['sensitivity_upstream']:
                    sensitivity.append({'orthogroup':og,'window':w,'cds_threshold':ct,'upstream_threshold':ut,
                        'resolved':int(base_resolved),'coding_conserved':int(base_resolved and cd<=ct and independent_cd<=ct),
                        'candidate':int(base_resolved and cd<=ct and independent_cd<=ct and uv>=ut)})
    filtering=[dict(r,window='all') for r in rows(out/'orthogroups/pre_upstream_filtering_summary.tsv')]
    for w in c['upstream_windows_bp']:
        data=[r for r in master if r['window']==w]
        resolved=sum(r['candidate_status']!='E' for r in data)
        filtering.append({'stage':'resolved_locus_classification','window':w,
            'remaining':resolved,'removed':len(data)-resolved})
    return [table(out/'candidates/master_locus_table.tsv',master),
            table(out/'statistics/sensitivity_locus.tsv',sensitivity),
            table(out/'orthogroups/filtering_summary.tsv',filtering,['stage','window','remaining','removed'])]
