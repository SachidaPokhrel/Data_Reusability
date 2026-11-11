"""Candidate evidence from syntenic assemblies; conservative repeat-associated indels."""
import itertools
from collections import defaultdict
from .common import path,rows,table,read_json
from .core import fasta,DNA
from .sequences import file_id


def gap_blocks(a,b,min_size):
    """Find long contiguous gap runs opposite callable DNA in an aligned pair."""
    blocks=[];start=None;direction=None
    for i,(x,y) in enumerate(list(zip(a,b))+[('?', '?')]):
        d='a_gap' if x=='-' and y in DNA else 'b_gap' if y=='-' and x in DNA else None
        if d!=direction:
            if start is not None and i-start>=min_size:blocks.append((start,i,direction))
            start=i if d else None
        direction=d
    return blocks


def flank_identity(a,b,start,end,flank=30):
    """Check callable identity on both sides of a putative gap event."""
    if start<flank or end+flank>len(a):return None
    values=[]
    for lo,hi in [(start-flank,start),(end,end+flank)]:
        pairs=[(x,y) for x,y in zip(a[lo:hi],b[lo:hi]) if x in DNA and y in DNA]
        if len(pairs)<flank*.9:return None
        values.append(sum(x==y for x,y in pairs)/len(pairs))
    return min(values)


def validate(c):
    """Record local synteny and gap evidence without assigning insertion direction or cause."""
    out=path(c,'results');selection=read_json(out/'orthogroups/selection.json')
    manifest={r['isolate']:r for r in rows(out/'qc/isolate_manifest.tsv')}
    # Evidence for all selected loci is retained even when not nominated for detailed review.
    selected=rows(out/'orthogroups/selected_genes.tsv')
    coords=rows(out/'upstream/upstream_coordinates.tsv')
    repeated={(r['orthogroup'],r['isolate'],int(r['window'])):r for r in rows(out/'repeats/upstream_repeats.tsv')}
    evidence=[];events=[]
    status={(r['orthogroup'],r['region']):r['alignment_status'] for r in rows(out/'upstream/alignment_status.tsv')}
    for og in selection['orthogroups']:
        for w in c['upstream_windows_bp']:
            state=status[(og,f'upstream_{w}')];n_events=flank_supported=0
            if not state.startswith('unresolved'):
                aligned=fasta(out/'upstream/alignments'/f'{file_id(og)}.upstream_{w}.fa')
                for a,b in itertools.combinations(selection['isolates'],2):
                    for start,end,direction in gap_blocks(aligned[a],aligned[b],c['large_indel_bp']):
                        ident=flank_identity(aligned[a],aligned[b],start,end)
                        supported=ident is not None and ident>=.9
                        n_events+=1;flank_supported+=supported
                        events.append({'orthogroup':og,'window':w,'isolate_a':a,'isolate_b':b,
                            'alignment_start':start,'alignment_end':end,'gap_in':a if direction=='a_gap' else b,
                            'indel_bp':end-start,'flanking_identity_min':ident,
                            'flanking_supported':int(supported),
                            'repeat_association':'not_event_validated',
                            'interpretation':'assembly_supported_indel' if supported else 'uncertain_alignment_indel',
                            'ancestral_direction':'unknown'})
            cr=[r for r in coords if r['orthogroup']==og and int(r['window'])==w]
            clean=all(r['extractable']=='1' for r in cr)
            evidence.append({'orthogroup':og,'window':w,'synteny_status':'conserved_immediate_core_neighbors',
                'extraction_complete':int(clean),'large_pairwise_gap_blocks':n_events,
                'flank_supported_gap_blocks':flank_supported,'graph_validation':'not_run',
                'TE_insertion_polymorphism':'not_established',
                'assembly_evidence':'moderate' if clean and not state.startswith('unresolved') else 'uncertain',
                'manual_review':'required','reason':'sequence and synteny screen; ancestral direction and TE events unvalidated'})
    outputs=[table(out/'synteny/locus_evidence.tsv',evidence),
        table(out/'synteny/structural_events.tsv',events,['orthogroup','window','isolate_a','isolate_b',
            'alignment_start','alignment_end','gap_in','indel_bp','flanking_identity_min','flanking_supported',
            'repeat_association','interpretation','ancestral_direction'])]
    return outputs
