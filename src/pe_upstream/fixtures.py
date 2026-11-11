"""Deterministic synthetic data, exclusively under tests/work."""
import itertools
import random
from .common import ROOT,table,write_json,path
from .core import write_fasta,revcomp


def create(c):
    """Generate deterministic TEST inputs with known variation and extraction failures."""
    if c['mode']!='test':raise ValueError('Fixture generation is TEST-only')
    root=ROOT/'tests/work';root.mkdir(parents=True,exist_ok=True)
    rng=random.Random(2026);base={}
    for chrom in ['chr1','chr2','chr3']:
        base[chrom]=''.join(rng.choice('ACGT') for _ in range(3500))
    samples=[];orth=[];presence=[]
    coding='ATG'+'GCT'*18+'TAA'
    for index in range(6):
        s=f'S{index+1}';genome={k:list(v) for k,v in base.items()};gff=['##gff-version 3'];repeat=['##gff-version 3']
        for ci,chrom in enumerate(base):
            for gi in range(8):
                start=200+gi*400;end=start+len(coding);strand='-' if gi%2 else '+'
                gid=f'{s}_g{ci}_{gi}';tx=gid+'.t1';og=f'OG{ci}_{gi}'
                genome[chrom][start:end]=coding if strand=='+' else revcomp(coding)
                # Intentionally variable upstream sites for an internal locus.
                if ci==1 and gi==3 and index>=3:
                    for pos in range(end+20,end+40):
                        old=genome[chrom][pos];genome[chrom][pos]=next(x for x in 'ACGT' if x!=old)
                if ci==2 and gi==4 and index==5:
                    genome[chrom][start-20:start-10]=list('N'*10)
                # Multi-exon first gene, preserving the exact coding sequence.
                if gi==0:
                    end=start+len(coding)+20
                    seq=coding[:30]+'G'*20+coding[30:]
                    genome[chrom][start:end]=seq
                    segments=[(start,start+30),(start+50,end)]
                else:segments=[(start,end)]
                gff += [f'{chrom}\tsynthetic\tgene\t{start+1}\t{end}\t.\t{strand}\t.\tID={gid}',
                    f'{chrom}\tsynthetic\tmRNA\t{start+1}\t{end}\t.\t{strand}\t.\tID={tx};Parent={gid}']
                for a,b in segments:
                    gff.append(f'{chrom}\tsynthetic\tCDS\t{a+1}\t{b}\t.\t{strand}\t0\tParent={tx}')
                orth.append({'orthogroup':og,'isolate':s,'gene_id':gid})
                repeat.append(f'{chrom}\tsynthetic\trepeat_region\t{max(1,start-60)}\t{start-30}\t.\t+\t.\tID=r{ci}_{gi};Name=synthetic_LTR;Classification=LTR/Gypsy')
        fa=root/f'{s}.fa';gf=root/f'{s}.gff3';rp=root/f'{s}.repeats.gff3'
        write_fasta(fa,{k:''.join(v) for k,v in genome.items()});gf.write_text('\n'.join(gff)+'\n');rp.write_text('\n'.join(repeat)+'\n')
        samples.append({'isolate':s,'haplotype':'synthetic_consensus','assembly_role':'synthetic',
            'fasta':str(fa),'gff':str(gf),'repeats':str(rp)})
    table(ROOT/c['samples'],samples)
    raw=path(c,'data')/'raw';raw.mkdir(parents=True,exist_ok=True)
    table(raw/'Pe_orthogroups.tsv',orth)
    table(root/'distances.tsv',[{'isolate_a':a['isolate'],'isolate_b':b['isolate'],'distance':.01*(j-i)}
        for (i,a),(j,b) in itertools.combinations(enumerate(samples),2)])
    return [ROOT/c['samples'],raw/'Pe_orthogroups.tsv',root/'distances.tsv']
