#!/usr/bin/env python3
"""Explicit GTF-to-GFF3 conversion retaining source gene boundaries when present."""
import argparse
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input');p.add_argument('output');a=p.parse_args()
    genes={};transcripts={};cds=[]
    with open(a.input) as f:
        for line in f:
            if line.startswith('#') or not line.strip():continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9:raise ValueError('Expected nine-column GTF')
            attr=dict(re.findall(r'(\S+)\s+"([^"]*)"',fields[8]))
            if 'gene_id' not in attr:raise ValueError('GTF feature lacks gene_id')
            gid=attr['gene_id'];chrom=fields[0];start,end=int(fields[3]),int(fields[4]);strand=fields[6]
            if gid in genes and (genes[gid][0],genes[gid][3])!=(chrom,strand):raise ValueError('Inconsistent GTF gene location')
            if gid not in genes:genes[gid]=[chrom,start,end,strand,False]
            if fields[2]=='gene':genes[gid]=[chrom,start,end,strand,True]
            elif not genes[gid][4]:
                genes[gid][1]=min(genes[gid][1],start);genes[gid][2]=max(genes[gid][2],end)
            if 'transcript_id' in attr:
                tx=attr['transcript_id']
                if tx not in transcripts:transcripts[tx]=[gid,chrom,start,end,strand]
                elif transcripts[tx][0]!=gid:raise ValueError('Transcript assigned to multiple genes')
                transcripts[tx][2]=min(transcripts[tx][2],start);transcripts[tx][3]=max(transcripts[tx][3],end)
                if fields[2]=='CDS':cds.append((fields,tx))
            elif fields[2]=='CDS':raise ValueError('CDS lacks transcript_id')
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True);tmp=dest.with_suffix('.part')
    with open(tmp,'w') as f:
        f.write('##gff-version 3\n')
        for gid,(chrom,start,end,strand,original) in sorted(genes.items()):
            origin='source_gene_feature' if original else 'derived_from_GTF_feature_extents'
            f.write(f'{chrom}\tGTF_conversion\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={quote(gid)};boundary_origin={origin}\n')
        for tx,(gid,chrom,start,end,strand) in sorted(transcripts.items()):
            f.write(f'{chrom}\tGTF_conversion\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={quote(tx)};Parent={quote(gid)}\n')
        for fields,tx in cds:
            fields[8]='Parent='+quote(tx);f.write('\t'.join(fields)+'\n')
    tmp.replace(dest)
    print(f'Converted {len(genes)} genes; inspect derived boundary_origin annotations before upstream analysis')

if __name__=='__main__':main()
