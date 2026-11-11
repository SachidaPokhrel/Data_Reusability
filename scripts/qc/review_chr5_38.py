#!/usr/bin/env python3
"""Rebuild or inspect the diagnostic CDS-to-neighbor alignment for Chr5_38."""
import argparse
import csv
from itertools import combinations
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from pe_upstream.common import ROOT, table
from pe_upstream.core import fasta, load_genes, revcomp, write_fasta


def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT / 'results/pilot')
    parser.add_argument('--rebuild', action='store_true', help='Extract sequences and run MAFFT')
    parser.add_argument('--mafft', default='mafft')
    args = parser.parse_args()
    result = args.results
    review = result / 'candidates/Chr5_38_review'
    windows = {r['isolate']: r for r in read_tsv(result / 'upstream/upstream_coordinates.tsv')
               if r['orthogroup'] == 'Chr5_38' and r['window'] == '500'}
    selected = {r['isolate']: r for r in read_tsv(result / 'orthogroups/selected_genes.tsv')
                if r['orthogroup'] == 'Chr5_38'}
    if len(windows) != 6 or set(windows) != set(selected):
        raise ValueError('Expected the six pilot isolates for Chr5_38')
    if args.rebuild:
        review.mkdir(parents=True, exist_ok=True)
        samples = {r['isolate']: r for r in read_tsv(ROOT / 'config/isolate_files.zenodo_v7.tsv')}
        sequences = {}
        for isolate, r in selected.items():
            sample = samples[isolate]
            genes = load_genes(ROOT / sample['gff'])
            gene = genes[r['gene_id']]
            if gene.strand != '-':
                raise ValueError(f'Unexpected strand for {isolate}')
            first = max(gene.cds, key=lambda segment: segment.end)
            flank = min(60, first.end - first.start)
            neighbor = min((n for n in genes.values()
                            if n.chrom == gene.chrom and n.start >= gene.end),
                           key=lambda n: (n.start, n.id))
            genome = fasta(ROOT / sample['fasta'])
            sequences[isolate] = revcomp(genome[gene.chrom][first.end-flank:neighbor.start])
        source = review / 'CDS_to_neighbor.fasta'
        write_fasta(source, sequences)
        tmp = review / 'CDS_to_neighbor.aligned.part'
        with tmp.open('w') as handle:
            subprocess.run([args.mafft, '--auto', str(source)], stdout=handle, check=True)
        tmp.replace(review / 'CDS_to_neighbor.aligned.fasta')

    aligned = fasta(review / 'CDS_to_neighbor.aligned.fasta')
    if set(aligned) != set(windows) or len({len(s) for s in aligned.values()}) != 1:
        raise ValueError('Alignment does not match the coordinate table')
    columns, projections = {}, []
    for isolate, seq in aligned.items():
        r = windows[isolate]
        boundary = int(r['gene_end']) + int(r['upstream_gene_distance'])
        start, end = boundary - int(r['end']), boundary - int(r['start'])
        available = [i for i, base in enumerate(seq, 1) if base != '-']
        if r['strand'] != '-' or not 0 <= start < end <= len(available) or end-start != 500:
            raise ValueError(f'Invalid projection for {isolate}')
        span = available[start:end]
        columns[isolate] = set(span)
        projections.append(dict(isolate=isolate, first_alignment_column=span[0],
                                last_alignment_column=span[-1]))
    pairs = []
    for a, b in combinations(aligned, 2):
        if (a in {'Pe5', 'Pe14'}) == (b in {'Pe5', 'Pe14'}):
            continue
        calls = [(x, y) for x, y in zip(aligned[a], aligned[b]) if x in 'ACGT' and y in 'ACGT']
        pairs.append(dict(isolate_a=a, isolate_b=b,
                          shared_original_window_columns=len(columns[a] & columns[b]),
                          callable_bases=len(calls),
                          larger_interval_callable_identity=sum(x == y for x, y in calls)/len(calls)))
    table(review / 'window_projection.tsv', projections)
    table(review / 'between_group_identity.tsv', pairs)
    print(f'Projected {len(projections)} windows; wrote {len(pairs)} between-group comparisons.')


if __name__ == '__main__':
    main()
