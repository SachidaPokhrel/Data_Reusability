"""Core utilities. Coordinates internally are zero-based, half-open."""
from __future__ import annotations
import bisect
import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

DNA = frozenset('ACGT')

def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def fasta(path):
    result = {}
    name, parts = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    result[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
                if name in result:
                    raise ValueError(f'Duplicate FASTA identifier: {name}')
            else:
                if name is None:
                    raise ValueError(f'Sequence before FASTA header in {path}')
                parts.append(line)
    if name is not None:
        if name in result:
            raise ValueError(f'Duplicate FASTA identifier: {name}')
        result[name] = ''.join(parts).upper()
    return result

def write_fasta(path, sequences):
    with open(path, 'w') as handle:
        for name, seq in sequences.items():
            handle.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + '\n')

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn',
                                     'TGCAYRKMVHDBNtgcayrkmvhdbn'))[::-1]

@dataclass
class Feature:
    chrom: str
    kind: str
    start: int
    end: int
    strand: str
    phase: str
    attrs: dict


def features(path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Expected 9-column GFF3 in {path}: {line[:80]}')
            attrs = {}
            for item in fields[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k] = unquote(v)
            a, b = int(fields[3])-1, int(fields[4])
            if a < 0 or b <= a:
                raise ValueError(f'Invalid GFF coordinates in {path}')
            yield Feature(fields[0], fields[2], a, b, fields[6], fields[7], attrs)

@dataclass
class Gene:
    id: str
    chrom: str
    start: int
    end: int
    strand: str
    cds: list
    transcript: str = ""


def load_genes(path):
    rows = list(features(path))
    genes, tx_to_gene, segments = {}, {}, defaultdict(list)
    for f in rows:
        if f.kind == 'gene':
            ident = f.attrs.get('ID')
            if not ident or ident in genes:
                raise ValueError(f'Missing/duplicate gene ID in {path}: {ident}')
            genes[ident] = Gene(ident, f.chrom, f.start, f.end, f.strand, [])
    for f in rows:
        if f.kind in ('mRNA', 'transcript'):
            parent = f.attrs.get('Parent', '')
            if parent in genes:
                if f.attrs['ID'] in tx_to_gene:
                    raise ValueError(f'Duplicate transcript ID: {f.attrs["ID"]}')
                tx_to_gene[f.attrs['ID']] = parent
    for f in rows:
        if f.kind == 'CDS':
            for parent in f.attrs.get('Parent', '').split(','):
                if parent in tx_to_gene or parent in genes:
                    segments[parent].append(f)
    # One deterministic representative transcript: greatest total CDS length, then ID.
    candidates = defaultdict(list)
    for tx, cds in segments.items():
        gid = tx_to_gene.get(tx, tx)
        candidates[gid].append((sum(x.end-x.start for x in cds), tx, cds))
    for gid, choices in candidates.items():
        choice = sorted(choices, key=lambda x: (-x[0], x[1]))[0]
        genes[gid].cds = choice[2]
        genes[gid].transcript = choice[1]
    if not genes:
        raise ValueError(f'No gene features found in {path}; convert to GFF3 with gene/CDS hierarchy')
    return genes


def coding_sequence(gene, genome):
    if not gene.cds or gene.strand not in ('+', '-'):
        raise ValueError('missing_cds_or_strand')
    ordered = sorted(gene.cds, key=lambda x: x.start, reverse=gene.strand == '-')
    if ordered[0].phase != '0':
        raise ValueError('partial_cds')
    accum, parts = 0, []
    previous = None
    for f in ordered:
        if (f.chrom != gene.chrom or f.strand != gene.strand or
            f.start < gene.start or f.end > gene.end or f.end > len(genome[gene.chrom])):
            raise ValueError('inconsistent_cds_coordinates')
        if previous and max(previous.start, f.start) < min(previous.end, f.end):
            raise ValueError('overlapping_cds_segments')
        expected_phase = (3 - accum % 3) % 3
        if f.phase != str(expected_phase):
            raise ValueError('inconsistent_cds_phase')
        seq = genome[f.chrom][f.start:f.end]
        parts.append(revcomp(seq) if gene.strand == '-' else seq)
        accum += len(seq)
        previous = f
    seq = ''.join(parts)
    if len(seq) % 3:
        raise ValueError('incomplete_cds')
    if any(seq[i:i+3] in {'TAA', 'TAG', 'TGA'} for i in range(0, len(seq)-3, 3)):
        raise ValueError('internal_stop')
    return seq


class IntervalIndex:
    def __init__(self, intervals):
        self.by_chrom = {}
        for chrom, entries in intervals.items():
            merged = []
            for a, b in sorted(entries):
                if merged and a <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
                else:
                    merged.append((a, b))
            self.by_chrom[chrom] = (merged, [a for a, b in merged], [b for a, b in merged])

    def overlap(self, chrom, start, end):
        merged, starts, ends = self.by_chrom.get(chrom, ([], [], []))
        i = bisect.bisect_right(ends, start)
        total = 0
        while i < len(merged) and merged[i][0] < end:
            a, b = merged[i]
            total += max(0, min(b, end)-max(a, start))
            i += 1
        return total

    def distance(self, chrom, start, end):
        merged, starts, ends = self.by_chrom.get(chrom, ([], [], []))
        if not merged:
            return None
        if self.overlap(chrom, start, end):
            return 0
        i = bisect.bisect_left(starts, end)
        distances = []
        if i < len(merged):
            distances.append(merged[i][0]-end)
        if i:
            distances.append(start-merged[i-1][1])
        return max(0, min(distances))


def gene_index(genes):
    intervals = defaultdict(list)
    for g in genes.values():
        intervals[g.chrom].append((g.start, g.end))
    return IntervalIndex(intervals)


def load_repeats(path):
    """Repeat GFF3; all supplied repeat features, not necessarily all TEs."""
    intervals = defaultdict(list)
    n = 0
    for f in features(path):
        intervals[f.chrom].append((f.start, f.end))
        n += 1
    if not n:
        raise ValueError(f'No repeat features in {path}')
    return IntervalIndex(intervals)


def upstream(gene, genome, window, other_genes):
    if gene.strand == '+':
        start, end = gene.start-window, gene.start
    elif gene.strand == '-':
        start, end = gene.end, gene.end+window
    else:
        raise ValueError('unknown_strand')
    if start < 0 or end > len(genome[gene.chrom]):
        raise ValueError('upstream_contig_edge')
    if other_genes.overlap(gene.chrom, start, end):
        raise ValueError('upstream_overlaps_gene')
    seq = genome[gene.chrom][start:end]
    return (revcomp(seq) if gene.strand == '-' else seq), start, end


def synteny_signatures(genes, gene_to_og):
    """Immediate neighbors must be mapped single-copy core genes.

    Signatures use transcription-relative order and relative strand, allowing
    complete locus inversion. Missing or unmapped adjacent genes fail closed.
    """
    by_chrom = defaultdict(list)
    for gene in genes.values():
        by_chrom[gene.chrom].append(gene)
    result = {}
    for entries in by_chrom.values():
        entries.sort(key=lambda g: (g.start, g.end, g.id))
        for i in range(1, len(entries)-1):
            g = entries[i]
            left, right = entries[i-1], entries[i+1]
            if g.strand not in ('+', '-') or any(n.strand not in ('+', '-') for n in (left,right)):
                continue
            if left.end > g.start or g.end > right.start:
                continue
            if left.id not in gene_to_og or right.id not in gene_to_og:
                continue
            neighbors = (left, right) if g.strand == '+' else (right, left)
            result[g.id] = tuple((gene_to_og[n.id], n.strand == g.strand) for n in neighbors)
    return result


def pair_metrics(a, b):
    if len(a) != len(b):
        raise ValueError('Alignment rows have unequal length')
    callable_sites = mismatches = gap_bases = comparable = 0
    for x, y in zip(a.upper(), b.upper()):
        if x == '-' and y == '-':
            continue
        if x in DNA and y in DNA:
            callable_sites += 1
            mismatches += x != y
            comparable += 1
        elif (x == '-' and y in DNA) or (y == '-' and x in DNA):
            gap_bases += 1
            comparable += 1
    ungapped_max = max(sum(x in DNA for x in a.upper()), sum(x in DNA for x in b.upper()))
    return {
        'callable_bases': callable_sites,
        'callable_fraction': callable_sites / ungapped_max if ungapped_max else 0,
        'snp_divergence': mismatches / callable_sites if callable_sites else None,
        'gap_fraction': gap_bases / comparable if comparable else None,
        'combined_divergence': (mismatches+gap_bases)/comparable if comparable else None,
    }


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def load_orthology(path, isolates):
    by_og = defaultdict(lambda: defaultdict(list))
    seen = set()
    for row in read_tsv(path):
        og, isolate, gene = (row[k] for k in ('orthogroup', 'isolate', 'gene_id'))
        if not all((og, isolate, gene)):
            raise ValueError('Blank orthology identifier')
        key = (isolate, gene)
        if key in seen:
            raise ValueError(f'Gene assigned more than once: {key}')
        seen.add(key)
        by_og[og][isolate].append(gene)
    return {og: {s: members[s][0] for s in isolates}
            for og, members in by_og.items()
            if all(len(members[s]) == 1 for s in isolates)}
