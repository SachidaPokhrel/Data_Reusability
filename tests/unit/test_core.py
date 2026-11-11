import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from pe_upstream.core import (Gene,Feature,IntervalIndex,revcomp,upstream,gene_index,
    coding_sequence,pair_metrics,load_orthology,write_tsv,synteny_signatures,fasta)

class CoreTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(revcomp('ACGTNRY'),'RYNACGT')

    def test_repeat_union(self):
        ix=IntervalIndex({'c':[(2,8),(4,10),(10,12),(20,25)]})
        self.assertEqual(ix.overlap('c',0,30),15)
        self.assertEqual(ix.overlap('c',5,9),4)
        self.assertEqual(ix.distance('c',14,16),2)
        self.assertIsNone(ix.distance('absent',0,1))

    def test_upstream_strands_and_bounds(self):
        genome={'c':'ACGT'*100}
        genes={'p':Gene('p','c',40,60,'+',[]),'m':Gene('m','c',200,220,'-',[])}
        index=gene_index(genes)
        self.assertEqual(upstream(genes['p'],genome,20,index),(genome['c'][20:40],20,40))
        self.assertEqual(upstream(genes['m'],genome,20,index),(revcomp(genome['c'][220:240]),220,240))
        with self.assertRaisesRegex(ValueError,'contig_edge'):
            upstream(genes['p'],genome,50,index)
        with self.assertRaisesRegex(ValueError,'overlaps_gene'):
            upstream(genes['m'],genome,20,IntervalIndex({'c':[(225,230)]}))

    def test_cds_phase_split_codon(self):
        # ATGG + CCTAA -> ATG GCC TAA. The phase=2 bases complete a split codon.
        segments=[Feature('c','CDS',0,4,'+','0',{}),Feature('c','CDS',8,13,'+','2',{})]
        gene=Gene('g','c',0,13,'+',segments)
        self.assertEqual(coding_sequence(gene,{'c':'ATGGNNNNCCTAA'}),'ATGGCCTAA')
        segments[1].phase='0'
        with self.assertRaisesRegex(ValueError,'phase'):
            coding_sequence(gene,{'c':'ATGGNNNNCCTAA'})

    def test_negative_cds(self):
        gene=Gene('g','c',0,9,'-',[Feature('c','CDS',0,9,'-','0',{})])
        self.assertEqual(coding_sequence(gene,{'c':revcomp('ATGGCCTAA')}),'ATGGCCTAA')

    def test_internal_stop(self):
        g=Gene('g','c',0,9,'+',[Feature('c','CDS',0,9,'+','0',{})])
        with self.assertRaisesRegex(ValueError,'internal_stop'):
            coding_sequence(g,{'c':'ATGTAAGCC'})

    def test_pair_metrics_callable_and_gaps(self):
        m=pair_metrics('ACGT-N-','ATGTAN-')
        self.assertEqual(m['callable_bases'],4)
        self.assertAlmostEqual(m['snp_divergence'],.25)
        self.assertAlmostEqual(m['gap_fraction'],.2)
        self.assertAlmostEqual(m['combined_divergence'],.4)
        self.assertAlmostEqual(m['callable_fraction'],.8)
        self.assertEqual(pair_metrics('NN--','--NN')['callable_fraction'],0)
        with self.assertRaises(ValueError): pair_metrics('A','AA')

    def test_core_and_duplicate_memberships(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'og.tsv'
            rows=[{'orthogroup':og,'isolate':s,'gene_id':og+s}
                  for og in ['core','accessory'] for s in (['a','b'] if og=='core' else ['a'])]
            write_tsv(p,rows,['orthogroup','isolate','gene_id'])
            self.assertEqual(set(load_orthology(p,['a','b'])),{'core'})
            write_tsv(p,rows+[rows[0]],['orthogroup','isolate','gene_id'])
            with self.assertRaisesRegex(ValueError,'more than once'):load_orthology(p,['a','b'])

    def test_synteny_inversion(self):
        gs={f'g{i}':Gene(f'g{i}','c',i*100,i*100+20,'+',[]) for i in range(3)}
        mapping={f'g{i}':f'og{i}' for i in range(3)}
        sig=synteny_signatures(gs,mapping)['g1']
        inverted={f'g{i}':Gene(f'g{i}','c',(2-i)*100,(2-i)*100+20,'-',[]) for i in range(3)}
        self.assertEqual(sig,synteny_signatures(inverted,mapping)['g1'])
        self.assertNotIn('g0',synteny_signatures(gs,mapping))

    def test_duplicate_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'a.fa';p.write_text('>a\nA\n>a\nC\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):fasta(p)

if __name__=='__main__':unittest.main()
