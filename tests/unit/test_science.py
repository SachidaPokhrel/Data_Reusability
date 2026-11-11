"""Scientific edge cases independent of downloaded data or third-party packages."""
import itertools
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from pe_upstream.core import Gene,load_genes,coding_sequence,revcomp,load_orthology,write_tsv
from pe_upstream.sequences import extract_one,translate
from pe_upstream.variation import metrics
from pe_upstream.validation import gap_blocks,flank_identity
from pe_upstream.repeats import repeat_records,indices
from pe_upstream.selection import farthest_isolates
from pe_upstream.graph import load_gfa


class ScientificTests(unittest.TestCase):
    def test_coordinate_conversion_and_transcript_choice(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'genes.gff3'
            p.write_text('##gff-version 3\nc\tx\tgene\t11\t39\t.\t+\t.\tID=g\n'
                'c\tx\tmRNA\t11\t39\t.\t+\t.\tID=short;Parent=g\n'
                'c\tx\tmRNA\t11\t39\t.\t+\t.\tID=long;Parent=g\n'
                'c\tx\tCDS\t11\t13\t.\t+\t0\tParent=short\n'
                'c\tx\tCDS\t11\t16\t.\t+\t0\tParent=long\n'
                'c\tx\tCDS\t34\t39\t.\t+\t0\tParent=long\n')
            g=load_genes(p)['g'];self.assertEqual((g.start,g.end),(10,39));self.assertEqual(g.transcript,'long')
            genome={'c':'N'*10+'ATGGCC'+'N'*17+'GCCTAA'}
            self.assertEqual(coding_sequence(g,genome),'ATGGCCGCCTAA')
            self.assertEqual(translate(coding_sequence(g,genome)),'MAA')

    def test_negative_multiexon(self):
        from pe_upstream.core import Feature
        coding='ATGGCCTAA';genome={'c':revcomp(coding[6:])+'N'*4+revcomp(coding[:6])}
        g=Gene('g','c',0,13,'-',[Feature('c','CDS',0,3,'-','0',{}),Feature('c','CDS',7,13,'-','0',{})])
        self.assertEqual(coding_sequence(g,genome),coding)

    def test_missing_upstream_and_edges(self):
        genome={'c':'A'*50};g=Gene('g','c',0,10,'+',[])
        seq,info=extract_one(g,genome,{'g':g},20)
        self.assertEqual(seq,'');self.assertEqual(info['actual_length'],0)
        self.assertIn('contig_boundary',info['flags']);self.assertIn('short_upstream',info['flags'])

    def test_upstream_neighbor_clip(self):
        g=Gene('g','c',40,50,'+',[]);n=Gene('n','c',10,30,'+',[])
        seq,info=extract_one(g,{'c':'A'*100},{'g':g,'n':n},20)
        self.assertEqual(len(seq),10);self.assertEqual(info['upstream_gene_distance'],10)
        self.assertEqual(info['overlaps_neighbor'],1)

    def test_gap_and_ambiguity(self):
        g=Gene('g','c',20,40,'-',[])
        seq,info=extract_one(g,{'c':'A'*40+'NNNNNA'*10},{'g':g},10)
        self.assertEqual(info['ambiguous_bases'],9)
        self.assertGreater(info['gap_runs'],0)

    def test_indel_events_not_bases(self):
        m=metrics('AAAA-----CCCC--TTTT','AAAAGGGGGCCCCAATTTT',large_indel=4)
        self.assertEqual(m['indel_events'],2);self.assertEqual(m['indel_bases'],7)
        self.assertEqual(m['large_indel_events'],1)
        self.assertLess(m['callable_fraction'],1)

    def test_low_coverage_not_conserved(self):
        m=metrics('A'*10+'-'*90,'A'*100)
        self.assertEqual(m['identity'],1)
        self.assertEqual(m['callable_fraction'],.1)

    def test_real_te_insertion_fixture_requires_flanks(self):
        # Constructed truth: one sequence contains an internal 60-bp insertion.
        a='A'*40+'-'*60+'C'*40;b='A'*40+'G'*60+'C'*40
        blocks=gap_blocks(a,b,50);self.assertEqual(blocks,[(40,100,'a_gap')])
        self.assertEqual(flank_identity(a,b,40,100),1)
        self.assertIsNone(flank_identity(a,b,5,65))
        # Gap evidence alone does not assign TE status in the implementation.

    def test_repeat_union_and_te_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'r.gff3';p.write_text('c\tx\trepeat_region\t1\t20\t.\t+\t.\tClassification=LTR/Gypsy\n'
                'c\tx\trepeat_region\t10\t30\t.\t+\t.\tClassification=Simple_repeat\n')
            rr=repeat_records(p,['LTR']);self.assertTrue(rr[0]['is_TE']);self.assertFalse(rr[1]['is_TE'])
            self.assertEqual(indices(rr).overlap('c',0,30),30)
            self.assertEqual(indices(rr,True).overlap('c',0,30),20)

    def test_missing_and_unequal_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'o.tsv'
            data=[dict(orthogroup='core',isolate=s,gene_id=s+'1') for s in ['a','b','c']]
            data += [dict(orthogroup='missing',isolate='a',gene_id='a2'),
                dict(orthogroup='dup',isolate='a',gene_id='a3'),dict(orthogroup='dup',isolate='a',gene_id='a4')]
            write_tsv(p,data,['orthogroup','isolate','gene_id'])
            self.assertEqual(set(load_orthology(p,['a','b','c'])),{'core'})

    def test_genotype_diversity_selection(self):
        isolates=['a','b','c','d'];dist={k:abs(isolates.index(k[0])-isolates.index(k[1]))/3 for k in itertools.combinations(isolates,2)}
        self.assertEqual(farthest_isolates(isolates,dist,2),['a','d'])

    def test_graph_walk_and_missing_segments(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.gfa';p.write_text('S\t1\tACGT\nS\t2\tGGG\nP\tp\t1+,2-\t*\n')
            segments,paths=load_gfa(p);self.assertEqual(paths['p'],[('1','+'),('2','-')])
            p.write_text('S\t1\t*\nP\tp\t1+\t*\n')
            with self.assertRaisesRegex(ValueError,'absent'):load_gfa(p)

if __name__=='__main__':unittest.main()
