import json
import tempfile
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from pe_upstream.common import ROOT,load_config
from pe_upstream.quality import normalize_orthology
from pe_upstream.core import load_genes

class ContractTests(unittest.TestCase):
    def test_synthetic_aligner_forbidden_for_real_mode(self):
        c=json.loads((ROOT/'config/pilot.yaml').read_text());c['test_stub_alignment']=True
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'c.yaml';p.write_text(json.dumps(c))
            with self.assertRaisesRegex(ValueError,'forbidden'):load_config(p)

    def test_unknown_gene_id_fails_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);(tmp/'raw').mkdir()
            (tmp/'raw/Pe_orthogroups.tsv').write_text('orthogroup\tisolate\tgene_id\nOG1\tA\tunknown\n')
            gff=tmp/'a.gff3';gff.write_text('c\tx\tgene\t1\t3\t.\t+\t.\tID=known\n')
            c=load_config(ROOT/'config/test.yaml');c['paths']['data']=str(tmp)
            with self.assertRaisesRegex(ValueError,'Unmapped gene IDs'):
                normalize_orthology(c,[{'isolate':'A','gff':str(gff)}])

    def test_unrecognized_schema_does_not_guess(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);(tmp/'raw').mkdir()
            (tmp/'raw/Pe_orthogroups.tsv').write_text('group\tannotation\tcount\nOG1\tunknown\t2\n')
            c=load_config(ROOT/'config/test.yaml');c['paths']['data']=str(tmp)
            with self.assertRaisesRegex(ValueError,'Unrecognized orthology schema'):
                normalize_orthology(c,[{'isolate':'A','gff':'unused'}])

    def test_duplicate_transcript_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'a.gff3';p.write_text('c\tx\tgene\t1\t9\t.\t+\t.\tID=g\n'
                'c\tx\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n'
                'c\tx\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n')
            with self.assertRaisesRegex(ValueError,'Duplicate transcript'):load_genes(p)
