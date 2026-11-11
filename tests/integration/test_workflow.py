import json
import csv
import os
import shutil
import subprocess
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

class WorkflowIntegration(unittest.TestCase):
    def test_offline_full_synthetic_stages_and_resume(self):
        command=[sys.executable,'scripts/pe.py','--config','config/test.yaml','test']
        first=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(first.returncode,0,first.stderr+first.stdout)
        master=ROOT/'tests/work/results/test/candidates/master_locus_table.tsv'
        expected=master.read_bytes()
        with master.open() as handle:
            observed=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual({r['orthogroup_id'] for r in observed if r['candidate_status']=='B'},{'OG1_3'})
        self.assertEqual({r['orthogroup_id'] for r in observed if r['candidate_status']=='E'},{'OG2_4'})
        second=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(second.returncode,0,second.stderr+second.stdout)
        self.assertEqual(master.read_bytes(),expected)
        report=ROOT/'tests/work/reports/test/test_report.md'
        self.assertIn('SYNTHETIC DATA',report.read_text())
        audit=json.loads((ROOT/'tests/work/reports/test/provenance/audit.json').read_text())
        self.assertEqual(audit['status'],'pass')
        # Corrupt an output: the audit must fail instead of accepting stale results.
        original=master.read_bytes();master.write_text(master.read_text()+'\n')
        checked=subprocess.run([sys.executable,'scripts/pe.py','--config','config/test.yaml','audit'],
            cwd=ROOT,capture_output=True,text=True)
        self.assertNotEqual(checked.returncode,0)
        master.write_bytes(original)
        restored=subprocess.run([sys.executable,'scripts/pe.py','--config','config/test.yaml','audit'],
            cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(restored.returncode,0,restored.stderr)

    @unittest.skipUnless(shutil.which('mafft'),'MAFFT unavailable')
    def test_real_mafft_preserves_indel_sequences(self):
        import tempfile
        sys.path.insert(0,str(ROOT/'src'))
        from pe_upstream.core import fasta,write_fasta
        from pe_upstream.sequences import run_alignment
        with tempfile.TemporaryDirectory() as tmp:
            src=Path(tmp)/'in.fa';dest=Path(tmp)/'out.fa'
            seq={'a':'ACGT'*30,'b':'ACGT'*15+'GGG'+'ACGT'*15}
            write_fasta(src,seq);run_alignment(src,dest,{'test_stub_alignment':False,'threads':1})
            result=fasta(dest)
            self.assertEqual({k:v.replace('-','') for k,v in result.items()},seq)

if __name__=='__main__':unittest.main()
