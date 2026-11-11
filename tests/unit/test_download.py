import hashlib
import io
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from pe_upstream.acquisition import download_one,safe_extract

class Response(io.BytesIO):
    def __init__(self,data,status=200,headers=None):
        super().__init__(data);self.status=status;self.headers=headers or {}

class DownloadTests(unittest.TestCase):
    def test_checksum_atomic_and_skip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'file';data=b'verified';h='sha256:'+hashlib.sha256(data).hexdigest()
            with patch('urllib.request.urlopen',return_value=Response(data)) as url:
                self.assertEqual(download_one('https://example.org/file',p,h,len(data),1),'verified_downloaded')
                self.assertEqual(download_one('https://example.org/file',p,h,len(data),1),'verified_existing')
                self.assertEqual(url.call_count,1)
            self.assertFalse(p.with_name('file.part').exists())

    def test_checksum_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'file'
            with patch('urllib.request.urlopen',return_value=Response(b'bad')):
                with self.assertRaisesRegex(ValueError,'Checksum mismatch'):
                    download_one('https://example.org/file',p,'sha256:'+'0'*64,3,1)
            self.assertFalse(p.exists());self.assertFalse(p.with_name('file.part').exists())

    def test_incomplete_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'file';data=b'abcdef';h='sha256:'+hashlib.sha256(data).hexdigest()
            with patch('urllib.request.urlopen',return_value=Response(b'abc')):
                with self.assertRaisesRegex(ValueError,'Incomplete'):download_one('https://example.org/file',p,h,6,1)
            self.assertFalse(p.exists());self.assertEqual(p.with_name('file.part').read_bytes(),b'abc')
            with patch('urllib.request.urlopen',return_value=Response(b'def',206,{'Content-Range':'bytes 3-5/6'})) as url:
                download_one('https://example.org/file',p,h,6,1)
                self.assertEqual(url.call_args[0][0].get_header('Range'),'bytes=3-')
            self.assertEqual(p.read_bytes(),data)

    def test_ignored_range_restarts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'file';p.with_name('file.part').write_bytes(b'abc');data=b'abcdef'
            with patch('urllib.request.urlopen',return_value=Response(data)):
                download_one('https://example.org/file',p,'md5:'+hashlib.md5(data).hexdigest(),6,1)
            self.assertEqual(p.read_bytes(),data)

    def test_zip_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'a.zip'
            with zipfile.ZipFile(p,'w') as z:z.writestr('../escape','bad')
            with self.assertRaisesRegex(ValueError,'Unsafe'):safe_extract(p,Path(tmp)/'out')
            self.assertFalse((Path(tmp)/'escape').exists())
