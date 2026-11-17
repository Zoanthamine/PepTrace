import hashlib,json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from check_release import check

class PortableReleaseTests(unittest.TestCase):
    def test_entry_works_from_an_unrelated_working_directory(self):
        work=ROOT/'validation/work';work.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=work) as folder:
            p=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/diagnose.py'),'--stats'],cwd=folder,capture_output=True,text=True,encoding='utf-8')
        self.assertEqual(p.returncode,0,p.stderr);stats=json.loads(p.stdout)
        self.assertEqual(stats['counts'],{'rules':450,'branches':308,'appendix':82})
    def test_integrity_guard_detects_changed_or_missing_files(self):
        work=ROOT/'validation/work';work.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=work) as folder:
            r=Path(folder);p=r/'core.py';p.write_text('accepted',encoding='utf-8')
            manifest={'version':'test','library_hash':'test','files_sha256':{'core.py':hashlib.sha256(p.read_bytes()).hexdigest()}}
            (r/'release_manifest.json').write_text(json.dumps(manifest),encoding='utf-8');self.assertTrue(check(r)['pass'])
            p.write_text('changed',encoding='utf-8')
            with self.assertRaises(ValueError):check(r)
            p.unlink()
            with self.assertRaises(ValueError):check(r)
    def test_frozen_package_has_no_absolute_runtime_dependencies(self):
        manifest=json.loads((ROOT/'release_manifest.json').read_text(encoding='utf-8'))
        self.assertTrue(all(not Path(name).is_absolute() and ':' not in name for name in manifest['files_sha256']))
        runtime=json.loads((ROOT/'engine/runtime/manifest.json').read_text(encoding='utf-8'))
        self.assertEqual(runtime['canonical_root'],'knowledge')
        self.assertTrue(runtime['original_source_recheck'])
        self.assertTrue(all('original_not_bundled' in x or x=='hash_verified' for x in runtime['original_source_recheck'].values()))
