"""Portable content-integrity check; no external project path or credentials."""
import hashlib,json,os,sys
from pathlib import Path
os.environ['PYTHONUTF8']='1'
os.environ['PYTHONIOENCODING']='utf-8'
for stream in (sys.stdout,sys.stderr):
    if hasattr(stream,'reconfigure'):stream.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(root=ROOT):
    manifest=json.loads((root/'release_manifest.json').read_text(encoding='utf-8'))
    failures=[]
    for name,expected in manifest['files_sha256'].items():
        p=root/name
        if not p.is_file() or digest(p)!=expected:failures.append(name)
    if failures:raise ValueError('Release files changed or missing; run scripts/build.py and revalidate: '+', '.join(failures[:8]))
    return {'version':manifest['version'],'checked_files':len(manifest['files_sha256']),
        'library_hash':manifest['library_hash'],'pass':True}
