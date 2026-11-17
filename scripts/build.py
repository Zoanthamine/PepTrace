"""Rebuild all runtime/index files from package-local canonical JSON, then freeze hashes."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
from check_release import digest
ROOT=Path(__file__).resolve().parents[1];ENGINE=ROOT/'engine'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def freeze():
    def accepted(p):
        rel=p.relative_to(ROOT)
        if not p.is_file() or p.name in {'release_manifest.json'} or '__pycache__' in rel.parts:return False
        if rel.parts[0]=='engine':
            if 'rebuild' in rel.parts:
                return rel.parts[:3]==('engine','rebuild','holdouts') and p.name in {'registration.json','answers.json','question-1.json','question-2.json','question-3.json'}
            return p.suffix in {'.py','.json','.sqlite3'}
        if rel.parts[0] in {'scripts','references','agents'}:return p.suffix in {'.py','.json','.md','.yaml'}
        if rel.parts[0] in {'examples','tests','.github'}:return p.suffix in {'.py','.json','.md','.yaml','.yml'}
        return len(rel.parts)==1 and p.name in {'SKILL.md','README.md','LICENSE','requirements.txt','requirements-msms.txt','.gitattributes','.gitignore'}
    paths=sorted(p for p in ROOT.rglob('*') if accepted(p))
    runtime=json.loads((ENGINE/'runtime/manifest.json').read_text(encoding='utf-8'))
    manifest={'version':'1.0.0','engine':'codex-dsh-fused-1','schema':'peptrace-portable-release-1',
        'library_hash':runtime['logical_sha256'],'counts':runtime['counts'],
        'files_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        'original_book_archive':'separate GitHub Release asset, not a runtime dependency'}
    (ROOT/'release_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    return manifest
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source-root',type=Path);a=parser.parse_args()
    (ENGINE/'rebuild').mkdir(exist_ok=True)
    args=['--source-root',str(a.source_root.resolve())] if a.source_root else []
    for script,extra in [('build_knowledge.py',args),('build_knowledge_index.py',[])]:
        subprocess.run([sys.executable,'-B',str(ENGINE/'tools'/script),*extra],cwd=ENGINE,check=True)
    manifest=freeze();print(json.dumps({'version':manifest['version'],'library_hash':manifest['library_hash'],
        'files':len(manifest['files_sha256']),'counts':manifest['counts']},ensure_ascii=False))
if __name__=='__main__':main()
