"""Read-only ZIP identity check against the book source catalog; never extracts files."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from check_release import digest
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('archive',type=Path);parser.add_argument('--output',type=Path);a=parser.parse_args()
    expected=json.loads((ROOT/'references/source-archive.json').read_text(encoding='utf-8'))
    catalog=json.loads((ROOT/'references/source-catalog.json').read_text(encoding='utf-8'))
    actual=digest(a.archive)
    if actual!=expected['sha256'] or a.archive.stat().st_size!=expected['bytes']:raise SystemExit('Original archive differs from registered bytes/SHA-256')
    rows=[]
    with zipfile.ZipFile(a.archive) as z:
        members={Path(name).name:name for name in z.namelist() if not name.endswith('/')}
        for source in catalog['source_files']:
            if source['title'] not in members:raise SystemExit('Missing registered source: '+source['title'])
            sha=hashlib.sha256(z.read(members[source['title']])).hexdigest()
            if sha!=source['sha256']:raise SystemExit('Source PDF hash mismatch: '+source['doc_id'])
            rows.append({'doc_id':source['doc_id'],'filename':source['title'],'sha256':sha,'pages':source['pages'],'pass':True})
    result={'archive':a.archive.name,'bytes':expected['bytes'],'sha256':actual,'verified_source_pdfs':len(rows),
        'source_pages':sum(r['pages'] for r in rows),'rows':rows,'pass':True}
    if a.output:a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False))
if __name__=='__main__':main()
