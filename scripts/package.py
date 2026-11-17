"""Make a deterministic complete skill ZIP and its SHA-256 for a release attachment."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from check_release import check
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path);a=parser.parse_args();check()
    output=a.output.resolve();output.parent.mkdir(parents=True,exist_ok=True)
    files=[]
    for p in ROOT.rglob('*'):
        rel=p.relative_to(ROOT)
        if not p.is_file() or any(x in {'.git','__pycache__','.venv'} for x in rel.parts):continue
        if p.suffix in {'.zip','.pyc'} or (len(rel.parts)>1 and rel.parts[:2]==('validation','work')):continue
        files.append(p)
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for p in sorted(files):
            name='peptide-side-product-diagnosis/'+p.relative_to(ROOT).as_posix()
            info=zipfile.ZipInfo(name,date_time=(2026,11,17,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o644<<16;archive.writestr(info,p.read_bytes())
    digest=hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix+'.sha256').write_text(digest+'  '+output.name+'\n',encoding='utf-8')
    print(json.dumps({'file':str(output),'files':len(files),'bytes':output.stat().st_size,'sha256':digest}))
if __name__=='__main__':main()
