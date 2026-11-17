"""PepTrace shared basic experimental entry, using one explicit JSON input and one report."""
import argparse,json,sys
from pathlib import Path
from check_release import check,ROOT
from _experimental import spps_review,hplc_review,cyclic_review
from _ms1 import extract_scans,mzml_scans
from _protocols import normalize_snapshot,fetch_snapshot

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['spps','hplc','cyclic','ms1','protocols'])
    parser.add_argument('input',nargs='?',type=Path);parser.add_argument('--targets',type=Path)
    parser.add_argument('--fetch');parser.add_argument('--execute',action='store_true');parser.add_argument('--output',type=Path);parser.add_argument('--report',type=Path)
    a=parser.parse_args()
    try:
        check()
        if a.mode=='protocols' and a.fetch:
            if a.input:raise ValueError('Provide a local snapshot or an exact fetch URL, not both')
            result=fetch_snapshot(a.fetch,a.execute)
        else:
            if not a.input:raise ValueError('A local input file is required')
            if a.fetch or a.execute:raise ValueError('Network flags apply only to protocols mode')
            if a.mode=='ms1':
                if not a.targets:raise ValueError('MS1 extraction requires --targets with explicit m/z windows')
                spec=json.loads(a.targets.read_text(encoding='utf-8-sig'))
                scans=json.loads(a.input.read_text(encoding='utf-8-sig'))['scans'] if a.input.suffix.lower()=='.json' else mzml_scans(a.input)
                result=extract_scans(scans,spec)
            else:
                data=json.loads(a.input.read_text(encoding='utf-8-sig'))
                result={'spps':spps_review,'hplc':hplc_review,'cyclic':cyclic_review,'protocols':normalize_snapshot}[a.mode](data)
        result['peptrace_version']='1.1.0';result['source_input_sha256']=None
        registry=json.loads((ROOT/'references/experimental-sources.json').read_text(encoding='utf-8'))['sources']
        mapping={'spps':{'SPPS-REVIEW','YANG-2016','BEHRENDT-2016'},'hplc':{'HPLC-METHOD','HPLC-BASICS','YANG-2016'},
            'cyclic':{'CYCLIC-EVAL'},'ms1':{'PYOPENMS-SKILL','OPENMS-DOCS'},'protocols':{'PROTOCOLS-SKILL','PROTOCOLS-API'}}
        result['source_notes']=[r for r in registry if r['id'] in mapping[a.mode]]
        if a.input:
            import hashlib
            digest=hashlib.sha256()
            with a.input.open('rb') as stream:
                for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
            result['source_input_sha256']=digest.hexdigest()
        text=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
        if a.output:a.output.write_text(text,encoding='utf-8',newline='\n')
        else:print(text,end='')
        if a.report:
            lines=['# 肽迹基础实验辅助','',f"模块：{result['module']}。",'',str(result.get('scope','只读请求计划。')),'','```json',text.rstrip(),'```']
            a.report.write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    except (ValueError,KeyError,TypeError,ImportError,OSError) as error:parser.error(str(error))
if __name__=='__main__':main()
