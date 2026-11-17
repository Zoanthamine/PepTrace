"""PepTrace portable entry: verify frozen files, then run the accepted fused engine."""
import json,runpy,sys
from pathlib import Path
from check_release import check
ROOT=Path(__file__).resolve().parents[1]
if sys.argv[1:]==['--version']:
    print('PepTrace 1.0.0');raise SystemExit(0)
try:verified=check()
except (ValueError,OSError) as error:raise SystemExit(str(error))
if sys.argv[1:]==['--verify-release']:
    print(json.dumps(verified,ensure_ascii=False));raise SystemExit(0)
sys.path.insert(0,str(ROOT/'engine/tools'))
runpy.run_path(str(ROOT/'engine/tools/fused_diagnosis.py'),run_name='__main__')
