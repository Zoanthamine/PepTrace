"""Smoke tests for diagnose.py — synthetic, clearly fictional sequences/observed masses.

Each case asserts that the expected candidate rule appears in the ranked list.
Exit code 0 means all pass; failures are printed with the actual candidates.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from diagnose import diagnose, RULES, rule_delta, mz_of
from mass_ledger import peptide_mass, residue_mass, formula_mass

CASES = []

def case(seq, note, expected, **kw):
    base = peptide_mass(seq)
    if kw.get("nterm_ac"):
        base += formula_mass({"C": 2, "H": 2, "O": 1})
    delta = kw.pop("delta")
    neutral = base + delta
    CASES.append({"seq": seq, "note": note, "expected": expected,
                  "kw": {k: v for k, v in kw.items() if k != "delta"},
                  "neutral": neutral})

case("GDSGKS", "Asn脱酰胺 +0.984016", "DEAMID",
     delta=formula_mass({"O": 1, "N": -1, "H": -1}))
case("GSAC", "Cys→Dha −H2S", "CYS_DHA", delta=formula_mass({"H": -2, "S": -1}))
case("GSM", "Met单次氧化 +O", "OX_1O", delta=15.99491461957)
case("GSAG", "删除Ala1次", "deletion", delta=-residue_mass("A"))
case("RKR", "Arg→Orn", "ARG_ORN", delta=formula_mass({"C": -1, "H": -2, "N": -2}))
case("GSC", "脱酰胺 vs 无匹配（应为空）", None, delta=4.5)

def run() -> int:
    fails = 0
    for c in CASES:
        out = diagnose(c["seq"], neutral_obs=c["neutral"], tol_ppm=20, **c["kw"])
        ids = [x.get("rule_id") for x in out["candidates"]]
        if c["expected"] == "deletion":
            ok = any(x["source"].startswith("del") for x in out["candidates"])
        elif c["expected"]:
            ok = c["expected"] in ids
        else:
            ok = len(out["candidates"]) == 0
        print(f"[{'PASS' if ok else 'FAIL'}] {c['note']:36s} seq={c['seq']:8s} "
              f"candidates={[i for i in ids if i][:4]}{'(deletion)' if any(x['source']!='rule' for x in out['candidates']) else ''}")
        if not ok:
            fails += 1
            print("      got:", [(x['rule_id'], x['name_cn'] if 'name_cn' in x else x['name']) for x in out['candidates']])
    print(f"\n{fails} failure(s) of {len(CASES)}")
    return 1 if fails else 0

if __name__ == "__main__":
    raise SystemExit(run())
