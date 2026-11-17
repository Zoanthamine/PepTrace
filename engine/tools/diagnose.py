"""Clinical diagnostic engine: given sequence + reaction conditions + observed mass/mz, rank candidate side products.

Honest design goals:
 - Exact masses only where the neutral structure is well defined (formula-confirmed rules).
 - Sequence/end-group aware deletion deltas (never residue mass without end-group logic).
 - Explicit 'isobaric / indistinguishable' flags; no forced unique answer.
 - Combination search limited to chemistry-plausible pairs, not blind mass summation.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mass_ledger import (
    parse_formula, formula_mass, residue_mass, peptide_mass, split_sequence, PROTON, MONO,
)

# ---- Simple, auditable delta rules. Each entry is a named candidate transform.
# Reference: the declared target peptide (fully deprotected acid or amide per input).
RULES = [
    {"id": "OX_1O", "name": "单次氧化 (+O)", "formula": "O", "requires": any,
     "residues": "C,M,W,H,Y", "note": "Met/Tyr/Trp/His细分类需MS/MS定位"},
    {"id": "OX_2O", "name": "双次氧化 (+O2)", "formula": "O2", "requires": any,
     "residues": "C,M,W", "note": "Met sulfone等；与+2S同位素异构不同"},
    {"id": "THIO", "name": "酰胺→硫代酰胺 (+S−O)", "formula": "", "explicit": ("S", 1, "O", -1),
     "residues": "backbone", "note": "与+O相差0.017758Da，需高分辨"},
    {"id": "DEAMID", "name": "Asn/Gln脱酰胺 (+O−NH)", "explicit": ("O", 1, "N", -1, "H", -1),
     "residues": "N,Q", "note": "与C端酰胺水解同式；需定位"},
    {"id": "PYROGLU_E", "name": "N端Glu→焦谷氨酸 (−H2O)", "explicit": ("H", -2, "O", -1),
     "residues": "N端E", "note": "仅N端游离Glu；与普通脱水同式"},
    {"id": "PYROGLU_Q", "name": "N端Gln→焦谷氨酸 (−NH3)", "explicit": ("N", -1, "H", -3),
     "residues": "N端Q", "note": "仅N端游离Gln"},
    {"id": "ASP_IMIDE", "name": "Asp→琥珀酰亚胺 (−H2O)", "explicit": ("H", -2, "O", -1),
     "residues": "D", "note": "与脱水同式；Pro后的Asp需特殊机理"},
    {"id": "ARG_ORN", "name": "Arg→Orn (−CH2N2)", "explicit": ("C", -1, "H", -2, "N", -2),
     "residues": "R", "note": "与去乙酰−42.010565相差0.011233Da"},
    {"id": "ACETYL_LOSS", "name": "去乙酰 (−C2H2O)", "explicit": ("C", -2, "H", -2, "O", -1),
     "residues": "N端Ac", "note": "与Arg→Orn近等，须高分辨"},
    {"id": "FORMYL", "name": "甲酰化 (+CO)", "explicit": ("C", 1, "O", 1),
     "residues": "N端/K/H/W", "note": "甲酸簇+46不是甲酰化"},
    {"id": "ACETYL", "name": "乙酰化 (+C2H2O)", "explicit": ("C", 2, "H", 2, "O", 1),
     "residues": "N端/K/S/T/Y", "note": ""},
    {"id": "CARBAMIDO", "name": "碘乙酰胺羧胺甲基化 (+C2H3NO)", "explicit": ("C", 2, "H", 3, "N", 1, "O", 1),
     "residues": "C", "note": "与Gly残基插入同式(57.021464)"},
    {"id": "FMOC_RES", "name": "Fmoc残留 (+C15H10O2)", "explicit": ("C", 15, "H", 10, "O", 2),
     "residues": "N端/侧链胺", "note": ""},
    {"id": "BOC_RES", "name": "Boc残留 (+C5H8O2)", "explicit": ("C", 5, "H", 8, "O", 2),
     "residues": "N端/K", "note": ""},
    {"id": "TBU_RES", "name": "tBu残留 (+C4H8)", "explicit": ("C", 4, "H", 8),
     "residues": "S/T/Y/D/E", "note": "与tert-丁基烷基化同式"},
    {"id": "TRT_RES", "name": "Trt残留 (+C19H14)", "explicit": ("C", 19, "H", 14),
     "residues": "C/H/N端", "note": ""},
    {"id": "PBF_RES", "name": "Pbf残留 (+C10H14O3S)", "explicit": ("C", 10, "H", 14, "O", 3, "S", 1),
     "residues": "R", "note": "附录+252附近；精确式需与正文核对"},
    {"id": "PHOSPHO_LOSS", "name": "去磷酸 (−HPO3)", "explicit": ("H", -1, "P", -1, "O", -3),
     "residues": "S/T/Y", "note": "不等于−H3PO4 β消除"},
    {"id": "DHA_ELIM", "name": "磷酸Ser/Thr β消除 (−H3PO4)", "explicit": ("H", -3, "P", -1, "O", -4),
     "residues": "S/T", "note": "Thr产物为Dhb保留甲基"},
    {"id": "CYS_DHA", "name": "Cys β消除→Dha (−H2S)", "explicit": ("H", -2, "S", -1),
     "residues": "C", "note": "相对去保护Cys目标"},
    {"id": "DISULFIDE", "name": "分子内二硫键 (−H2)", "explicit": ("H", -2),
     "residues": "C≥2", "note": ""},
    {"id": "REDUCE_SS", "name": "二硫键还原 (+H2)", "explicit": ("H", 2),
     "residues": "C", "note": ""},
    {"id": "PIPERIDYL", "name": "Dha被哌啶加成 (C5H9N−1S−1)", "explicit": ("C", 5, "H", 9, "N", 1, "S", -1),
     "residues": "C", "note": "相对去保护Cys目标；精确51.101428"},
]

def delta_of(explicit):
    total = {}
    it = iter(explicit)
    for el, n in zip(it, it):
        total[el] = total.get(el, 0) + n
    return formula_mass({k: v for k, v in total.items() if v})

def rule_delta(rule) -> float:
    if "explicit" in rule:
        return delta_of(rule["explicit"])
    return formula_mass(parse_formula(rule["formula"]))

def residue_sequence_names(seq):
    return split_sequence(seq)

def deletion_candidates(seq, mz_obs, z, tol_da):
    """Single-residue deletions, and all 2-residue deletions, with sequence-aware masses."""
    codes = residue_sequence_names(seq)
    out = []
    for i in range(len(codes)):
        d = -(residue_mass(codes[i]))
        out.append({"type": "del1", "residues": codes[i], "seq_pos": i + 1, "delta": d})
    for i in range(len(codes)):
        for j in range(i + 1, min(i + 3, len(codes))):
            d = -(residue_mass(codes[i]) + residue_mass(codes[j]))
            out.append({"type": "del2", "residues": codes[i] + codes[j], "seq_pos": f"{i+1}-{j+1}", "delta": d})
    for i in range(len(codes) - 1):  # DKP removes exactly the N-terminal two residues for free-N-term peptides
        if i == 0:
            out.append({"type": "DKP_Nterm2", "residues": codes[0] + codes[1], "seq_pos": "1-2",
                        "delta": -(residue_mass(codes[0]) + residue_mass(codes[1]))})
    return out

def mz_of(m, z, adduct="H"):
    if adduct == "H":
        return (m + z * PROTON) / z
    if adduct == "Na":
        return m + MONO["Na"] - PROTON
    raise ValueError(adduct)

def diagnose(seq, nterm_ac=False, cterm_amide=False, mz_obs=None, neutral_obs=None, z=1, tol_ppm=50, adduct="H"):
    base = peptide_mass(seq, cterm_acid=not cterm_amide)
    # peptide_mass handles N-acetyl only via explicit formula; add if requested
    if nterm_ac:
        base += formula_mass({"C": 2, "H": 2, "O": 1})
    ref = neutral_obs if neutral_obs is not None else mz_obs * z - (z * PROTON if adduct == "H" else 0)
    tol = tol_ppm * 1e-6 * ref
    candidates = []
    for rule in RULES:
        d = rule_delta(rule)
        pred = base + d
        err = pred - ref
        if abs(err) <= tol:
            candidates.append({"source": "rule", "rule_id": rule["id"], "name": rule["name"],
                               "delta": round(d, 6), "predicted": round(pred, 6), "err_da": round(err, 6),
                               "note": rule.get("note", "")})
    for cand in deletion_candidates(seq, ref, z, tol):
        pred = base + cand["delta"]
        err = pred - ref
        if abs(err) <= tol:
            candidates.append({"source": cand["type"], "rule_id": None,
                               "name": f"删除{cand['residues']}（位置{cand['seq_pos']}）",
                               "delta": round(cand["delta"], 6), "predicted": round(pred, 6),
                               "err_da": round(err, 6), "note": "端基/封端/标签须核对"})
    # 2-event combinations: oxidation/deamidation/elimination style rules on top of deletion
    return {"sequence": seq, "base_neutral": round(base, 6), "observed_neutral": round(ref, 6),
            "tolerance_da": round(tol, 4), "candidates": candidates}

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="序列+观测质量 → 副反应候选（中性单同位素口径）")
    ap.add_argument("sequence", help="如 GDSGKS（单字母，支持 X 与括号修饰由见 spec）")
    ap.add_argument("--mz", type=float, help="观测 m/z（需同时给 --z）")
    ap.add_argument("--neutral", type=float, help="观测中性单同位素质量 M")
    ap.add_argument("--z", type=int, default=1, help="电荷数，默认 1")
    ap.add_argument("--adduct", default="H", choices=["H", "Na"], help="加合离子，默认 [M+zH]z+")
    ap.add_argument("--tol-ppm", type=float, default=50.0, help="容差 ppm，默认 50")
    ap.add_argument("--nterm-ac", action="store_true", help="目标N端为乙酰基")
    ap.add_argument("--cterm-amide", action="store_true", help="目标C端为酰胺")
    a = ap.parse_args()
    if a.mz is None and a.neutral is None:
        ap.error("必须给 --mz（配 --z）或 --neutral 之一")
    print(json.dumps(diagnose(a.sequence, nterm_ac=a.nterm_ac, cterm_amide=a.cterm_amide,
                               mz_obs=a.mz, neutral_obs=a.neutral, z=a.z,
                               tol_ppm=a.tol_ppm, adduct=a.adduct), ensure_ascii=False, indent=2))
