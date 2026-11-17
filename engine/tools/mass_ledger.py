"""Neutral monoisotopic mass ledger for peptide side-reaction diagnostics.

All deltas are M(candidate neutral product) - M(declared neutral reference),
in Da, monoisotopic. Never mix a charge/ion-adduct mass into a covalent delta.
"""
from __future__ import annotations
from knowledge_source import mass_constants, components, residue_table
import json, math, re
from pathlib import Path

MONO = mass_constants()['monoisotopic_atoms']
PROTON = mass_constants()['proton_da']


def has_neutral_lookup_delta(rule):
    """Eligibility for mass lookup, before checking the rule's reference and substrate."""
    if rule.get("automatic_mass_event") is False:
        return False
    formula = rule.get("delta_formula")
    return ((rule["mass_status"] == "formula_confirmed" and isinstance(formula, dict) and bool(formula)
             and all(k in MONO and type(n) is int for k, n in formula.items())
             or rule["mass_status"] == "isomeric" and rule.get("delta_mass_mono_da") == 0)
            and rule["rule_id"] != "SUP-039"
            and isinstance(rule.get("delta_mass_mono_da"), (int, float)))

# Residue formulas (free amino acid - H2O) for 20 standard + specials
RESIDUE = residue_table()
NAME2CODE = components()['residue_aliases']
FORMULA_RE = re.compile(r"([A-Z][a-z]?)(-?\d*)")

def parse_formula(formula: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    pos = 0
    while pos < len(formula):
        m = FORMULA_RE.match(formula, pos)
        if not m:
            raise ValueError(f"bad formula at {pos}: {formula!r}")
        el, n = m.group(1), m.group(2)
        count = int(n) if n not in ("", "+", "-") else (1 if n == "" else int(n + "1"))
        counts[el] = counts.get(el, 0) + count
        pos = m.end()
    return {k: v for k, v in counts.items() if v != 0}

def formula_mass(counts: dict[str, int]) -> float:
    return sum(MONO[e] * n for e, n in counts.items())

def combine(*formulas: str, signs=(1,)) -> dict[str, int]:
    total: dict[str, int] = {}
    for f, s in zip(formulas, signs):
        for e, n in parse_formula(f).items():
            total[e] = total.get(e, 0) + n * s
    return {e: n for e, n in total.items() if n != 0}

def split_sequence(seq: str) -> list[str]:
    tokens = []
    i = 0
    while i < len(seq):
        three = seq[i:i + 3]
        if three.capitalize() in NAME2CODE and len(three) == 3 and not seq[i:i+1].isupper():
            pass
        if three in NAME2CODE:
            tokens.append(NAME2CODE[three]); i += 3; continue
        ch = seq[i].upper()
        if ch in RESIDUE:
            tokens.append(ch); i += 1; continue
        raise ValueError(f"unparsed residue at {i}: {seq[i:i+4]!r}")
    return tokens

def residue_mass(code: str) -> float:
    return formula_mass(parse_formula(RESIDUE[code]))

def peptide_mass(seq: str, nterm_free: bool = True, cterm_acid: bool = True) -> float:
    """Linear peptide: sum residues + H2O for free acid/free amine.
    N-acetyl: add C2H2O. C-terminal amide: -O + NH? use explicit modifiers."""
    codes = split_sequence(seq)
    m = sum(residue_mass(c) for c in codes) + formula_mass(parse_formula("H2O"))
    if not cterm_acid:
        m = m + formula_mass(parse_formula("NH")) - formula_mass(parse_formula("O"))  # COOH -> CONH2 : +NH - O
    return m

def deletion_delta(residues: list[str]) -> float:
    return -sum(residue_mass(r) for r in residues)

def mz(neutral_mass: float, z: int, adduct: str = "H") -> float:
    """[M+zH]z+ if adduct H; Na for [M+Na]+ etc."""
    if adduct == "H":
        return (neutral_mass + z * PROTON) / z
    if adduct == "Na":
        return neutral_mass + MONO["Na"] - PROTON
    raise ValueError(adduct)

def write_table(out: Path):
    """Emit a machine-checkable residue and delta table."""
    table = {
        "residues": {c: {"formula": RESIDUE[c], "mono": round(residue_mass(c), 6)} for c in RESIDUE},
        "water": round(formula_mass(parse_formula("H2O")), 6),
        "proton": PROTON,
    }
    out.parent.mkdir(exist_ok=True, parents=True)
    out.write_text(json.dumps(table, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    print(json.dumps({
        "gly": residue_mass("G"),
        "arg2orn": round(residue_mass("R") - residue_mass("O"), 6),
    }, ensure_ascii=False))
