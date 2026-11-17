"""质量组件库加载与自校验。

组件库 `多肽质量组件.json` 是既有资产（由原始项目生成），本模块只**读取**、不改写。
加载时用本包独立的原子量表重算每一个组件的质量，与文件中记录值比对；
不一致就报错——既不静默采用文件的数，也不静默采用我算的数。
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .formula import Formula, FormulaError


def kb_root() -> Path:
    """知识库根目录：默认取本包上两级（vnext_dsh 的父目录），可用 DSH_KB_ROOT 覆盖。"""
    env = os.environ.get("DSH_KB_ROOT")
    if env:
        return Path(env).resolve()
    return Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Component:
    name: str
    formula: Formula
    mass: float
    note: str = ""

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.name}({self.formula.name()}, {self.mass:.6f})"


class ComponentError(RuntimeError):
    pass


class ComponentLibrary:
    """残基、自由氨基酸、共价净增量、中性小分子的登记表。"""

    def __init__(self, data: dict, tol: float = 1e-3, source_path: Path | None = None) -> None:
        self.tol = tol
        self._raw_data = data
        self.source_path = source_path
        self.residues: dict[str, Component] = {}
        self.free_amino_acids: dict[str, Component] = {}
        self.deltas: dict[str, Component] = {}
        self.neutral: dict[str, Component] = {}
        self.sequences_free_ends: dict[str, dict] = {}
        self._load(data)
        self._verify()

    # ---------- 加载 ----------
    def _load(self, data: dict) -> None:
        for key, name, table in (
            ("residues", "residue", self.residues),
            ("residues", "free_amino_acid", self.free_amino_acids),
        ):
            for aa, rec in data.get(key, {}).items():
                part = rec.get(name)
                if not part:
                    continue
                f = Formula.of(part["formula"])
                table[aa] = Component(aa, f, f.mass(), rec.get("note", ""))
        for gname, rec in data.get("covalent_deltas", {}).items():
            f = Formula.of(rec["formula"])
            self.deltas[gname] = Component(gname, f, f.mass(), rec.get("reference", ""))
        for nname, rec in data.get("neutral_molecules", {}).items():
            f = Formula.of(rec["formula"])
            self.neutral[nname] = Component(nname, f, f.mass(), rec.get("note", ""))
        self.sequences_free_ends = data.get("sequences_free_ends", {})

    # ---------- 自校验 ----------
    def _verify(self) -> None:
        """用独立原子量表重算全部组件质量；与文件记录值偏差超容差即报错。"""
        problems: list[str] = []
        raw = self._raw
        for aa, comp in self.residues.items():
            rec = raw["residues"][aa]["residue"]
            if "monoisotopic_da" in rec:
                d = abs(rec["monoisotopic_da"] - comp.mass)
                if d > self.tol:
                    problems.append(f"残基 {aa}: 文件 {rec['monoisotopic_da']:.6f} 与重算 {comp.mass:.6f} 差 {d:.2e}")
        for aa, comp in self.free_amino_acids.items():
            rec = raw["residues"][aa]["free_amino_acid"]
            if "monoisotopic_da" in rec:
                d = abs(rec["monoisotopic_da"] - comp.mass)
                if d > self.tol:
                    problems.append(f"自由氨基酸 {aa}: 文件 {rec['monoisotopic_da']:.6f} 与重算 {comp.mass:.6f} 差 {d:.2e}")
        for table, key in ((self.deltas, "covalent_deltas"), (self.neutral, "neutral_molecules")):
            for name, comp in table.items():
                rec = raw[key][name]
                if "monoisotopic_da" in rec:
                    d = abs(rec["monoisotopic_da"] - comp.mass)
                    if d > self.tol:
                        problems.append(f"{key} {name}: 文件 {rec['monoisotopic_da']:.6f} 与重算 {comp.mass:.6f} 差 {d:.2e}")
        if problems:
            raise ComponentError("质量组件库自校验失败：\n  " + "\n  ".join(problems))

    # ---------- 查询 ----------
    def residue(self, aa: str) -> Component:
        if aa not in self.residues:
            raise ComponentError(f"不支持残基 {aa!r}；已登记：{''.join(sorted(self.residues))}")
        return self.residues[aa]

    def delta(self, group: str) -> Component:
        if group not in self.deltas:
            raise ComponentError(f"未登记共价净增量 {group!r}")
        return self.deltas[group]

    def has_delta(self, group: str) -> bool:
        return group in self.deltas

    def neutral_molecule(self, name: str) -> Component:
        if name not in self.neutral:
            raise ComponentError(f"未登记中性分子 {name!r}")
        return self.neutral[name]

    @property
    def _raw(self) -> dict:
        return self._raw_data

    def stats(self) -> dict:
        return {
            "residues": len(self.residues),
            "free_amino_acids": len(self.free_amino_acids),
            "covalent_deltas": len(self.deltas),
            "neutral_molecules": len(self.neutral),
            "source": str(self.source_path),
        }


@lru_cache(maxsize=4)
def load_library(path: str | None = None) -> ComponentLibrary:
    p = Path(path) if path else kb_root() / "多肽质量组件.json"
    if not p.exists():
        raise ComponentError(f"找不到质量组件库：{p}")
    data = json.loads(p.read_text(encoding="utf-8"))
    return ComponentLibrary(data, source_path=p)