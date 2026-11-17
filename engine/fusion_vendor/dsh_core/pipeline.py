"""T06 流水线：单事件匹配（T04/T03）→ 组合搜索（T06）→ 报告（T07 雏形）。

这样 CLI、示例运行器、验收脚本共用同一条路径，避免三处各自拼装导致口径不一致。
"""
from __future__ import annotations

from dataclasses import dataclass

from .candidates import Diagnosis, diagnose
from .components import ComponentLibrary
from .report import render_json, render_markdown
from .rules import RuleSet
from .schema import CaseInput
from .search import SearchResult, search_combinations


@dataclass
class RunResult:
    diagnosis: Diagnosis
    search: SearchResult | None
    markdown: str
    payload: dict


def run_case(case: CaseInput, ruleset: RuleSet, lib: ComponentLibrary,
             combine: bool = True, top_k: int | None = None,
             max_events: int = 2, max_nodes: int = 200_000,
             max_results: int = 5_000) -> RunResult:
    """跑一条病例：单事件诊断 + 可选的组合搜索，返回报告与结构化结果。"""
    diag = diagnose(case, ruleset, lib)
    search = None
    if combine:
        search = search_combinations(case, ruleset, lib, max_events=max_events,
                                     max_nodes=max_nodes, max_results=max_results)
        # 组合层排除的条目并入诊断排除清单的可读说明（不改变单事件口径的统计）
        if search.excluded:
            diag.notes.append(
                f"组合层另有 {len(search.excluded)} 条单事件因你明确否定的条件被排除；"
                "单事件口径的统计未计入这些条目。")
    md = render_markdown(diag, ruleset, top_k=top_k, search=search)
    payload = render_json(diag, ruleset, top_k=top_k, search=search)
    return RunResult(diagnosis=diag, search=search, markdown=md, payload=payload)