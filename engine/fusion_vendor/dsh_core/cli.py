"""命令行入口。

用法
----
单条病例（最常用）::

    python -m dsh_core.cli case.json                # 输出中文 Markdown 报告
    python -m dsh_core.cli case.json --json out.json --top-k 8
    python -m dsh_core.cli - < case.json            # 从 stdin 读

自然语言（受限解析，失败会明确报错而不会猜）::

    python -m dsh_core.cli --text "序列 MELK，实测比理论多出 15.99 Da，用 DMSO 处理过"

按序列/条件查（不打质量）::

    python -m dsh_core.cli --lookup --residue C --chapter 2
    python -m dsh_core.cli --lookup --reagent Piperidine
    python -m dsh_core.cli --rule CH02-R002         # 看单条规则全文

其他::

    python -m dsh_core.cli --stats                  # 知识库统计与自检
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .candidates import diagnose
from .pipeline import run_case
from .components import ComponentError, load_library
from .report import render_json, render_markdown
from .rules import load_rules
from .schema import InputError, parse_case


def _read_json(path: str) -> dict:
    if path == "-":
        return json.load(sys.stdin)
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"找不到输入文件：{p}")
    return json.loads(p.read_text(encoding="utf-8"))


def cmd_case(args: argparse.Namespace, lib, ruleset) -> int:
    if getattr(args, "text", None):
        from .nl_input import parse_text

        raw = parse_text(args.text, case_id=args.case_id)
    else:
        raw = _read_json(args.case)
    case = parse_case(raw, lib)
    res = run_case(case, ruleset, lib, combine=not args.no_combine, top_k=args.top_k,
                   max_events=args.max_events, max_nodes=args.max_nodes,
                   max_results=args.max_results)
    print(res.markdown)
    if args.json:
        Path(args.json).write_text(json.dumps(res.payload, ensure_ascii=False, indent=2),
                                   encoding="utf-8")
        print(f"\n[已写出 JSON：{args.json}]", file=sys.stderr)
    return 0


def cmd_lookup(args: argparse.Namespace, lib, ruleset) -> int:
    hits = list(ruleset.rules.values())
    if args.chapter:
        hits = [r for r in hits if r.chapter == args.chapter]
    if args.residue:
        aa = args.residue.strip()
        hits = [r for r in hits if aa in r.named_residues or aa.upper() in r.susceptible.join("").upper()]
    if args.reagent:
        key = args.reagent.strip().lower()
        hits = [r for r in hits if any(key == x.lower() for x in r.required_any + r.contextual)]
    if args.status:
        hits = [r for r in hits if r.mass_status == args.status]
    if args.keyword:
        kw = args.keyword.lower()
        hits = [r for r in hits if kw in (r.name_cn + r.name_en + r.mechanism + r.conditions).lower()]
    hits.sort(key=lambda r: r.rule_id)
    print(f"# 命中 {len(hits)} 条规则\n")
    print("| 规则 | 名称 | 章 | 页码 | 质量差 | 条件要求 | 相关残基 |")
    print("|---|---|---|---|---|---|---|")
    for r in hits[: args.limit]:
        print(
            f"| {r.rule_id} | {r.name_cn} | {r.chapter} | {'/'.join(map(str, r.print_pages))} | "
            f"{r.mass_label()} | {'/'.join(r.required_any) or '—'} | {'/'.join(r.named_residues) or '—'} |"
        )
    if len(hits) > args.limit:
        print(f"\n（还有 {len(hits) - args.limit} 条，用 --limit 放开）")
    return 0


def cmd_rule(args: argparse.Namespace, lib, ruleset) -> int:
    r = ruleset.rules.get(args.rule)
    if not r:
        raise SystemExit(f"没有这条规则：{args.rule}")
    print(f"# {r.rule_id} {r.name_cn}（{r.name_en}）\n")
    print(f"- 出处：第 {r.chapter} 章 {r.section}；印刷页 {list(r.print_pages)}；PDF 页 {list(r.pdf_pages)}"
          + (f"；图 {list(r.figures)}" if r.figures else ""))
    print(f"- 来源：{'书中' if r.origin == 'book' else '编辑补充（不得当原书数据引用）'}")
    print(f"- 阶段：{'；'.join(r.stage) or '—'}")
    print(f"- 条件：{r.conditions or '—'}")
    print(f"- 敏感序列：{'；'.join(r.susceptible) or '—'}")
    print(f"- 产物：{r.products or '—'}")
    print(f"- 机理：{r.mechanism or '—'}")
    print(f"- 质量状态：{r.mass_status}")
    print(f"- 质量差：{r.mass_label()}" + (f"（净式 {r.delta.name()}）" if r.delta else ""))
    print(f"- 质量表达：{r.mass_expression or '—'}")
    print(f"- 诊断证据：{r.diagnostic_evidence or '—'}")
    print(f"- 竞争解释：{r.competing or '—'}")
    print(f"- 缓解：{r.mitigation or '—'}")
    print(f"- 条件要求（机器判定用）：需要 {'/'.join(r.required_any) or '（无权威要求）'}；"
          f"旁证 {'/'.join(r.contextual) or '—'}；相关残基 {'/'.join(r.named_residues) or '—'}")
    return 0


def cmd_stats(args: argparse.Namespace, lib, ruleset) -> int:
    st = ruleset.stats()
    print("# 知识库统计（DSH 版 vnext_dsh）\n")
    print(f"- 规则总数：{st['total']}（书中 {st['book']}，编辑补充 {st['supplement']}）")
    print(f"- 可计算精确质量差：{st['exact_mass']}；仅整数 amu：{st['nominal_only']}；无单一质量差：{st['no_mass']}")
    print(f"- 质量状态分布：{json.dumps(st['mass_status'], ensure_ascii=False)}")
    print(f"- 组件库：{json.dumps(lib.stats(), ensure_ascii=False)}")
    print(f"- 含权威条件要求的规则：{sum(1 for r in ruleset.rules.values() if r.required_any)}")
    print(f"- 数据自检问题：{st['data_problems']}")
    for p in ruleset.problems[:20]:
        print(f"  - {p['rule_id']}: {p['detail']}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="多肽副产物诊断（DSH 版）")
    ap.add_argument("case", nargs="?", help="病例 JSON 路径，'-' 表示 stdin")
    ap.add_argument("--json", help="把结构化结果写到该路径")
    ap.add_argument("--top-k", type=int, default=None, help="报告里最多展示多少条候选（不影响 JSON 完整性判断）")
    ap.add_argument("--no-combine", action="store_true",
                    help="关闭 T06 组合搜索（默认开启：在声明范围内枚举 1–2 个事件）")
    ap.add_argument("--max-events", type=int, default=2, help="组合事件上限（当前引擎只支持 2）")
    ap.add_argument("--max-nodes", type=int, default=200_000, help="搜索节点预算（触顶会显式报截断）")
    ap.add_argument("--max-results", type=int, default=5_000, help="结果数上限（触顶会显式报截断）")
    ap.add_argument("--text", help="中文自然语言病例描述（受限解析；解析不了会报错，不猜）")
    ap.add_argument("--case-id", help="给自然语言病例指定 case_id")
    ap.add_argument("--lookup", action="store_true", help="按序列/条件/试剂查规则（不打质量）")
    ap.add_argument("--rule", help="打印单条规则全文")
    ap.add_argument("--stats", action="store_true", help="知识库统计与自检")
    ap.add_argument("--chapter", type=int)
    ap.add_argument("--residue")
    ap.add_argument("--reagent")
    ap.add_argument("--status")
    ap.add_argument("--keyword")
    ap.add_argument("--limit", type=int, default=40)
    args = ap.parse_args(argv)

    try:
        lib = load_library()
    except ComponentError as exc:
        print(f"[组件库错误] {exc}", file=sys.stderr)
        return 2
    ruleset = load_rules(lib)

    if args.stats:
        return cmd_stats(args, lib, ruleset)
    if args.rule:
        return cmd_rule(args, lib, ruleset)
    if args.lookup:
        return cmd_lookup(args, lib, ruleset)
    if not args.case and not args.text:
        ap.print_help()
        return 1
    try:
        return cmd_case(args, lib, ruleset)
    except InputError as exc:
        print(f"输入错误 [{exc.code}] {exc.field}：{exc.message}", file=sys.stderr)
        if exc.hint:
            print(f"  → {exc.hint}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    sys.exit(main())