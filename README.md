# 肽迹 · PepTrace

**多肽副反应与副产物诊断。** 输入多肽序列、端基／保护态、反应条件以及质量或LC–MS m/z，得到并列候选、位点、条件冲突、组合路径和来源报告。

PepTrace combines a peptide state-search engine with an independent DSH mass/formula checker. It runs locally in Python and ships with its rule library, mass components, index, examples, and tests. No API key is required for the deterministic CLI.

## 完整包包含什么

- Codex／DSH可发现的技能入口：[SKILL.md](SKILL.md)，名称为“肽迹”，技术ID保留 `peptide-side-product-diagnosis`。
- 可执行状态搜索、有限组合、结构化／受限自然语言输入、中文与JSON报告。
- 450条规范来源规则、240模板、308编译分支、82条附录去向，残基与保护基质量组件。
- 可从规范JSON重建的运行包、兼容导出和SQLite索引。
- 五个公开练习、参考答案修订、独立枚举与三个随机质量家族保留题。
- 回归测试、逐分支检查、哈希清单、可迁移安装和ZIP验证脚本。

原书ZIP保留在[v1.0.0 Release](https://github.com/Zoanthamine/PepTrace/releases/tag/v1.0.0)附件中，文件身份见 `references/source-archive.json`。它不作为运行时依赖。原书及其他第三方材料的权利归各自权利人，仓库MIT许可仅适用于本项目自编代码和说明。

下载后可运行 `python scripts/verify_source_archive.py <原书ZIP>`，独立核对压缩包及其中21份PDF的SHA-256。

## 安装

需要Python 3.11或更新版本；发布验证使用Python 3.12。先克隆到技能目录，也可以下载Release的完整ZIP后解压为同名目录：

```text
git clone https://github.com/Zoanthamine/PepTrace.git peptide-side-product-diagnosis
cd peptide-side-product-diagnosis
python -m pip install -r requirements.txt
python scripts/verify.py
```

要执行MS/MS和全部回归测试，再安装：

```text
python -m pip install -r requirements-msms.txt
python scripts/verify.py --full
```

Codex：将整个仓库目录放到 `<CODEX_HOME>/skills/peptide-side-product-diagnosis`。DSH：放到其技能目录的同名文件夹。两者读取同一份 `SKILL.md`；不要只复制说明文件而漏掉 `engine/`。已有科研Python环境可以直接使用，不必另建环境。依赖缺失时仅相应运行分支不可执行。

## 使用

用户可以提供普通实验记录；agent负责整理复杂病例，简单文本可直接运行：

```text
python scripts/diagnose.py --text "序列 MELK，N端游离，C端羧酸，比理论多15.994915 Da，容差0.001 Da，用DMSO处理过" --output result.json --report report.md
python scripts/diagnose.py examples/five_cases/case-1.json --max-events 1 --output result.json --report report.md
python scripts/diagnose.py --lookup CH13-R014
python scripts/diagnose.py --find 甲酯
python scripts/diagnose.py --stats
```

全部路径保存在JSON中，中文报告解释观测参照、误差、位点、条件、出处与边界。默认最多2次规则应用；可明确设置1–6次，预算耗尽会报告截断。低分辨整数线索与精确质量分开，离子电荷及加合不自动猜测。

输入例子与字段见 [病例格式](references/case-format.md)，MS/MS见 [碎片辅助](references/msms.md)。

## 如何理解候选

`candidates` 是实例化的质量匹配假设；`condition_conflicts` 是明确条件冲突的路径。`source_leads` 和 `dsh_source_review` 是来源反查，不能冒充已执行结构。同式位点、异构体、固定锍电荷及不同反应参照均保留。质量吻合需要其他谱证或工艺对照才能形成结构结论。

五题涉及额外Phe的N端／Lys侧链、Asp/Glu/C端甲酯、二硫→三硫、DCC阴性信息、MTBE叔丁基与Met锍。原参考答案的还原口径、+2O归属、分辨率和CH₂算术修订见 [案例说明](examples/five_cases/案例说明与答案修订.md)。公开答案例题不称盲测；真实样品准确率尚未验收。

## 重建与核验

```text
python scripts/build.py
python scripts/verify.py --full
python scripts/verify.py --full --benchmark
python scripts/package.py --output PepTrace-v1.1.0.zip
```

构建只需包内的规范JSON，不需要原书PDF。希望重新核对原书文件哈希时，先将原书ZIP解压并按来源目录组织为 `doc_XX/source.pdf`，然后使用 `python scripts/build.py --source-root <目录>`；缺失或哈希不一致会报错。重建数据不等于重新目视核验原图。

`--full` 执行原193项回归、打包迁移测试、308分支检查、独立枚举、三题随机保留检索及五个公开病例。`--benchmark` 执行已登记的性能测试，10秒门槛来自发布验证机器的短序列工作负载，不能当作其他机器或长序列的保证。每个门槛与验证结果保存在 `validation/`。

## 目录

```text
SKILL.md / agents/         技能入口与展示信息
scripts/                  诊断、重建、核验与打包入口
references/               病例、范围、来源、MS/MS说明
engine/knowledge/         唯一规范数据
engine/runtime/           冻结运行包
engine/index/             可重建SQLite索引
engine/tools/             执行核心与原回归
engine/fusion_vendor/     DSH独立核算模块与来源记录
examples/five_cases/      五个公开测试输入
tests/ / validation/      迁移测试与验收证据
```

基于 Yi Yang, *Side Reactions in Peptide Synthesis* (2016, Academic Press) 编制事实规则与编辑模型；规则中保留印刷页、图号及证据层级。原书报告、编辑推算和当前诊断假设分开。详细来源见 [范围与证据](references/scope-and-provenance.md) 和 [文件目录](references/source-catalog.json)。

代码与自编说明采用 [MIT](LICENSE)。第三方Python依赖遵循其自身许可；原书及原书ZIP不适用MIT许可。

## 基础实验辅助（1.1.0）

新增 `scripts/assist.py`，共用现有规范质量与证据口径：

|模式|基础能力|
|---|---|
|spps|组成、重复片段、Asp-X/Asp-Pro及实际过程缺项；不预测成功率|
|hplc|梯度/柱尺寸换算和TFA、酸/甲醇、pH兼容性提示|
|ms1|已索引质心mzML或JSON的TIC/BPC/XIC、明确窗口及积分|
|cyclic|给定SMILES的映射环键、分子式、形式电荷与酸前体失水关系|
|protocols|本地规程快照版本/材料/步骤；官方GET显式执行与token缺项|

格式与示例见 [基础实验辅助](references/experimental-basics.md)。上游方法、许可、实际读取范围与可加入的书籍/网站见 [来源索引](references/experimental-sources.json) 和 [知识扩充计划](references/knowledge-expansion.md)。所有新增基础脚本独立编写；未纳入未标许可的SPPS扫描器、非商业环肽运行器或GPU模型。现有诊断核心的450来源规则和308编译分支保持原口径。
