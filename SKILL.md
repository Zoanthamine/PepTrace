---
name: peptide-side-product-diagnosis
description: 肽迹（PepTrace）：用多肽序列、端基或保护态、反应条件和质量／LC–MS m/z 查找副反应家族、并列位点和组合候选，生成带来源的诊断报告。
---

# 肽迹 · 多肽副反应与副产物诊断

本技能包含完整Python引擎、规范规则、质量组件和索引。正式入口是技能目录下 `scripts/diagnose.py`；以该目录为基准解析路径，不依赖某台电脑的工程目录。使用已安装Pyteomics、JSON Schema、RDKit的Python；MS/MS辅助另外需要PyOpenMS。依赖与安装方式见仓库README。

从用户原话提取序列、端基、修饰／连接、处理步骤、观测口径和实际误差。简单自然语言可以直接交给入口；多链、双峰或保护态由agent按[病例格式](references/case-format.md)整理内部JSON，保留原文。未给端基或误差的质量差可先粗查来源，并明确缺项。不要猜电荷、加合或第二个峰。

在技能目录执行：

```text
python scripts/diagnose.py --text "序列 MELK，N端游离，C端羧酸，比理论多15.994915 Da，容差0.001 Da，用DMSO处理过" --output result.json --report report.md
python scripts/diagnose.py case.json --output result.json --report report.md
python scripts/diagnose.py --lookup CH13-R014
python scripts/diagnose.py --find 甲酯
python scripts/diagnose.py --stats
```

生成文件使用用户指定目录或宿主的输出目录。默认最多2次规则应用，`--max-events`可明确配置1–6次；`--max-nodes`限制状态预算。检查结果里的搜索范围，截断时说明未搜完。

解释结果时区分：

- `candidates`：实际实例化的质量匹配结构假设；同净式不同位点、连接、固定电荷和支持路径保留，质量吻合不能单独确证结构。
- `condition_conflicts`：被明确否定前提的路径，解释具体排除原因；未知条件继续保留。
- `source_leads`、`dsh_source_review`：全书／附录或独立算法的来源线索；参照、结构和保护态未核对时，不能升级成已执行候选。
- `fusion_audit`：规范数据版本与DSH独立分子式核算。本体状态核心负责结构和组合，DSH旧加质量组合和旧条件表不接管推理。

原料和目标、平均／单同位素／原书整数、双峰误差及固有电荷分别处理。书中观察、编辑模型和当前推断按证据等级陈述；输出全部候选JSON，正文可以优先解释条件支持的少事件假设，并给出MS/MS、标品或工艺对照等区分证据。

本技能采用反应家族、净式、必要位点和条件的实用粒度。全部450条来源规则与82条附录有去向；240个模板、308个编译分支用于状态搜索，复杂或未定义结构继续作为来源线索。真实样品准确率、全部原子价态／几何及无限反应深度尚未验收。范围与来源见[能力和证据](references/scope-and-provenance.md)。

有MS2峰表或质心mzML时按[碎片辅助](references/msms.md)使用 `engine/tools/compare_msms.py`；缺PyOpenMS时明确该模块不可执行，不能假称算过。

修改规范 `engine/knowledge/*.json` 后执行 `python scripts/build.py` 重建运行包、索引和版本清单，再执行 `python scripts/verify.py --full`。入口核对文件哈希，避免混用不同数据版本。

常规诊断直接读取冻结规则，无需反复翻图；需要翻图时使用局部图并清理可再生缓存，避免单次请求累计接近50张图片，保留原始证据与文字记录。
