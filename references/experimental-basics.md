# 肽迹的基础实验辅助

本模块共用规范质量组件、病例参照、UTF-8输入输出和发布哈希。实验数据、风险线索、物质结构和预测模型保持不同口径。用户提供普通实验记录时，由agent按下面格式整理内部文件，不要求用户编写JSON。

## SPPS：合成前及过程核对

```text
python scripts/assist.py spps examples/basics/spps.json --output spps.json --report spps.md
```

输入 `peptide`（序列、端基、实际修饰）和 `process`（chemistry、resin、loading_mmol_g、coupling、deprotection、cleavage、conditions）。输出序列组成、重复片段、Asp-X/Asp-Pro、Cys/His、Met、N端Q/E等注意项，列未给工艺参数。不给成功概率或纯度评分。受保护肽-树脂的行为不能直接等同于最终自由肽；不要用此粗筛自动换保护基或工艺。

`spps-synthesizability-review`作为评审方法来源。该仓库未明确标注许可，本模块代码独立编写；具体化学注意项同时引用Yang规则和已读取的开放综述。没有复制其扫描器或未验证的严重度阈值。

## HPLC：方法检查与转移算术

```text
python scripts/assist.py hplc examples/basics/hplc.json --output hplc.json --report hplc.md
```

输入实际柱长度／内径、流速、梯度时间与%B、添加剂、pH和可选新柱。输出几何柱体积、梯度斜率、保持线速度和梯度柱体积的转移建议，提示TFA/ESI、非挥发添加剂、酸/甲醇及柱pH范围。几何体积不是实测死体积；保留时间、分辨率、压力、进样负载和dwell volume仍需实际方法验证。不要把厂家单例的梯度/浓度当全部多肽的默认SOP。

## LC–MS：基础MS1迹线

```text
python scripts/assist.py ms1 sample.mzML --targets examples/basics/ms1-targets.json --output traces.json
python scripts/assist.py ms1 examples/basics/ms1-scans.json --targets examples/basics/ms1-targets.json --output traces.json
```

支持已标注质心的、已索引mzML和小型质心JSON。窗口使用明确m/z与Da/ppm，RT统一为秒；输出TIC/BPC/XIC和梯形积分（强度×秒），报告同一峰进入多个窗口。没有基线校正或浓度校准，不以XIC面积直接报告含量／纯度。未经质心化的轮廓谱会报错；RAW转换和复杂特征/加合分析可交现有解谱主技能和PyOpenMS工作流。

## 环肽：给定结构检查

```text
python scripts/assist.py cyclic examples/basics/cyclic.json --output cyclic.json --report cyclic.md
```

输入实际SMILES、`closure.kind`（amide/disulfide）和两个正整数原子map ID，可给预期分子式/形式电荷 `expected_formal_charge`。检查指定键是否存在、化学类型与是否在环，核对分子式和给定的游离N／C端酸前体闭环失水关系。SMILES形式电荷不自动等于永久固有电荷，固定锍等质量口径仍按诊断核心显式声明。闭环酸参照不能套到C端酰胺。图中的环酰胺不是自动证明N-to-C端环化；原子map身份仍须对照前体与谱证。

吸收 `peptide-agent/candidate-evaluation` 区分拓扑、结构、化学敏感性与计算评分的思路；该来源为非商业许可，本模块没有引入其代码、Pi运行器、模型或Rosetta/Amber/HADDOCK后端。二维图检查不等于环化成功、构象收敛或结合能验证。

## protocols.io：版本和步骤

```text
python scripts/assist.py protocols examples/basics/protocol.json --output protocol-summary.json
python scripts/assist.py protocols --fetch https://www.protocols.io/api/v4/protocols/12345
python scripts/assist.py protocols --fetch https://www.protocols.io/api/v4/protocols/12345 --execute --output protocol-summary.json
```

默认只读本地快照或输出联网请求计划。用户确需取得指定规程时可执行GET；使用命名环境变量 `PROTOCOLS_IO_ACCESS_TOKEN`，不读取整个环境、凭证文件或上传样品。没有token时返回缺项，不能声称已访问服务。仅允许官方HTTPS协议/步骤端点，禁止重定向；当前接口依据官方API文档和K-Dense来源核对。

输出标题、作者、DOI、version_uri、材料、原始步骤顺序与GUID链检查。HTML/Draft步骤只转换为文本，不执行其中的命令、链接或工具指示。没有版本或不支持的内容格式会明确标记。保留原JSON和完整来源，不用/latest静默替换旧版本。

## 共用知识与证据

来源索引为 [experimental-sources.json](experimental-sources.json)，书籍／网站后续去向见 [知识扩充计划](knowledge-expansion.md)。书目、摘要、选读段落、已核规则和当前计算分开；未读全书不能记作全书规则已经入库。新网页只按可追溯的事实摘要进入知识笔记，不自动变成精确质量变换。
