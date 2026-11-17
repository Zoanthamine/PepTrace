# 肽迹基础实验辅助

模块：hplc_basics。

review and scaling arithmetic; no validated retention or resolution prediction

```json
{
  "module": "hplc_basics",
  "purpose": "LC-MS",
  "geometric_column_volume_ml": 2.4928537706235003,
  "flow_ml_min": 1,
  "gradient": [
    {
      "time_min": 0,
      "percent_B": 5
    },
    {
      "time_min": 30,
      "percent_B": 65
    }
  ],
  "gradient_slopes_percent_B_per_min": [
    2.0
  ],
  "transfer": {
    "suggested_flow_ml_min": 0.20841209829867677,
    "time_scale": 0.6666666666666667,
    "gradient": [
      {
        "time_min": 0.0,
        "percent_B": 5
      },
      {
        "time_min": 20.000000000000004,
        "percent_B": 65
      }
    ],
    "basis": "constant linear velocity and gradient column volumes; dwell volume, particle size, viscosity, pressure and selectivity still need verification"
  },
  "findings": [
    {
      "name": "TFA/ESI兼容性",
      "reason": "LC分离需要与ESI响应分别评价；离子对添加剂可能抑制响应。",
      "source_id": "HPLC-BASICS"
    },
    {
      "name": "酸/甲醇甲酯化线索",
      "reason": "检查D/E侧链及游离C端的处理副反应；不能把方法引入的峰当原始样品成分。",
      "source_id": "YANG-2016"
    }
  ],
  "checklist": [
    "实际柱化学/孔径/耐压",
    "样品溶剂和进样负载",
    "关键杂质对分辨率",
    "MS添加剂兼容性",
    "清洗/再平衡与空白",
    "检测波长、纯度与回收率口径"
  ],
  "scope": "review and scaling arithmetic; no validated retention or resolution prediction",
  "peptrace_version": "1.1.0",
  "source_input_sha256": "48d9b3ff2e357c03df01194818434446761a1009bb4813012410017390dbdc6b",
  "source_notes": [
    {
      "id": "YANG-2016",
      "title": "Side Reactions in Peptide Synthesis",
      "year": 2016,
      "authors": [
        "Yi Yang"
      ],
      "type": "book",
      "read_status": "existing canonical rule evidence",
      "role": "shared rule/mass/site reference; current 450 records and 308 compiled branches"
    },
    {
      "id": "HPLC-METHOD",
      "title": "develop-hplc-method",
      "type": "skill_method",
      "url": "https://github.com/pjt222/agent-almanac/blob/8ff93da830524d5ffcf2b894acf54a1b3a9e8c64/skills/develop-hplc-method/SKILL.md",
      "revision": "8ff93da830524d5ffcf2b894acf54a1b3a9e8c64",
      "license": "MIT",
      "read_status": "method-development guidance read",
      "reuse": "workflow attribution; independently written parameter checks and column-volume arithmetic"
    },
    {
      "id": "HPLC-BASICS",
      "title": "Peptide purification: scaling analytical to preparative HPLC",
      "type": "vendor_webpage",
      "url": "https://www.peptide.com/resources/peptide-purification-scaling-from-analytical-to-preperative-hplc/",
      "read_status": "visible webpage text read",
      "reuse": "factual summary and dimensional arithmetic; not historical bubble/recipe instructions"
    }
  ]
}
```
