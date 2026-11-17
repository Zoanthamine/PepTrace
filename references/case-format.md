# 融合版病例格式

Agent将用户原文写入`input_text`，整理下列内部JSON。保留原始观测和假设；不要要求用户自行编写JSON。

```json
{
  "case_id": "sample-1",
  "input_text": "用户原话",
  "intended_product": {"sequence": "MELK", "nterm": "free", "cterm": "acid", "modifications": []},
  "reaction": {"description": "DMSO处理", "stage_tags": ["dmso"], "stage_tags_absent": []},
  "observations": [{"kind": "delta_mono", "value": 15.994915, "reference": "intended_product", "tolerance": 0.001, "tolerance_unit": "Da"}],
  "max_events": 2,
  "max_nodes": 100000
}
```

端基为`nterm=free/acetyl`和`cterm=acid/amide/hydrazide`等已登记状态。位点1起计；端基位点`N`、`C`。常用修饰示例`{"site":2,"group":"Carbamidomethyl"}`；CAM净增57.021464 Da，Acm净增71.037114 Da，两者不同。保护基和已知局部净变化使用登记组件；不同原子同位点确需区分时可补`atom`。

多链写`{"chains":[{"id":"a","sequence":"CG","nterm":"free","cterm":"acid"},{"id":"b","sequence":"GC","nterm":"free","cterm":"acid"}],"bonds":[{"kind":"disulfide","endpoints":[{"chain":"a","site":1},{"chain":"b","site":2}]}]}`。链内桥也可写单链`bonds:[{"kind":"disulfide","sites":[1,3]}]`。连接不明确时保留缺项，不按普通线性肽代算。

观测：

- 中性单同位素：`kind=neutral_mono`；平均质量：`kind=neutral_average`；容差`Da`或`ppm`。
- `delta_mono`和`delta_approx`允许负数或零，误差使用`Da`；`reference`可为`intended_product`或`starting_material`，后者必须有起始结构。
- m/z：`{"kind":"mz","value":500.1,"charge":2,"ions":{"H":2},"tolerance":5,"tolerance_unit":"ppm"}`。显式固定电荷用`intrinsic_charge`，裸固定+1离子用`charge=1,intrinsic_charge=1,ions={}`，不另加质子。
- 未知离子可明确给`ion_hypotheses:[{"charge":1,"ions":{"H":1}},{"charge":1,"ions":{"Na":1}},{"charge":2,"ions":{"H":2}}]`。这是有限假设范围。
- 双峰：`{"kind":"peak_pair","reference":"intended_product","peaks":[观测1,观测2],"error_model":"bounded"}`。两峰各自有电荷／加合和容差。默认传播边界误差；明确独立高斯标准不确定度时可用`independent_gaussian`。

若要按带电片段计算，链上给`fixed_charge_sites:[{"site":"N","atom":"Nalpha","charge":1}]`或侧链具体位置。未标明位置的固定电荷不会在断链时被猜配到某片段。

条件使用明确已有事实。例如`piperidine`、`acid_exposure`、`acid_global_deprotection`、`edt`、`dmso`、`formaldehyde`、`formic_acid`。只提“用了TFA”可支持酸接触，不必推断已经完成全局脱保护。未知历史保留；确证未使用的条件写`stage_tags_absent`。必要位点历史可放`context.features_present/features_absent`。

结果中的`candidates`、`condition_conflicts`和`source_leads`分开解释。`complete_for_practical_mass_scope`按本次预算和支持的模型说明完整性；`full_book_capabilities_complete=false`表示没有声称所有复杂书源化学都已精确执行。

额外氨基酰偶联使用`reaction.acyl_donors:[{"sequence":"F","nterm":"free","cterm":"acid"}]`，此处供体是最终脱保护的Phe；若供体仍带Fmoc，必须保留其保护基。N端与Lys侧链酰化以显式酰胺连接分别枚举，不仅仅加147。

三硫形成中，DCC路径为`cystine_coupling+dcc`，BOP/HBTU分支要求相应偶联剂；明确未加DCC写`stage_tags_absent:["dcc"]`。此阴性事实不能否定所有硫插入终点。MTBE仅支持`mtbe`，不自动推定`strong_acid_workup`。

显式旧DSH JSON可经适配入口读取；端基、容差、电荷、加合均须明确。含双峰、固定电荷、旧tags/stage或二硫修饰请整理为上面的状态格式，避免继承旧默认值。`dsh_source_review`仅为独立来源反查。
