# MS2区分辅助

技能目录的`engine/tools/compare_msms.py`使用PyOpenMS，比较标准线性肽的完整b/y片段和局部元素质量移位。需要至少两个标准残基；未定义连接和非局部变化不作为普通残基移位处理。

内部输入包含`candidates`（候选名到序列／端基／局部修饰的对象）、`peaks`（m/z、强度二列列表）、`tolerance`、`tolerance_unit`和可选`max_fragment_charge`。

```text
python engine/tools/compare_msms.py ms2-case.json --output 碎片比较.json
```

也可用`--mzml`与`--native-id`指定一个质心MS2谱，此时输入不同时放`peaks`。结果列出吻合、未吻合及仅在受测候选中专属的片段。用作区分证据，不作为结构确证或谱图准确率声明；真实样品归属需另有独立依据。
