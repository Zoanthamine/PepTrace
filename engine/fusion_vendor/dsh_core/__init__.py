"""DSH 版多肽副产物诊断引擎（vnext_dsh）。

与项目内既有的 `vnext`（另一套实现）并列、互不写入：本包只**读**知识库资产
（规则库_总.json、SUPPLEMENTS.json、多肽质量组件.json），把全部产物写在 vnext_dsh/ 下。
"""

__version__ = "0.1.0"
__all__ = [
    "formula",
    "components",
    "state",
    "schema",
    "rules",
    "matching",
    "candidates",
    "report",
]