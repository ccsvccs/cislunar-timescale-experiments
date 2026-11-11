# 固定 LTE440 发布产品：独立验收参考

本轮固定使用 Lu Xu（作者随包代码标注 Copyright © 2025, Lu Xu）发布的
[20250620 Release](https://github.com/xlucn/LTE440/releases/tag/20250620)，
资产为 [lte440-tclmtdb.zip](https://github.com/xlucn/LTE440/releases/download/20250620/lte440-tclmtdb.zip)。
该发布早于 Lu、Yang、Xie 的 [arXiv:2509.18511v1](https://arxiv.org/pdf/2509.18511v1)，
不能把论文描述的每个生产细节无条件归于此固定产品。

2026-11-11 从已下载的 `/tmp/lte440-reference-20250620.zip` 实际验证 SHA 并解压：

| 文件 | 字节数 | SHA256 |
|---|---:|---|
| 原 ZIP | 10072926 | `baa74b47ef42633297617a1ed507daeb2157eb7a9a53bda26b5d5cc951b7b6dc` |
| `data/lte440.bsp` | 35361792 | `f48c36b0c1b872f51616790db753538b4895bc2ceae73386a72eb0150e44f1e8` |
| `data/lte440.tpc` | 278 | `aa897dd4aeba15dfeb6f80883ad6b3e2d7c4f098690b02cdeb021167f5734c92` |

大型 BSP 不提交 Git；以固定 ZIP、路径和哈希复现。没有改写作者的 BSP/TPC，
没有把参考值送入核心积分、零点校准或修正表。`lunar_model.py` 不导入参考模块。

## 产品读取规则

BSP 只有一个 SPK Type 2 段：center `1000000000`、target `1000000005`、
frame `1`，TDB JD 覆盖 **2287184.5–2688976.5**。X 是以秒表示的 TCL−TDB
周期部分，Y/Z 不作物理位置解释；不作 km→m 换算。TPC 明确给出：

```text
BODY1000000005_RATE = +6.798355238724282D-10
```

完整发布函数由作者 ZIP 内 `lte440.py`/`lte440_2.py` 的实际运算确定：

\[
Q_{\rm ref}(s)=X_{\rm BSP}(s)+R_{\rm TPC}(s-s_0),\qquad
Q'_{\rm ref}=\frac{dX_{\rm BSP}}{ds}+R_{\rm TPC},
\]

其中 `s0=T0+TDB0`，`T0=JD2443144.5003725`，`TDB0=−65.5e−6 s`。
jplephem 的解析导数为每 JD 日，实际除以 86400 转成每秒。输入两部分 JD
用于状态和相对原点的间隔，线性恢复使用 NumPy longdouble，避免几十秒函数
之差引入不必要的舍入。没有取代发布的 TPC rate，尤其没有用论文式(22)的
LTE441 校正率替换原发布值。

`src/tt_tdb/lte_reference.py` 默认拒绝 SHA 不匹配的 BSP/TPC、错误的段身份、
缺失/多重 TPC rate 以及覆盖外时刻。参考输入与模型输入的哈希分别进入结果 JSON。

## 已执行独立验证

`tests/test_lte_reference.py` 核查：

- 作者随包样例 1600.0/J2000.0/2200.0 的完整函数；J2000.0 TDB 为
  `TCL−TDB=+0.49330749643254862 s`。
- 原始月心锚点给 `TCL−TDB≈+65.5 μs`，与定义相差约 3.6e−17 s；
  本项目模型仍按定义固定 +65.5 μs，不对齐这个参考舍入差。
- 实际 SPK 多项式记录边界、覆盖两个端点，与独立 CSPICE 的 X/VX 比较；
  保持秒和 /s 单位，不混入距离换算。
- 缺数据、错版本、覆盖外时刻及不配套 TPC 都显式报错。

作者未随此 Release 发布生产积分代码或输入轨道清单及哈希。因此本项目相对该产品
的原始 residual 只说明数值函数一致程度，不能据其小于某阈值证明原始 DE440
小体状态已取得或月球物理时间绝对准确度达到 ps。相关模型限制与实际结果见
`docs/full_span_lunar_report.md`。
