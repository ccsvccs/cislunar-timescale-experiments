# 参考论文与数据来源（地心时间星历项目）

> 版本与入库状态要明确。不要把官方链接、PDF 本体、转录文本和本地计算输出混为一谈。

## 1. Lu, Yang & Xie — Lunar Time Ephemeris LTE440

- 论文：*Lunar Time Ephemeris LTE440: definitions, algorithm and performance*
- 固定版本：[arXiv:2509.18511v1](https://arxiv.org/pdf/2509.18511v1)
- **已在 Git 中保存 PDF**：[sources/lte440_2509.18511v1.pdf](sources/lte440_2509.18511v1.pdf)
- **已在 Git 中保存提取文本**：[sources/lte440_v1.txt](sources/lte440_v1.txt)
- SHA256（原有本地审计）：`b6a620c0388f01b5349cbafdc2c7febc7ec8db43336b49dc65293ddaba2f72ac`
- 重点：Sec. 3 数值积分，半天分段/10 阶 Romberg、常数速率分离与恢复、地心 TT−TDB 对 DE440 的验证；Figure 1。
- 该文报告：分段步长实验差异约 10 fs，地心 TT−TDB 对照 JPL 在全时段小于 1 ps；这是数值复现一致性，**不等于** LTE440 月球模型的绝对物理精度为 1 ps。

## 2. Park, Folkner, Williams & Boggs (2021) — DE440/DE441

- 论文：*The JPL Planetary and Lunar Ephemerides DE440 and DE441*，*The Astronomical Journal* **161**, 105 (2021)
- DOI：[10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414)
- **JPL 官方完整 PDF**：[Park.2021.AJ.DE440.pdf](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf)（15 页，开放获取，CC BY 4.0；转载须保留作者、题目、期刊、DOI 归属）
- **已取得并验证实际 15 页 PDF**：[sources/park_2021_de440_de441.pdf](sources/park_2021_de440_de441.pdf)，SHA256 `96c3945deb12e0f0843c2bf0703d3542d950cf85ff903b1fe11a501634f68ee0`；见 [sources/park_2021_de440_de441.reference.md](sources/park_2021_de440_de441.reference.md)。
- 下载日期 2026-10-14；PDF 5,501,569 字节，`file`/`pdfinfo`/签名验证，页面 3–4 原式视觉核对；提取文本也已存档。
- 重点：Sec. 2.3、Eq. (3)–(7) 的 TT−TDB 模型（包含外部势、太阳 J2、矢量势、补充 PN 项）；Sec. 3 的轨道动力学；与 DE440/441 模型/数据适用范围相关内容。
- 勿把含站点位置项的完整表达式直接当成地心项；先置地心观测位置为零，再核对正负号、变量、单位、时间尺度缩放。

## 3. 原始研究基线

- Turyshev et al., *High-Precision Relativistic Time Scales for Cislunar Navigation*, 固定 [arXiv:2507.22145v2](https://arxiv.org/pdf/2507.22145v2)。
- 仓库存档：[sources/turyshev_2507.22145v2.pdf](sources/turyshev_2507.22145v2.pdf)；[model_audit.md](model_audit.md)；[tt_tdb_eq21_report.md](tt_tdb_eq21_report.md)。

## 4. 工程阶段

见 [phase2_jpl_time_ephemeris.md](phase2_jpl_time_ephemeris.md)。本阶段仅针对**地球质心 TT↔TDB 与 DE440 内置时间星历的一致性复现**，目标是在明确区间和未经调零/去趋势的原始差值上达到最大绝对残差 < 1 ps。此目标不是时钟固有时同步/绝对物理准确度认证。

**数据参考**：`data/de440.bsp` 提供几何轨道状态，`data/de440t.bsp` 的指定时间 SPK 段提供对照；常规 `de440.bsp` 不能假定包含时间段。详见 [data/README.md](../data/README.md)。
