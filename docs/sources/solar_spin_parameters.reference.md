# 太阳自转直接势：参数与实验来源

审计日期：2026-11-10。这是一项单独标记的物理诊断，尚未确认为 DE440
时间星历的组成项。Phase 1 模型与 Park 候选模型均未改动。参数、轴方向、
零点和趋势均未依据 JPL residual 选择。

## 原始公式与固定参数

- [Turyshev arXiv:2507.22145v2](https://arxiv.org/pdf/2507.22145v2),
  存档为 `turyshev_2507.22145v2.pdf`：p35 (A5)、(A8)、(A9)，p40
  (A33)、(A34)、(A39)。已视觉核查 PDF p35/p40。采用的角动量为 p40 的
  `S_sun = 1.8838e41 kg m^2/s`，该文将其归于 Ref.45；该段未给不确定度。
- [Park et al. 2021](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf),
  存档为 `park_2021_de440_de441.pdf`：**Sec.3.4**，PDF p7，Eq.(39)
  后给出太阳自转北极 `RA=286.13 deg`、`Dec=63.87 deg`，引用 Archinal
  et al. 2018。本诊断将该极轴固定于 J2000/ICRF。项目既有资料将该极轴
  指向 Sec.3.5，节号有误；既有文件保持原样。未以 Park 的惯量/转速参数
  重建的值替换本次采用的 Turyshev 角动量。
- [NIST SP961, CODATA 2022, May 2024](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
  单页表给出 `G=6.67430e-11 m^3 kg^-1 s^-2`，标准不确定度为
  `1.5e-15 m^3 kg^-1 s^-2`。审计当日通过浏览器 PDF 解析器读取了该表及
  括号不确定度的说明。这些是测量常数，不是精确定义值。尝试直接下载存档时，
  NIST 不带 www 的主机返回 HTTP403，带 www 的主机返回 HTTP502；
  未保存 NIST PDF，也未为其虚构 PDF 哈希。
- [Cang et al. 2016 出版社摘要](https://www.hanspub.org/journal/paperinformation?paperID=17615),
  DOI [10.12677/AAS.2016.42004](https://doi.org/10.12677/AAS.2016.42004):
  审计当日随后核查了原始出版社摘要。摘要说明，相同太阳角动量估值所用的
  太阳惯量由 Eggleton 恒星演化程序计算。只核查了摘要/元数据，未核查全文
  或不确定度分析。这是恒星结构估计，不能称为 DE440 拟合的自转参数。
  此项后续核查记录于此，未追溯修改实验中已经计算哈希的源码。

Turyshev 和 Park 原始 PDF 的 SHA256 保存在实验 JSON 的
`model_signature.paper_sha256` 与 `axis_source_sha256` 中。上述两项网页核查
通过浏览器读取，未归档原始字节快照；本 Markdown 记录的 SHA 不能称为这些
原始网页/PDF 来源的哈希。

## 标量势、矢量势与尺度约定

令 `r = x_E - x_sun`，原式给出

\[
\delta\mathbf w=-\frac{\mathbf r\times(G\mathbf S)}{2r^3},\qquad
\delta D=\frac{2\mathbf v_\odot\cdot[\mathbf r\times(G\mathbf S)]}{r^3}.
\]

各被积函数通道为

\[
\delta F_{vector}=-\frac{4\mathbf v_E\cdot\delta\mathbf w}{c^4},
\qquad \delta F_{scalar}=-\frac{\delta D}{c^4},
\qquad
\delta F_{total}=-\frac{2(G\mathbf S)\cdot[\mathbf r\times(\mathbf v_E-\mathbf v_\odot)]}{c^4r^3}.
\]

最后一式也单独计算并积分，用于检验其与标量/矢量通道之和的闭合性。按本项目
正负号约定，增加的速率为 `d(delta Delta)/dTDB = A_SCALE delta F_total`，
增加的时间差为该速率自原始 `T0_TDB` 锚点起的积分。锚点处修正严格为零，
不再次加入 TDB0。

将文献中的 SI 角动量解释为未缩放的物理参数。令 `k=1-L_B`、`t'=k t+const`、
`x'=k x`，则速度不变；`GS/r^2` 不变要求 `(GS)_TDB=k^2 G_SI S_SI`。
这是明确声明的尺度解释，并不声称原始恒星计算指定了 TCB/TDB 尺度。
独立的 unscaled-GS 通道用于衡量该解释的敏感性。未重新积分地球轨道：
既有 LT 轨道效应与新增的时间直接势是不同项。

## 复现方法与实测结果

命令：

```sh
.venv/bin/python -m pytest -q tests/test_solar_spin.py
.venv/bin/python scripts/compare_phase2_solar_spin.py > results/tables/phase2_solar_spin_run.log 2>&1
```

专项测试结果：5 passed。测试包括顺行/逆行及自转符号、径向运动、标量/矢量
闭合与共同速度抵消、一致的自转/坐标缩放、旋转不变性、拒绝无效输入，以及
真实 DE440 固定锚点案例。未修改主要物理来源。主 session 的全套验证另行记录。

正式计算采用 GL24/0.5day、Romberg6/0.5day 和 GL16/1day，三者均从
`T0_TDB` 独立积分至 2026-01-01。报告记录绘图前墙钟耗时 31.77s、
方法统计、软件版本及代码/数据 SHA256 身份。结果采样包括 2024-01-01
至 2026-01-01 的 732 个每日 TDB 历元，以及精确锚点和五个历史历元。
另取 17,898 个每日采样点描述实测速率包络；这些采样不能确立连续最大值或
物理不确定度上界。

| 量 | 实测诊断值 |
|---|---:|
| 2024-01-01 合并锚点修正 | -6.097405777979 ps |
| 2026-01-01 合并锚点修正 | -6.357068818062 ps |
| 2026-01-01 矢量部分 | -6.357053467970 ps |
| 2026-01-01 标量部分 | -0.000015350092 ps |
| 现代区间合并速率最小/最大值 | -4.323605154485e-21 / -3.909872201820e-21 |
| GL24/0.5day 与各独立方法的最大自转通道差 | 8.077935669463e-28 s |
| 标量+矢量积分与合并积分的闭合差 | 1.382811537094e-27 s |
| 合并积分与相对速度积分的闭合差 | 8.077935669463e-28 s |
| GS 尺度解释的最大敏感性 | 1.971352216226e-19 s |
| 仅 G 标准不确定度的敏感性 | 1.428704617277e-16 s |

很小的求积差与尺度差不能当作总物理误差上界。本项未评估来源角动量、固定太阳
极轴的不确定度及遗漏的其它外部天体自转。

将修正叠加到前次实测的混合 SB441 诊断之前，脚本核查当前 Park 完整模型指纹、
存档源码/数据 SHA256 身份、原始锚点、全部 732 个历元、主要 CSV 的来源，
以及独立 JPL 值/导数。还从相同锚点重新积分两种太阳 J2 轴的差，以验证既有
轴差表；没有改动该表。

| 前次实测诊断 | 原始最大绝对 residual | 加入太阳自转直接势 |
|---|---:|---:|
| SB441，太阳赤道 J2 轴 | 2.818949568673 ps | 9.176018386735 ps |
| SB441，Park 字面黄道 J2 轴 | 1.953519689910 ps | 8.308582589430 ps |

这些叠加结果说明，按独立来源定义的自转项**不能**修复既有约 2ps 的 DE440
residual。它们不能证明原始 DE440 时间积分器采用了哪些项，也不能确立所得
候选模型的绝对物理误差。

机器可读证据：`results/tables/phase2_solar_spin_diagnostic.json`、
`phase2_solar_spin_time_series.csv`, `phase2_solar_spin_history.csv`,
`phase2_solar_spin_overlay.csv`, and `phase2_solar_spin_run.log`.
图：`results/figures/phase2_solar_spin_history.png` 和
`results/figures/phase2_solar_spin_modern.png`.
