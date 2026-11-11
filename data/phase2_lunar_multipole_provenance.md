# 月球静态二阶外势：来源、姿态与锚定积分

审计与实算日期：2026-11-10。该诊断独立于默认 TT–TDB 模型，输出符号为
`Δ = TDB − TT` 的可加增量。没有使用参考时间数据构造 integrand、选轴或调参。
它测量静态月球形状项在地心的直接外势，不宣称重现 JPL 原始时间积分的缺项。

## 1. 原始来源与实际文件

- Park, Folkner, Williams & Boggs (2021), *The JPL Planetary and Lunar Ephemerides
  DE440 and DE441*, AJ 161:105, DOI
  [10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414)。
  [JPL 原始 PDF](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf)，本地
  `docs/sources/park_2021_de440_de441.pdf`。
  Eq. (8)–(11) 与 Sec. 2.4 定义 PA 姿态；p6 Eq. (28) 说明未归一化球谐约定；
  Eq. (30)–(34) 给出系数与惯量的关系；p8 Eq. (47)–(54) 定义未变形惯量及变形项。
  本轮实际视觉核对 p6 公式与正文。论文为 CC BY 4.0，保留作者、题名和 DOI 归属。
- [JPL DE440 ASCII header](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440)，
  本地固定副本 `docs/sources/de440_header.440`。从 GROUP 1040/1041 配对解析
  `AM,J2M,C21M,S21M,C22M,S22M,LBET,LGAM`，没有从 Turyshev 表格另取数值。
- [NAIF DE440 PA binary PCK](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc)，
  本地 `data/moon_pa_de440_200625.bpc`，12,863,488 bytes。官方原始注释
  [moon_pa_de440_200625.cmt](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.cmt)
  保存于 `docs/sources/moon_pa_de440_200625.cmt`。注释标注 2021-06-25 制作、
  来源为 DE440，并明确这是 PA 相对 ICRF 的真实姿态。
- [NAIF lunar FK 2025-04-16](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf)，
  本地 `data/moon_de440_250416.tf`，19,478 bytes。明确采用 `MOON_PA_DE440`，
  frame/class ID 31008、center 301、class 2；不用通用别名或 `IAU_MOON`。
  此 FK 仍指向相同 DE440 PCK，包含官方独立几何示例。

| 文件 | SHA256 |
|---|---|
| `moon_pa_de440_200625.bpc` | `60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f` |
| `moon_de440_250416.tf` | `a47c71e9c9f33796bdafb2c9d69a7ee447b6016ecad80f71cd6f3e479f9cf768` |
| `moon_pa_de440_200625.cmt` | `97d3bb1e75e7807a5e0b76f836a2d1d192e0d5bde0c0a87979b4c13255aef8c1` |
| `de440_header.440` | `0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359` |
| Park 2021 PDF | `96c3945deb12e0f0843c2bf0703d3542d950cf85ff903b1fe11a501634f68ee0` |

2026-11-10 下载经过环境网络审批，以 Python `urllib.request.urlopen` 从上述
三个 NAIF HTTPS 地址原样下载；核验 PCK 开头为 `DAF/PCK `，FK 为 `KPL/FK`，
实际用 jplephem/CSPICE 解码，不把 HTTP 成功当作内容有效性。PCK 二进制不应提交 Git；
这里记录下载地址与哈希以供重取。本轮没有下载大型轨道核，已有 `de440.bsp`
和 `gm_de440.tpc` 的实际 SHA 记录在结果 JSON。

## 2. 系数、归一化与适用范围

| header 名 | 原值 | 本轮用途 |
|---|---:|---|
| AM | 1738.0 km | 与 header 系数配套参考半径 |
| J2M | 2.03214360135e−4 | 主诊断静态 J2；未归一化 C20=−J2 |
| C22M | 2.2380850528171017e−5 | 主诊断静态 C22 |
| C21M | 2.277956919447986e−10 | 仅字面静态敏感度通道 |
| S21M | 1.334345035839262e−9 | 仅字面静态敏感度通道 |
| S22M | 0 | 仅字面静态敏感度通道 |

系数不作归一化变换。用 Park Eq. (47)–(51) 的 J2M/LBET/LGAM 重建
未变形惯量，所得 `(B−A)/(4mR²)` 与 C22M 相差 6.83e−18，支持主诊断使用
未变形 J2/C22；这一检查不是用参考时间数据校准参数。

**范围限制：** Park Eq. (30)–(34)、(54) 明确完整二阶系数随潮汐和转动变形而变化，
不能把 header 的静态数值冒称完整瞬时重力场。公开 header 中非零 C21M/S21M 的
确切静态物理角色未由本轮来源完全解释，因此主结果只采用有上述惯量对应的 J2+C22；
其余项独立列为“把 header 数值原样固定”的敏感度，不混入主结果。
没有猜测潮汐系数的时间规律、延迟姿态或平均修正。

几何和 GM 使用同一 DE440 的 TDB-compatible SI 输入，`GM_M=4902800118457.55 m³/s²`。
AM 按 header 的 km 直接换算为 m；若其长度约定存在一个 L_B 量级的缩放差，
本次 R² 因子造成约 3.1e−8 的相对变化，即累计结果约 2.4e−8 ps，远小于本次目标；
本轮仍明确保留此半径尺度约定，不把它混称完全核实的物理半径测量。

## 3. 地心外势与真实姿态

令 n 为月心指向地心、表达在 DE440 PA 帧的单位向量，r 为两中心距离：

\[
U_2=\frac{\mu_M R_M^2}{r^3}\left[
-\frac{J_2}{2}(3n_z^2-1)+3C_{21}n_xn_z+3S_{21}n_yn_z
+3C_{22}(n_x^2-n_y^2)+6S_{22}n_xn_y\right].
\]

这使用未归一化、无 Condon–Shortley 相位的地学约定，C20=−J2。
符号和系数也可直接由 Park Eq. (30)–(34) 的惯量关系以及
`U2=G[Tr(I)−3n·I·n]/(2r³)` 恢复，不依赖 Turyshev 的 C20/J2 记号。
`lunar_multipole.py` 按五项分别输出 `U2/c²`；积分结束后统一乘
`A=(1−L_G)/(1−L_B)=1.0000000148082688`，正向加入 TDB−TT。

Earth−Moon 向量直接用两个真实 EMB-relative SPK 段 `(3,399)−(3,301)`，
减少公共 AU 大数相减。二阶偶数势在反向 r 上不变，但实际方向仍按月心指地心定义。
PCK type-2 Euler 角矢量化解码，使用被动旋转
`R3(psi) R1(theta) R3(phi)` 从 ICRF/J2000 到 DE440 PA，
不引入同步自转、固定轨道法线、ME 或 `IAU_MOON` 替代姿态。

PCK 真实 descriptor 覆盖 JD 2287184.5–2688976.5，
TDB J2000 seconds `[-14200747200,20514081600]`，内部有两个连续段；
本轮积分位于覆盖内。130 个独立姿态检查遍布 1977–2026 并含 FK 官方示例：

- 矢量化矩阵与 `spice.pxform('J2000','MOON_PA_DE440',t)` 最大元素差 2.85948e−12；
- 正交误差与 determinant−1 最大均为 4.44089e−16；
- 地月向量与独立 `spice.spkpos` 最大差 0.00116051 m；
- 地心方向 PA 的 n_x 在 0.9845586–0.9998448，实际天平动明显非零；
- FK 示例 2022-09-30 TDB (`t=717768000`) 的 PA 向量为
  `(373997.028323,−23558.987290,10284.057074) km`，与官方示例的 1 m 舍入吻合。

这些是输入/数值一致性检查，不是重力参数或 DE440 绝对姿态的物理误差评估。

## 4. 已执行积分与实测结果

命令：

```sh
.venv/bin/python -m pytest -q tests/test_lunar_multipole.py
.venv/bin/python scripts/compare_phase2_lunar_multipole.py
```

新测试 **5 passed in 0.11s**。覆盖 J2 极点/赤道、C22 长短主轴、m=1/S22 的
解析符号和系数、外势域检查、header 惯量一致性、真实帧正交与方向、官方示例、
跨 PCK 段边界，以及缺核时报错。原默认模型没有修改。

锚点固定为 `T0 + TDB0`，即 TDB J2000 seconds `−725803167.8160655`，
初始所有新通道贡献为零；不是对照参考数据调零。连续积分至 2026-01-01，
输出 2024-01-01–2026-01-01 共 732 个每日 TDB 历元，以及锚点和
1978–2026 每年 1 月 1 日共 50 个历史检查点。

| 主结果 | 累计加入 TDB−TT |
|---|---:|
| 2024-01-01，J2 | +0.430303969169 ps |
| 2024-01-01，C22 | +0.284630635305 ps |
| 2024-01-01，J2+C22 | **+0.714934604474 ps** |
| 2026-01-01，J2 | +0.448629953754 ps |
| 2026-01-01，C22 | +0.296746688541 ps |
| 2026-01-01，J2+C22 | **+0.745376642295 ps** |

现代主增量 RMS 为 0.730212403394 ps，最大绝对瞬时速率为 6.09574666614e−22。
敏感度 C21/S21 的现代最大累计绝对值分别为 1.30813e−9 ps、2.21075e−9 ps；
S22 因 header 原值为零而恒为零。五通道字面静态总和与主结果几乎一致，
但两者的来源解释仍保持区别。

主求积采用半日步长 GL24，35,793 个完整分段、1,316,556 次节点求值；
独立半日 Romberg8 共 4,717,659 次求值，另用日步长 GL24 共 672,264 次求值。
所有输出历元主和的最大差为：

| 数值检查 | 最大绝对差 |
|---|---:|
| 半日 GL24 − 半日 Romberg8 | 1.06022905662e−27 s |
| 半日 GL24 − 日步 GL24 | 1.00974195868e−28 s |

最终源码版本这一轮总计算用时 8.612154932 s（不含绘图保存）；各方法壁时分别
1.594096497、5.988668304、0.923858002 s。独立求积一致仅验证数值积分，
不能证明静态系数足以表示完整物理势。统计 helper 另列斜率与去趋势 RMS
作为诊断，原始输出和所有主指标没有施加拟合修正。

结果：`results/tables/phase2_lunar_multipole_diagnostic.json`、
`phase2_lunar_multipole_all_epochs.csv`、`phase2_lunar_multipole_modern.csv`，
图为 `results/figures/phase2_lunar_multipole_diagnostic.png`。
JSON 内包括所有输入核与源码依赖 SHA、C_LIGHT/A、固定历元、版本、方法及实际耗时。
没有叠加旧 SB441 残差，没有改变既有验收结论；主 session 可做明确标记的后处理。

## 5. 尚未解决的物理项

完整潮汐/转动变形的二阶重力系数、degree≥3，以及多极势进入 c⁻4 标量/向量/混合项
未实现；本轮未判断 JPL 原时间星历是否含上述直接月球形状项。不能据新增增量的大小
或某个 residual 的改善反推这些项应该进入复现模型。没有做 1100 年累计验收。
