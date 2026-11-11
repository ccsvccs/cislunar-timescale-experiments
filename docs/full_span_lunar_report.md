# 月心 TCL−TDB：固定 LTE440 产品的全覆盖比较

本报告的比较量为同一 TDB 标签下月球质心处的 `Q=TCL−TDB`，单位为 SI 秒。
没有引入 TL 或具体月面/航天器钟的固有时。参考是固定 20250620 发布产品，
不是 DE440 内的 TT−TDB 段，也不是把论文中经 LTE441 校正的平均率替换进旧产品。

本轮实际结果与数值检查完成后，第 5 节结果已由实际输出补全。公式、输入与参数选择
在全区间运行前已冻结；任何原始残差均未用于选择常数、轴方向或修正趋势。

## 1. 原点、尺度与数值恢复

采用 IAU 月心原点：在 TCB=T0 时 TCL=T0，其中
`T0=JD2443144.5003725`。TDB 输入原点为 `s0=T0+TDB0`，
`TDB0=−65.5 μs`。因此模型严格固定 `Q(s0)=+65.5 μs`。

令 `k=1−L_B`，`L_B=1.550519768e−8`。以 TDB-compatible 的位置、速度、GM
计算月心速率函数 F，则

\[
Q'(s)=\frac{L_B-F_M(s)}{k},\qquad
Q(s)=-\mathrm{TDB}_0+\frac{L_B}{k}(s-s_0)
-\frac1k\int_{s_0}^s F_M(u)\,du.
\]

完整 `1/k` 覆盖 c⁻²、c⁻⁴ 两组，不使用 Turyshev 式(44)中省略比例的周期近似。
为减少大线性项抵消，预先固定纯计算常数 `F*=1.4825e−8`，实际积分 `F−F*`，
再按同一个浮点常数精确恢复：

\[
Q(s)=-\mathrm{TDB}_0+
\frac{(L_B-F_*)(s-s_0)-\int_{s_0}^s[F_M(u)-F_*]du}{k}.
\]

F* 不是来自 LTE BSP/TPC、拟合残差或重新估计平均率。恢复使用 longdouble，
锚点左右分别累计，禁止先从区间左端累计再减两个大总和来构造原点。

## 2. 明确实现的物理项

主要体的表达复用 Park et al. (2021) 式(3)–(7) 的公开代数，在月心重新计算：

\[
F_M=\frac{v_M^2/2+w_0+w_l}{c^2}
+\frac{v_M^4/8+3v_M^2w_0/2-w_0^2/2-4\boldsymbol v_M\cdot\boldsymbol w-\Delta_M}{c^4}.
\]

本轮 c⁻⁴ 中的 w0、vector、Δ 取主天体集合；SB441 与环只另外加入 c⁻² 直接势，
尚未并入 c⁻⁴ 的混合、平方、向量和 PN 项，不能把同一记号误读为已实现全部小体修正。

外部体排除 Moon、加入 Earth；其余为 Sun、Mercury、Venus 和 Mars–Pluto
系统质心。Moon GM 仍进入 PN 内层体间势，不能把外势排除自势误用为内层也排除 Moon。
真实 DE440 SPK Type 2 多项式解析得到重心位置、速度和加速度，未用有限差分或玩具轨道。

- leading 标量外势包含上述主天体、真实 SB441 的 343 asteroids/30 KBO，以及公开
  44 AU/36 等质量点环的独立静态几何势。
- 主天体的向量势与 PN 速度、体间势、径向速度、加速度分量分别计算。
- Sun J2 使用现有模型预先固定的太阳自转轴（RA=286.13°、Dec=63.87°），
  半径、J2 来自同一 DE440 header；没有根据 LTE 残差挑选另一个轴。
- Earth 外部 J2 在月心计算，参数为同一 header 的 `J2E=0.00108262539`、
  `RE=6378.1366 km`；Moon 自身重力势没有加入月心坐标变换。
- 太阳直接自旋项按本轮范围排除。Moon 固有时、自势、月面站几何不进入 Q。

实现位于 `src/tt_tdb/lunar_model.py`。原 geocentric 核心文件没有修改。
Moon 观察中心兼容适配器显式把外 Earth 命名为 `earth_external`，避免原地心模块
禁止 Earth 自势的 guard；这些旧字段名只供代数接口复用，不改变本轮观察中心。

## 3. Earth J2 轴与尚不能宣称完整的部分

轴源为 Park (2021) Sec.2.5 式(20)–(26)：Vondrák 等 (2011，2012更正)
长期岁差，叠加明确的 18.6 年章动项。独立实现两组 pole 系数；ERFA 只在测试中
作公式核对，未提供参考时间或核心速率输入。date 基轴以 inertial 分量为行的矩阵，
转回 inertial 使用转置；这一方向经独立 reviewer 与 ERFA/CSPICE 旋转核查。

主轴还包含 DE440 header `ROTEX=+0.005343991604414805`、
`ROTEY=−0.01712882484031753`，按 arcsec 解释、按 Park 式(25)
`Rx(−Φx) Ry(−Φy)` 顺序作用。单位证据链是 JPL Folkner et al. (2014)
p50 Table10 明确列出 Φx/Φy 的 arcsec 单位，并与官方 DE430 header 同字段数值
逐字对应；将同字段单位约定用于 DE440 是**明确记录的版本推断**，不是见到了原始
DE440 极轴生产代码。主参数选择在读本轮最终残差前固定。

另用第 11 个通道单独积分“零 offset J2 势−主型 J2 势”，只作轴约定敏感度；
`rate` 和 `restore` 明确只取前 10 个物理通道，诊断通道不会悄悄改变主结果。

Earth J2EDOT/J2ET2 的时间规律和参考历元未由本轮一级资料核清，因此保持
**constant J2**，没有从字段名猜时间律。高阶 Earth 球谐、潮汐、原 DE440 环动力学、
SB441/环的 PN 与向量项未完整取得。SB441 是以 DE441 行星另行拟合/积分的产品，
不能称作原始 DE440 小体状态。上述限制不会因残差较小而自动消失。

## 4. 发布参考、覆盖与验证

固定文件身份、SHA、读取规则及作者样例见
[LTE440 来源记录](../data/lte440_reference_provenance.md)。reader 默认要求固定
BSP/TPC 哈希，X 为秒、解析导数由 /day 转为 /s，并恢复 TPC 的完整线性项。
参考锚点微小插值/舍入差保持原样，没有被拿来调本项目原点。

实际共同覆盖为 TDB JD **2287184.5–2688976.5**。每天采样、保留两个覆盖端点，
再额外加入精确 1977 原点；完整每日表采用 gzip 文本，数字以 20 位科学计数输出，
保留至少 float64 的往返精度。原始 residual 在扩展精度相减后另列，验收直接用该列，
不从导出的两个大 Q 值相减重建。401,793 个日点加原点共 401,794 个时刻；
这不是连续全域误差上界。

主要体、PN、两个 J2 与环按日 GL8 独立积分；SB441 leading 势按日 GL4。
另将前者完整重算为日 Romberg5、后者完整重算为日 Romberg4。全覆盖共同点直接
比较两条完整积分函数，不对原始模型速率作时间序列插值，也不把参考曲线当积分输入。
原点两侧分别用 longdouble 累计真实分段增量，避免百万 Python checkpoint 对象。

运行全覆盖前，现代 732 个区间已比较上述配置与半日 GL12：主要体主法/独立法对
半日高阶法最大累计差分别约 `6.17e−18 s`、`2.08e−18 s`；SB441 分别约
`4.00e−26 s`、`3.88e−26 s`。这仅是现代步长/阶数检查，全覆盖误差必须看下节
实际两种算法完整累计之差。局部阶数差的绝对和只是经验指标，不是物理误差界。

## 5. 实际运行结果

以下数值由实际运行产生，主 session 另从 gzip 每日表逐行独立重算复核，结果与
`results/tables/full_span_lunar_{full,modern}_diagnostic.json` 一致；以 JSON 原始值为准。
比较量 `raw_time_residual = model − LTE440`，单位秒，未拟合零点或趋势。

### 5.1 全覆盖（TDB JD2287184.5–2688976.5，1549-12-31 至 2650-01-25）

命令 `.venv/bin/python scripts/compare_full_span_lunar.py`，实际耗时约 888.6 s（绘图前）。
401,793 个每日网格点加 1977 原点，共 **401,794** 个时刻。

| 指标 | 数值 |
|---|---:|
| 原始时间 residual 最大绝对值 | **402.922349999 ps** |
| 原始时间 residual RMS | 220.829010871 ps |
| 原始时间 residual 最小/最大 | −328.676638857 / +402.922349999 ps |
| 最大绝对值历元 | TDB 2648-11-11 00:00 |
| 原始速率 residual 最大绝对值 | 4.07121133501e−19 |
| 原始速率 residual RMS | 1.00842509662e−19 |
| 两算法完整累计最大差 | `7.28583859910259e−17 s`（约7.29e−5 ps） |
| 锚点 | `Q_TDB0` 模型 `6.55e−05 s`，参考 `6.549999999999212e−05 s` |

每日样本小于 1 ps 的验收状态为 **False**。

### 5.2 现代 732 日点（1976-12-31 12:00 至 2026-01-01 00:00）

命令 `.venv/bin/python scripts/compare_full_span_lunar.py --modern-only`，约 56.8 s。

| 指标 | 数值 |
|---|---:|
| 原始时间 residual 最大绝对值 | **2.830296761 ps** |
| 原始时间 residual RMS | 2.804314140 ps |
| 原始时间 residual 最小/最大 | +2.792562188 / +2.830296761 ps |
| 最大绝对值历元 | TDB 2026-01-01 00:00 |
| 原始速率 residual 最大绝对值 | 1.64266695920e−20 |
| 两算法完整累计最大差 | `7.399679827213568e−17 s`（约7.40e−5 ps） |

现代残差集中在 +2.8 ps 且为同号，说明存在未被吸收的物理/输入缺项而非随机噪声；
本轮未据此外推或拟合修正。

结论：月心全覆盖的**数值求积一致性**得到支持（两独立算法差约 10⁻¹⁶ s），但原始
残差对固定 LTE440 产品最大约 402.9 ps、现代区间约 2.83 ps，**不能宣称满足 1 ps 一致性目标**。
小求积差只说明积分收敛，不能替代物理模型完整性与输入星历不确定度的评估。

图：`results/figures/full_span_lunar_full_raw_residuals.png`（全覆盖）、
`results/figures/full_span_lunar_modern_raw_residuals.png`（现代 732 日点）。

## 6. 复现与审阅入口

```sh
.venv/bin/python -m pytest -q tests/test_lunar_model.py tests/test_lte_reference.py
.venv/bin/python scripts/compare_full_span_lunar.py --modern-only
.venv/bin/python scripts/compare_full_span_lunar.py
```

新测试覆盖 Vondrák equator/ecliptic poles 对独立 ERFA、Park 轴方向与
CSPICE elementary 旋转顺序/arcsec、J2 极点/赤道、改变计算减常率后结果不变、
原始边界、Moon/Earth 身份与真实 SPICE 状态、诊断通道排除，以及原参考样例与
SPK 记录边界/覆盖端点。源码、输入 SHA、参数、两个方向累计、各方法节点数与耗时
均写入 `results/tables/full_span_lunar_full_diagnostic.json`。

图标题、坐标标签与解释采用中文，时间尺度术语仍用 TCL/TDB；图显示原始速率残差、
原始时间残差与独立求积差，没有拟合零点或趋势。任何小于1 ps的结论只允许针对
实际采样与固定发布参考，不代表绝对物理准确度达到1 ps。
