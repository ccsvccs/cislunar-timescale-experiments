# 地心 TT↔TDB 式 (21) 模型审计

审计日期：2026-10-13。此文件先于积分器实现建立；公式实现以固定版本 **arXiv:2507.22145v2 的 PDF** 为准。用户提供的 `paper.pdf` 实际为 v3，首页标记 15 Nov 2025，正文日期 November 18, 2025。已下载并逐式阅读 [官方 v2 PDF](https://arxiv.org/pdf/2507.22145v2)，本地见 `docs/sources/turyshev_2507.22145v2.pdf`，SHA256 为 `dbbd9b74ee55f8bd1383d39d4b9aec505dca4a8fbb31531319ae5e4b110ece01`。v2 首页标记 20 Sep 2025，正文日期 September 23, 2025。在线 HTML 的日期及附录编号与 PDF 不一致，故只作为交叉检查。

已阅读 v2 Sec. II.B.1–3、II.C（式 10–22）以及 Appendix A.1、A.3（式 A31–A44）。对照本地 v3 的同一部分，式 (13)–(22) 的核心代数、常数与问题相同。本项目不覆盖其月球章节。

## 1. 符号、适用位置与原式

论文的 $\boldsymbol r_E=\boldsymbol x-\boldsymbol x_E(t)$ 是**观测点相对地心位置**；本任务 $\boldsymbol r_E=0$，两个 $\boldsymbol v_E\cdot\boldsymbol r_E$ 项严格为零。$\boldsymbol x_E$ 是地球在 SSB/ICRF 中的位置，$\boldsymbol v_E$ 是其质心速度，均不为零。地球自身势不进入外部势，月球势进入。

令所有时间间隔以秒计，$t=\mathrm{TDB}$，$t_0=T_0+\mathrm{TDB}_0$，并定义

$$U=\sum_{B\ne E}\frac{GM_B}{|\boldsymbol x_E-\boldsymbol x_B|},\qquad F_2=\frac{v_E^2/2+U}{c^2},\qquad F_4=\frac{v_E^4/8+3v_E^2U/2-U^2/2}{c^4},\qquad F=F_2+F_4.$$

直接视觉检查 v2 PDF **第 7 页**：式 (21) 是一个线性项**减去**整个大括号（其中两个积分均为正号）。这避免 PDF 文本提取的换行误读。地心原式为

$$\mathrm{TT}-t=\frac{L_B-L_G}{1-L_B}(t-T_0)-\frac{1-L_G}{1-L_B}\left[\mathrm{TDB}_0+\int_{t_0}^{t}F(s)\,ds\right].$$

此处 $T_0$ 是 JD 数值定义的历元，不能直接与秒相加：实现应保留 $2443144+0.5003725$ 的日期分拆，再单独加 $-65.5\,\mu\mathrm s$。$T_0$ 的日内数值 32.184 s 是 TT/TCG/TCB 读数，对应事件是 **1977-01-01 00:00:00 TAI**，不是论文文字中误写的 00:00:32.184 TAI；[IAU 2006 B3 的注 1](https://www.iau.org/static/resolutions/IAU2006_Resol3.pdf) 与 [BIPM CCTF 2017 Recommendation 3](https://www.bipm.org/en/-/cctf-recommendation-2017-3) 支持此解释。代码只采用定义 JD，不根据论文的 TAI 文字另加 32.184 s。

## 2. 稳定积分形式的严格代数证明

取 $A=(1-L_G)/(1-L_B)$，$C=(L_B-L_G)/(1-L_G)$，故 $AC=(L_B-L_G)/(1-L_B)$ 且 $A(1-C)=1$。反转符号得到

$$\Delta(t)\equiv t-\mathrm{TT}=A\left[\mathrm{TDB}_0+\int_{t_0}^{t}F(s)\,ds-C(t-T_0)\right].$$

因为 $t-T_0=(t-t_0)+\mathrm{TDB}_0$，所以

$$\boxed{\Delta(t)=\mathrm{TDB}_0+A\int_{t_0}^{t}[F(s)-C]\,ds.}$$

于是 $\Delta(t_0)=\mathrm{TDB}_0=-65.5\,\mu\mathrm s$，$d\Delta/dt=A(F-C)$。此式与式 (21) **精确代数等价**，未拟合常数或强制去趋势。应将它作为主计算链；原式只作数值交叉检查。所有 $c^{-4}$ 贡献都保留，不在构造 $F-C$ 时重新调平均常数。

## 3. TDB-compatible 状态与 GM

式 (13) 明示 TCB-compatible；式 (20) 同时缩放时间、空间与 GM。令 $k=1-L_B$，则

$$t_{\rm TDB}=kt_{\rm TCB}+\mathrm{const},\quad \boldsymbol x_{\rm TDB}=k\boldsymbol x_{\rm TCB},\quad (GM)_{\rm TDB}=k(GM)_{\rm TCB}.$$

因此 $\boldsymbol v_{\rm TDB}=d\boldsymbol x_{\rm TDB}/dt_{\rm TDB}=\boldsymbol v_{\rm TCB}$，$GM/r$、$v^2$ 与本次 $F_2,F_4$ 在一致缩放下不变，而 $dt_{\rm TCB}=dt_{\rm TDB}/k$，这正是式 (21) 外面的 $A$ 的来源。DE440 的 TDB-compatible 坐标和其配套 TDB-compatible GM 可直接使用，不应再缩放速度或额外乘一次 $k$。以 [JPL DE440/DE441 说明](https://ssd.jpl.nasa.gov/doc/de440_de441.html)、实际星历 header/GM 核验以及 [SPICE Time Subsystem](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html) 为数据来源。ET 是 J2000 TDB 起算秒；使用几何状态 `J2000`、`NONE`。

地心状态优先使用 SPK 路径 SSB→EMB→Earth：EMB(3) 与 Earth(399) 不同，月球(301) 独立。外部行星应把系统质心位置与系统 GM 配对，避免再加已含在系统 GM 中的卫星。SPICE km、km/s、km³/s² 必须在接口处乘 $10^3,10^3,10^9$ 转 SI。仅主行星/月球 SPK 不等于完整 DE440 势模型：小行星、30 个 KBO 与 KBO 环可能只通过轨道扰动进入状态，而没有进入显式势求和；须另行记录并估计遗漏，不能宣称完整物理精度。

## 4. 定义常数及 $L_C$ 一致性

| 常数 | 采用值 | 性质与来源 |
|---|---:|---|
| $c$ | 299792458 m/s | SI 定义；[BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure) |
| $L_G$ | $6.969290134\times10^{-10}$ | IAU 2000 B1.9 定义；[IAU NSFA fixed constants](https://iau-a3.gitlab.io/NSFA/NSFA_fixed) |
| $L_B$ | $1.550519768\times10^{-8}$ | IAU 2006 B3 定义，无拟合误差条 |
| $T_0$ | JD 2443144.5003725 | 同一 IAU B3 定义 |
| $\mathrm{TDB}_0$ | $-65.5\times10^{-6}$ s | 同一 IAU B3 定义 |
| $L_{C,\mathrm{paper}}$ | $1.48082685455\times10^{-8}$ | v2 式 (15) 后的经验长期平均，不能当作定义常数 |

用十进制 60 位算术独立计算：

$$C_{\rm exact}=1.48082686769203120791681964477547938429911331181162049089568\times10^{-8},$$
$$A=1.00000001480826889620513653425604379503005077667833883283518,$$
$$C_{\rm exact}-L_{C,\rm paper}=1.314203120791681964477547938429911331181162049089568\times10^{-16}.$$

这个差异在 $\Delta$ 中若被漏掉，幅度为 $11.3547151318$ ps/day，或 $4.1473097019$ ns/Julian year。直接把论文 $L_C$ 代入式 (19) 得到 $L_B=1.550519754857968801242243199030\times10^{-8}$，与定义值不符。**不能以式 (19) 来修改 $L_B$ 或 $L_C$。**

IAU B3 注 2 将固定 $L_B$ 描述为根据 DE405 经验 $L_C$ 的估计，并允许换用星历后存在小漂移；它没有定义所有将来星历的 $\langle d\mathrm{TDB}/d\mathrm{TT}\rangle$ 必须精确为 1。因此论文式 (18) 的逐时 $d\mathrm{TDB}/d\mathrm{TT}=1$ 与其周期项也不相容，应理解为理想平均关系。

本审计确认了数值矛盾，却**没有足够原始证据把特定差值归因于某个星历或勘误**。论文未给所列 $L_C$ 的计算星历、平均区间或足够推导。按文献中完整与截断模型的区别，遗漏的完整 $c^{-4}$ 项可能影响平均速率；本项目也不能把两者相近的数量级当作归因证明。

若需要用 $L_{C,\rm paper}$ 分离积分，应严格恢复

$$\Delta=\mathrm{TDB}_0+A\left\{\int_{t_0}^{t}(F-L_{C,\rm paper})\,ds+(L_{C,\rm paper}-C_{\rm exact})(t-t_0)\right\}.$$

## 5. 式 (22) 不是可直接替代式 (21) 的精确公式

v2 PDF 第 7 页视觉核对表明它确实写的是 $P(t)-P(t_0)$；文本提取表面上的连续两个负号来自大括号丢失，不能改成 $P(t)+P(t_0)$。

但从式 (14)–(16) 若定义 $\dot P=F-L_C$，地心式 (21) 严格给出

$$\Delta=\mathrm{TDB}_0+A[P(t)-P(t_0)]+A(L_C-C_{\rm exact})(t-t_0).$$

论文式 (22) 省去了 $A$ 与常数不一致造成的线性项。即使设 $L_C=C_{\rm exact}$，省去 $A$ 也在 1.66 ms 年周期上产生约 **24.6 ps** 振幅差；因此在亚 ps 分析中不可忽略。这是从原式推导得到的审计结论，而非已发布勘误。主积分严格执行 (21)，报告保留这些问题；不为了满足“无长期漂移”文字强制处理结果。

## 6. $c^{-4}$ 特定模型及余项含义

本次 $F_4$ 的三个项直接对应式 (21)：$+v^4/8$、$+3v^2U/2$、$-U^2/2$，不是完整 N-body 1PN 时间星历。Appendix A31–A42 中完整标量势还包含 PN 修正，向量势 $\boldsymbol w_{\rm ext}$ 进入 $4\boldsymbol v_E\cdot\boldsymbol w_{\rm ext}$；本次按原式截断，禁止把 LTE440/JPL 的完整速率无说明地换进来。

余项 $2.14\times10^{-19}(t-t_0)$ 来自 Appendix A41 对遗漏外部向量势**速率**项的数量级估计。它每日累积 $0.0184896$ ps，1 Julian year 约 6.753264 ps，50 年约 0.3376632 ns，不是整套转换的固定绝对误差。$1.91\times10^{-16}$ s 来自 A42 GPS 轨道位置相关项的估计，是 0.000191 ps；地心 $\boldsymbol r_E=0$ 时该项严格消失，不应当作地心全部误差。

同时，附录 A38 估计太阳 $J_2$ 的速率为 $4.80\times10^{-20}$，A40 还讨论被省略的混合势修正。有限天体集合、DE440 状态/GM 误差、数值积分与时间表示误差均不在前述两个数字里，必须独立预算。论文正文式 (1) 后举例中把 $v^4/8$ 写成 $v^2/8$，该例量纲错误；Appendix A41 的叙述估值 $2.07\times10^{-16}$ 与它列出的加减项也不一致。实现遵循有正确量纲的原始公式，而不复现这些粗略估值。正文给出的 $F_4\sim1.10\times10^{-16}$ 可作数量级诊断，不当作硬阈值或真值。

## 7. LTE440 算法借鉴与验证边界

已阅读 [Lu, Yang & Xie Section 3](https://arxiv.org/html/2509.18511#S3)：采用 10 阶 Romberg、0.5 day 分段，对速率减常数并最终恢复线性项；比较 0.5/1 day 的约 10 fs 是该论文的实验，不能写成我们的结果。其地球独立重建与内置 DE440 时间星历小于 1 ps 也依赖其完整模型与数据。

本地 LTE440 HEAD `76ea3a97fdbbd645352d3855b91dd2e878221351`，工作树仅 README；没有原始积分程序，也没有 LICENSE 文件。README 提供的是数据读取样例与 Releases 链接。网页 Release `20250620` 明确提供 `lte440-tclmtdb.zip` 与中文 PDF；它不是本任务的地球时间星历，不下载月球大数据。项目独立写积分代码，只借鉴算法；无明确许可时不把 README 样例复制成工程实现。来源详见 `source_notes.md`。

ERFA `dtdb.c` 只用于地心 Fairhead–Bretagnon 对照；输入 `u=v=0` 令地面项为零，输出符号是 TDB−TT 秒。其说明的 DE405 约 ±3 ns（1950–2050）不构成 DE440 精确真值。`tttdb.c` 只是将外部提供的时间差加到双部分 JD 的较小分量，并没有独立的核心时间模型。[ERFA dtdb 原始源码](https://github.com/liberfa/erfa/blob/master/src/dtdb.c)、[tttdb 源码](https://github.com/liberfa/erfa/blob/master/src/tttdb.c)。

JPL 参考数据必须核查原始 header 的 TT−TDB 分量和实际 SPK segment；标准 `de440.bsp` 不应被假定含该分量。若真实数据能得到，则符号、边界条件、完整模型差异必须分别报告；SOFA 不能替代这一缺失数据。

后续实际数据审计已取得官方 `de440t.bsp`：普通 `de440.bsp` 的14段没有时间目标，普通ASCII header末尾时间指针系数/子段数为0；`de440t.bsp` 则真实含 Type2段 `(1000000000,1000000001)`。按官方README，X是TT−TDB秒，Vx是其TDB导数；`jplephem`的导数按day给出，故分别取−X、−X′/86400。SPICE独立读取已核验。锚点JPL差值−65.50000000001775μs，与定义边界一致到约17.75as。调查与实际数据哈希见 [data/README.md](../data/README.md)，最终时间差/导数结果见 [教学报告](tt_tdb_eq21_report.md)。

## 8. 实现准入结论

理论审计已完成，准许实现上述残余速率积分。接口必须写明 TDB-compatible、SI、地心及有限外部天体集合；保存独立的 `L_C_paper` 与 `C_exact`。本任务复现的是式 (21) 固定边界条件的数值模型，和标准工程参考的常数偏差或漂移必须原样展示。结果不得用任意零点校准、线性去趋势或 SOFA 生成的检查点替代从 $t_0$ 的完整积分。
