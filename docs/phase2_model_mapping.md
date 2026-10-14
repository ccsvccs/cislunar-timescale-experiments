# Phase 2：Park 原式与独立实现的逐项映射

审计日期 2026-10-14。范围严格为地心。原始来源是已取得的 [Park et al. 2021 PDF](sources/park_2021_de440_de441.pdf)，Sec. 2.3，第 3–4 页式 (3)–(7)，并读 Sec. 2.2/3.1/3.2/3.5；[下载及 SHA256](sources/park_2021_de440_de441.reference.md)。保留 [Turyshev v2](sources/turyshev_2507.22145v2.pdf) 式 (13)/(21) 和 Appendix A 的 Phase 1 模型不变。两者均为坐标时变换；与 JPL 数值时间星历的一致性不代表绝对物理准确度。

## 1. 时间尺度、符号、边界与单位

统一 $t=\mathrm{TDB}$，$\Delta=t-\mathrm{TT}$，$\mu_i=(GM_i)_{\mathrm{TDB}}$，$r_{ij}=|\boldsymbol x_i-\boldsymbol x_j|$，$\boldsymbol d_i=\boldsymbol x_E-\boldsymbol x_i$，$\boldsymbol v_i=d\boldsymbol x_i/dt$，$\boldsymbol a_i=d^2\boldsymbol x_i/dt^2$。Park 原文用 JD；内部积分统一秒，位置 m，速度 m/s，加速度 m/s²，GM m³/s²；$w_0,w_L$ 为 m²/s²，$\boldsymbol w$ 为 m³/s³，$\mathcal D$ 为 m⁴/s⁴。地心条件是站点位移 $\boldsymbol r_S-\boldsymbol r_E=0$；地球的非零重心位置与速度保留。

$$A=\frac{1-L_G}{1-L_B},\quad C=\frac{L_B-L_G}{1-L_G},\quad t_0=T_0+\mathrm{TDB}_0.$$

$L_G=6.969290134\times10^{-10}$，$L_B=1.550519768\times10^{-8}$，$T_0=\mathrm{JD}\,2443144.5003725$，$\mathrm{TDB}_0=-65.5\,\mu\mathrm s$ 均采用 IAU 2006 B3/既有 `constants.py` 的固定定义。Park 式 (3) 的地心形式严格化为

$$\boxed{\Delta(t)=\mathrm{TDB}_0+A\int_{t_0}^{t}[F(s)-C]ds,\qquad\frac{d\Delta}{dt}=A(F-C).}$$

这个化简与既有式 (21) 边界处理相同，没有 JPL 参考偏置或漂移拟合。数值使用两部分 JD。原始锚点采用定义 JD，并且仅加一次 TDB0；不能将 1977 年锚点换成现代校准历元。论文列出的经验 $L_C$ 与固定 $C$ 的差异仍按 [Phase 1 审计](model_audit.md) 原样保留。

TCB→TDB 同时缩放 $t,\boldsymbol x,\mu$，缩放因子 $k=1-L_B$，故速度和 $\mu/r$ 不变；加速度缩放为 $a_{\mathrm{TDB}}=a_{\mathrm{TCB}}/k$，乘位移的 $\boldsymbol a_i\cdot\boldsymbol d_i$ 不变。DE440 状态与配套 GM 一致为 TDB-compatible，禁止再缩放速度或势。所有轨道是几何 ICRF/J2000，未加光行时。

## 2. 从原式读出的完整候选速率

Park 式 (3) 的 $c^{-4}$ 积分外有负号，括号内的动能和混合势项也有负号；因此

$$F=\frac{v_E^2/2+w_{0E}+w_{LE}}{c^2}+\frac{v_E^4/8+(3/2)v_E^2w_{0E}-w_{0E}^2/2-4\boldsymbol v_E\cdot\boldsymbol w_E-\mathcal D_E}{c^4}.$$

$$w_{0E}=\sum_{i\ne E}\frac{\mu_i}{r_{iE}}\tag{Park 4}$$

$$w_{LE}=-\frac{\mu_\odot J_{2\odot}R_\odot^2}{2r_{\odot E}^3}(3\sin^2\varphi_{E,\odot}-1)\tag{Park 5}$$

$$\boldsymbol w_E=\sum_{i\ne E}\frac{\mu_i\boldsymbol v_i}{r_{iE}}\tag{Park 6}$$

$$\mathcal D_E=\sum_{i\ne E}\frac{\mu_i}{r_{iE}}\left[-2v_i^2+\sum_{j\ne i}\frac{\mu_j}{r_{ij}}+\frac12\left(\frac{\boldsymbol v_i\cdot\boldsymbol d_i}{r_{iE}}\right)^2+\frac12\boldsymbol a_i\cdot\boldsymbol d_i\right].\tag{Park 7}$$

**式 (7) 内层求和包含地球**，外层排除地球；不能从外势排除自势的规则错误推出内层也排除地球。`ExtendedDE440Ephemeris` 单独提供 Earth GM。与 Turyshev 相比，Phase 1 $F_4$ 正好是前三个 $c^{-4}$ 项；新增的是太阳 J2、向量势、完整 PN 标量四项和更完整的显式质量源。

## 3. 逐项代码、来源、量纲、状态与验证策略

下列数量级是先验尺度或原文估计，不是实际验收结果。**所有实测最大速率、锚点累计量、现代变化和残差影响以 Phase 2 输出表为准**。每个开关只关掉该物理项，不重新确定 $C$。

| 物理项 / 原式 | 精确速率贡献（无量纲） | 数据及约定 | Phase 1 / Phase 2 与验证 |
|---|---|---|---|
| 地球重心动能，Turyshev21/Park3 | $v_E^2/(2c^2)$ | `(0,3)+(3,399)`，真实 Earth 非 EMB，m/s | 既有 `rate_c2`，约 $5\times10^{-9}$；`kinetic_c2` 开关；单位/缩放/速度测试 |
| 外部单极势，Turyshev21/Park3–4 | $\sum\mu_i/(r_{iE}c^2)$ | Sun、Moon、Mercury/Venus、Mars–Pluto 系统；配套 DE440 GM | 既有 `external_potential`，太阳约 $10^{-8}$；`monopole_c2`；地球不进入外层，卫星质量不重复 |
| 高阶动能，Turyshev21/Park3 | $v_E^4/(8c^4)$ | 同一速度 | 既有，约 $10^{-17}$；`kinetic_c4`；多项式/量纲测试 |
| 混合势，Turyshev21/Park3 | $3v_E^2w_0/(2c^4)$ | 同一完整单极势 | 既有，约 $10^{-16}$；`mixed_c4`；分解与原式一致 |
| 势平方，Turyshev21/Park3 | $-w_0^2/(2c^4)$ | 同一完整单极势 | 既有，约 $-5\times10^{-17}$；`potential_squared_c4`；负号测试 |
| 太阳扁率，Park3/5 | $w_L/c^2$ | `ASUN=696000 km`，`J2SUN=2.19613915165298248e-7`，原始 header；方向见下文 | 新增 `solar_j2_potential`，约数 $10^{-20}$；`solar_j2_c2`；赤道/极点解析正负号、显式轴测试 |
| 外部向量势，Park3/6 | $-4\boldsymbol v_E\cdot\sum\mu_i\boldsymbol v_i/r_{iE}/c^4$ | 真实重心速度；不能用相对 Earth 速度替代 | 新增 `vector_potential`，先验约 $10^{-19}$；`vector_c4`；多体向量和/负号测试 |
| PN 体速度，Park3/7 | $+2\sum\mu_i v_i^2/(r_{iE}c^4)$ | 同一质量体重心速度 | 新增 `pn_scalar_components['velocity']`，`pn_velocity_c4`；解析两体和正号测试 |
| PN 体间势，Park3/7 | $-\sum_i(\mu_i/r_{iE})\sum_{j\ne i}\mu_j/r_{ij}/c^4$ | 内层 Earth GM；相互距离；不能使用 Earth 外势统一替代 | 新增 `body_potential`，`pn_body_potential_c4`；显式内层 Earth 解析测试 |
| PN 径向速度，Park3/7 | $-\frac12\sum_i(\mu_i/r_{iE})(\boldsymbol v_i\cdot\boldsymbol d_i/r_{iE})^2/c^4$ | 同一方向/速度；平方不改径向符号 | 新增 `radial_velocity`，`pn_radial_velocity_c4`；垂直/平行解析测试 |
| PN 加速度，Park3/7 | $-\frac12\sum_i(\mu_i/r_{iE})\boldsymbol a_i\cdot\boldsymbol d_i/c^4$ | SPK Type2 同一位置 Chebyshev 的解析二阶 TDB 导数，SI | 新增 `acceleration`，`pn_acceleration_c4`；多项式一/二阶导数独立测试；不以 Newtonian toy force代替实际轨道 |
| 343 小行星直接势，Park4/Sec2.2/3.1 | $\sum_a\mu_a/(r_{aE}c^2)$ | DE440 header `MA####`；必须实际轨道 | 原核心缺失；SB441 支持独立诊断，见后文。质量总和 $1.7053212\times10^{11}$ m³/s²，不由 Earth 已受摄动轨道自动包含 |
| 30 KBO 直接势，Park4/Sec2.2/3.1 | $\sum_k\mu_k/(r_{kE}c^2)$ | header `MA8001..8030`；必须查明与 NAIF ID 对应 | 原核心缺失；同上；GM总和 $2.3154029\times10^{12}$ m³/s²。不能用一组假定距离作实测 |
| Kuiper 环直接势，Park4/Sec2.2 | $\sum_{n=0}^{35}(\mu_R/36)/(|\boldsymbol x_E-\boldsymbol x_\odot-44\,\mathrm{AU}\,\boldsymbol n_n|c^2)$ | header `MA8201..8236` 等质量；论文半径44AU、黄道平面 | 新增 `kuiper_ring_potential`，`kuiper_ring_c2`；约 $1.5\times10^{-17}$。公开几何候选不等于获取了原 JPL ring 轨道；说明见下文 |
| 地心位置项，Turyshev21/Park3 | $\boldsymbol v_E\cdot(\boldsymbol r_S-\boldsymbol r_E)/c^2$ 和 $(3w_0+v_E^2/2)\boldsymbol v_E\cdot(\boldsymbol r_S-\boldsymbol r_E)/c^4$ | 站点位移严格0 | 两项严格0；诊断明确0；不实现地面钟 |

Code: [jpl_model.py](../src/tt_tdb/jpl_model.py) 的 `JPLTermSwitches`/`diagnostics`/`correction_rates` 对应上述贡献，[extended_ephemeris.py](../src/tt_tdb/extended_ephemeris.py) 提供同一 SPK 多项式导数。对照工具保留 Phase 1 baseline 速率通道的算术；新 reader 内部独立计算同一 Chebyshev 多项式的状态与导数，完整锚点直接积分检查测定两种算术顺序的舍入差。附加贡献在减去 $C$ 后小量相加，避免把 $10^{-20}$ 量级重新舍入进 $10^{-8}$ 总速率。缓存签名包含模型代码、BSP/GM/header 内容、开关、太阳轴和环参数；物理改变拒绝复用旧缓存。

## 4. 太阳方向、质量变化与未确认的模型细节

Park 式 (5) 的文字确实称 $\varphi$ 为 **heliocentric ecliptic latitude**。但轴对称太阳扁率的球谐纬度通常相对于太阳赤道；Sec. 3.5 给太阳自转轴 RA286.13°、Dec63.87°，Sec. 3.2 说太阳二阶带谐与所有天体相互作用。候选选择这个太阳轴，**明确为物理上有根据的解释，不是原 JPL time code 轴的证明**，可显式传入黄道轴比较。进一步独立读取的 [Folkner et al. 2014, IPN 42-196](https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf) Sec. II.C，第6–7页式(7)明确采用 **相对太阳赤道的纬度**，为当前太阳轴解释提供直接 JPL 原始依据；它仍不是 DE440 原time代码执行方式的证明。原 header 没有太阳轴常数；`J3SUN=J4SUN=0`，因此不添加猜测的太阳高阶系数。地球自己的 multipoles 不进入地心外势；Park 式(5)本身只列太阳J2，不能擅自加入其他行星 multipoles 来“改善残差”。行星轨道已包含的 LT、辐射压力和地球/月球 figure 摄动通过实际状态进入；它们不自动赋予新直接时间势项。

Header 另有 `GMSDOT=-9e-14` 和 `JDEPOC=2440400.5`。Park Sec.2.3 未规定如何把这两个量应用到时间积分的 GM，也没有给此导数的单位/时间零点规则；不能仅凭动力学参数名就默认 $μ_\odot(t)$。候选使用配套 GM 常数，质量变化约定列为待核验问题。这项即使微小也可能在几十年积累超过ps，应与实际DE源代码/更完整注释核对，不能用残差推导系数。

## 5. Kuiper 环公开几何与限制

使用 SPICE 内建 `pxform('ECLIPJ2000','J2000',0)`，保留明确的 ECLIPJ2000 平面假设；这是 [SPICE Frames Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html) 的标准黄道J2000框架，不声称 Park 原 JPL 内部环恰用此定义。36点等间隔，整体相位选0，不根据时间参考选择相位。$μ_R=8.916783186240352\times10^{12}$ m³/s²。

若环保持公开描述的等距圆形，位势 Legendre 展开中整体相位依赖从第36阶开始。对 $q=|\boldsymbol x_E-\boldsymbol x_\odot|/R<1$，两个整体相位的势差保守界为 $4(μ_R/R)q^{36}/(1-q)$。Earth 的 $q\le1.02/44$ 下远低于ps目标。**这个界只约束相同圆形几何的整体相位**，不约束实际 JPL 环质点可能的径向演化、非均匀布局或平面差异。没有获取环动力学状态，所以其 vector/PN/mixed $c^{-4}$ 项仍缺失；不能合成环速度/加速度当作原模型。完整小体修正亦需相应 $w_0^2$、混合和 PN 项，与 leading-only 诊断分开。

## 6. 小体真实数据的获取与一致性

实际读取 JPL [small_bodies目录](https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/)、[asteroids_de441目录](https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/)，本地保存两个HTML证据。目录包含DE430/DE441小体文件，没有公开的DE440小体子目录。NAIF通用asteroids目录仅有2010年的300体集，不能匹配DE440343+30体。没有借助现代 Horizons 拟合轨道冒充 DE440。

已下载 [Farnocchia IOM392R-21-005](sources/sb441_perturbers.pdf)，SHA256 `124fe0fb48f8e71f3f1dba9607cfe5b47a016b95aa5429b61bf7cba301ed4b22`，12页，993,300字节；保持原作者与Caltech版权归属。该memo说373与DE440/441的**体集合**一致，但SB441是使用DE441 planetary states作参考、另行迭代轨道拟合的产品，使用1950–2021年光学、Gaia、雷达和Dawn资料。其收敛差不是轨道不确定度。故不能把“同一集合”误写成“同一DE440轨道”。Table1记录各实际NAIF ID和GM，KBO `MA800x` 是header顺序而非NAIF ID。

SB441-N373 full为14.13GB；实际获取短版 `sb441-n373s.bsp`（目录标936.6MB）以测量真实小体 leading potential，大小、SHA、体映射与实测coverage见 [small-body provenance](../data/phase2_small_body_provenance.md)。诊断使用真实DE440 header GM且逐项匹配IOM质量，但 **SB441 states + DE440 planets 是版本混合诊断**；无法单凭该诊断符合度证明精确恢复JPL原time model。若纳入积分，必须从同一个原始1977锚点累积，不能用现代区间局部调零。

## 7. 验证准则与未解决来源

速率 $R_{\rm rate}=A(F-C)-d\Delta_{\rm JPL}/dt$ 与时间 $R_{\rm time}=\Delta_{\rm model}-\Delta_{\rm JPL}$ 独立报告；时间参考只在对照工具读取，模型无导入。首先固定732个2024-01-01..2026-01-01每日TDB样本。只有实际原始最大时间差严格低于1ps才验收。固定定义常数、锚点，禁止参考零点拟合或趋势校准。诊断贡献累计量采用每项独立积分，保留Romberg/Gauss–Legendre交叉检查。已知Phase1约34.585ns误差不可与假设环23ns直接相等；实测表决定每项改变量。

尚未确认：DE440原始343+30小体积分轨道及ring真实动力学状态；时间积分太阳纬度轴；GM随时间变化规则；行星系统质心与真实resolved satellite势差；公开行星SPK插值与原time integrator状态的差别。缺失源是科学模型/输入误差，不等于积分误差。约1100年目标只能在现代初始门槛实测通过后且数据覆盖可行时单独测量，不预先宣称成功。

## 8. 两个原始太阳纬度定义的独立诊断

两份原始JPL文献给出的文字定义不同：Park2021式(5)第3页写 heliocentric ecliptic latitude；Folkner2014 Sec.II.C式(7)第7页明确相对太阳赤道。后者支持物理上通常采用的太阳对称轴，但不能因此悄悄删除前者的字面解释。`compare_phase2_solar_axis.py`以两种事先指定的轴从同一1977定义锚点分别积分太阳J2贡献：太阳赤道轴由Park Sec.3.5的RA/Dec确定，字面黄道轴由SPICE ECLIPJ2000确定；两个值完整并列，默认模型不改变，也不用时间参考选择轴。

实测732个现代日点中，黄道轴减太阳轴的累计差从0.832016992486ps到0.867435797305ps，最大瞬时速率差约1.14896e−21。这个差与1ps阈值相近，所以约定必须显式公布。太阳轴诊断输出只是在同样的轨道、GM、半径下比较原文给出的两个固定定义，不能把更好的参考符合度单独当作DE440内部代码所用纬度定义的证据。

实际加入真实SB441小体leading势后，默认太阳轴候选在指定732个日点的raw最大残差为2.818949568673ps，RMS2.801277755568ps。独立叠加事先规定的Park字面黄道轴差后，raw最大残差为1.953519689910ps，RMS1.950451572377ps；两者均未达到1ps。黄道轴raw现代峰峰仅0.005356609179ps，近似固定偏差，但不能用常量拟合消去它来宣称通过。太阳轴小体诊断最大rate差为7.951273638166e−22，而小体独立quadrature差约1.65e−24s；剩余ps误差不能归因于该已测量quadrature差。DE440原始小体轨道、环动态状态、all-bodyPN以及GMSDOT时间势规则仍未完整恢复，剩余近常量误差未能唯一归因到其中任何一项。两个轴的实测值都保留，不从参考残差挑选一个作为事后校准。1100年完整时间一致性仍须实际另行积分测量。
