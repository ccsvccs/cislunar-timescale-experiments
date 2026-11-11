# TT 与 TCL 全覆盖计算的独立科学审查

审查日期：2026-11-11。范围为地心 `Δ=TDB−TT` 对 DE440t，以及月心
`Q=TCL−TDB` 对独立 LTE440 发布产品。用户已授权扩展至月心 TCL；本轮不涉及
TL 的约定或具体月面钟固有时。审查与实现分工独立，未修改两实现 agent 的源码。

本记录核查原式、实现、来源与结果；地心完整累计结果已复核，月心最终结果已于
结果完成后补入第 5 节，并由审查侧独立重算核对。
“每日全覆盖采样”不是连续时间最大误差的数学证明。

## 1. 月心 TCL 公式、边界与参考

令 `s=TDB`、`k=1−LB`、`s0=T0+TDB0`。月心没有位置几何项，定义
`FM=F2+F4`，其中 F4 包含矢量势与补充 PN 标量势的正确符号，则

\[
Q(s)=-TDB_0+\frac{L_B}{k}(s-s_0)-\frac1k\int_{s_0}^{s}F_M(u)\,du.
\]

因此 `Q(s0)=−TDB0=+65.5 μs`，`Q'=(LB−FM)/k`。使用预先固定、与参考无关的
计算中心常数 `F*=1.4825e−8` 后，应恢复

\[
Q(s)=-TDB_0+\frac{(L_B-F_*)(s-s_0)-\int_{s_0}^{s}[F_M(u)-F_*]\,du}{k}.
\]

已视觉核查固定 LTE440 v1 PDF p3 式(6)–(7)、(11)：式(7)外层大方括号覆盖
整个 c−2 与 c−4 积分，二者均带 `1/(1−LB)`。文本抽取的换行容易掩盖括号范围；
本审查以 PDF 图形为准。

`lunar_model.py` 的 MoonObserverEphemeris 将观察中心改为真实月心，外势去掉 Moon、
加入 `earth_external`；PN 内层的观察体质量改为 Moon GM。外层禁止 Moon 自势，
不代表 PN 源间势内层可以删去 Moon。向原 Park 代数兼容层传入旧字段名的含义
已在源码显式记录。`rate` 与 `restore` 的正负号、缩放与计算中心常数恢复符合上式。

`lte_reference.py` 独立读取指定发布产品，完整函数为
`BSP_X(s)+RATE×(s−s0)`；导数为 `BSP_X_derivative_per_day/86400+RATE`。
`RATE=+6.798355238724282e−10`，不能省略或改为论文另行讨论的 LTE441 校正平均率。
BSP 的 X 单位是秒。产品哈希、版本与锚点详见 [月球验证方案](lunar_validation_plan.md)。

独立轻量核验使用新读取器与官方 CSPICE 对照，检查覆盖两端、1600/2000/2200 的
五个历元：完整 Q 值差与完整导数差在该计算中均为零；1977 锚点的
`Q+TDB0=2.710505431213761e−20 s`。这是读取与单位/边界验证，不能替代模型验证。

## 2. 地心模型与完整累计

已只读核查 `full_earth_model.py` 与 `compare_full_span_earth.py`：主组合为
Park 主要体、公开环 leading 势、实际 SB441 leading 势，以及月球静态 J2+C22；
未导入或加入太阳直接自转项。太阳 J2 采用预先声明的太阳赤道轴；黄道轴差单独保存。
12 个 fast 通道和2个小体通道的组合索引正确，主结果不混入副轴差通道。

地心脚本主方法对每日全部真实区间求积；独立 fast 方法同样逐日累计，独立小体
方法对连续4日区间全部求积。原始锚点至相邻数值网格点的短积分单独计算，再从该
点分别向覆盖两端作 Neumaier 累计；这与原始物理锚点积分等价，未使用参考零点。
独立小体结果在约100449个共同节点比较的是**完整累计积分**，不是稀疏速率采样。

月心脚本将原始锚点直接加入每日网格，各区间均以真实 integrand 求积，并以
long double 从锚点分别向前/向后累计。计算中心常数只作数值处理且完整恢复。
两脚本可预读参考身份与覆盖来固定比较域，参考值只在物理积分完成后计算。
已审查的物理 integrand 均未读取时间参考，也未拟合零点、趋势、轴或质量。

## 3. Earth J2 极轴与新归档来源

Park 2021 Sec.2.5 p4–5、式(20)–(26)给出 Vondrak 长期岁差与18.6年章动。
式(26)写成行基底时，其行是 date 坐标基轴在惯性系中的分量，故 date 极向量回到
惯性系需要该基底的转置。代码按此方向构造；测试另以 ERFA 的 equator/ecliptic
pole 函数验证独立系数实现，ERFA 未进入 integrand。静态 Earth J2 与拟合极偏置
须分别声明；不得从字段名称反推未确认的时间律或单位。

新增一级来源已下载、验证 PDF 身份，并保留原版权和署名：

| 来源 | 原始链接与本地存档 | 本次核查范围与身份 |
|---|---|---|
| Folkner, Williams, Boggs, Park & Kuchynka (2014), *The Planetary and Lunar Ephemerides DE430 and DE431*, IPN Progress Report 42-196 | [JPL 官方 PDF](https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf)；[存档](sources/folkner_2014_de430_de431.pdf) | 81页，3651368字节；实际核查 Sec.II F p10–11、p50 Table10，并视觉核 p50。SHA256 `76df76941ccc2602df03daa5bd1a3a0d5a54ba20e313c78f42f572fac3ed196f`。不是全文审读声明。 |
| JPL DE430 ASCII header `header.430_572` | [JPL 官方文件](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de430/header.430_572)；[存档](sources/header.430_572) | 20343字节；核查 GROUP1040/1041 配对字段及4个极偏置参数。SHA256 `ed4c6f777ee9d17b8f6f498e8a55538a12c9efee3acd21a0fe3cfc08f64f16a1`。 |

Table10 中 `Φx0=+5.6754203322893470e−3`、`Φy0=−1.7022656914989530e−2`
的单位明确为 arcsec；它们与该 DE430 header 的 ROTEX/ROTEY 原值逐字一致。
DROTEX/DROTEY 也与表中 arcsec/yr 的两个率值对应，p11 给出偏置旋转顺序。
这直接闭合 **DE430** 字段身份与单位。将相同字段解释延用到 DE440 是结合 Park
相同符号的跨版本推断，须如此标记，不能称为已找到 DE440 参数单位表。

本轮在正式全覆盖计算前预先固定 DE440 ROTEX/ROTEY 按上述 arcsec 解释的候选，
零偏置轴作为独立敏感性。已再次核查被动 `Ry(−Φy)` 后 `Rx(−Φx)` 的符号与顺序。
第11个通道仅存零偏置减主轴的差，实际 Q 与速率只恢复前10个物理通道；新增测试
以 ERFA 独立 date 基底和 CSPICE 基本旋转核查该顺序及诊断通道排除。
这不是按对 LTE440 的 residual 选择。J2EDOT 的定义历元和时间势使用范围仍不能
从上述映射推出。现阶段静态 J2 的范围限制必须保留。

附带来源冲突：同一 Folkner PDF p50 Table9 将太阳极轴 RA 印为268.13度，
而 Park2021 Sec.3.4 p7 为286.13度。本次只记录冲突，不据旧表替换既定 Park 极轴。

Vondrák、Capitaine 与 Wallace 岁差模型的主引用为
[2011原文 DOI](https://doi.org/10.1051/0004-6361/201117274) 与
[2012 Corrigendum DOI](https://doi.org/10.1051/0004-6361/201117274e)。本次直连
A&A 的 www/非 www 原文和更正 PDF 均返回 HTTP403，未取得官方 PDF 字节，
不为其编造哈希。实际归档的是[第三方托管合辑](https://www.archaeoastronomy.it/Archaeoastronomy%20Program/Vondrak/Vondrak.pdf)
的[镜像字节](sources/vondrak_2011_2012_annotated_mirror.pdf)：20页，2335737字节，
SHA256 `45429bccb663a19f9ee15732fdac3c5dcfa161318869eabd711475f2fbb7c5df`。
该合辑包含19页原论文、1页正式更正及镜像添加的意大利语勘误提示，不能称为未改动的
官方原版 PDF。原署名、DOI 与 ESO 版权声明保留；PDF 解析器报表单字段警告，
正文及20页均可读取。

实际核查 p4–5 式(8)/(9)及 Tables1/2 的两组多项式和周期系数，并视觉核查第20页
正式更正：`QA` 的 `C7` 应为 `198.296701`，不是 `198.296071`。实现采用更正值；
两组系数与核查范围内原文一致。未把主文搜索索引的个别转录误差当作系数来源，
也未声称全文物理模型均经审查。

## 4. 地心完整结果的独立复核

读取实际 gzip CSV 全行重新计算统计，不复用 summary 内的指标：

| 指标 | 独立复核值 |
|---|---:|
| 地心每日样本数 | 401793 |
| 原始时间 residual 最大绝对值 | 388.37358199220184 ps |
| 原始时间 residual RMS | 215.88374232348207 ps |
| 独立完整累计比较点数 | 100449 |
| 两方法完整累计最大差 | `4.597017211338539e−17 s`（约0.00004597 ps） |
| 现代732日点重组与前轮组件结果最大差 | `1.235990476633475e−17 s` |

每日表两端为 TDB J2000 秒 `−14200747200` 与 `20514081600`，所有步长均为
86400秒；独立累计表包含相同两端，步长为345600秒。每日表 `residual=model−JPL`
逐行闭合差为零。summary 保存的全部输出 SHA256 以及当前 producer、专项测试和
积分源码 SHA256 均匹配实际文件。已视觉核查中文图，标签与单位可读、上下两图的
速率/时间量及图例相符；图题使用约数科学计数，精确边界以表与报告的 JD 为准。

结论：地心全覆盖的**数值求积一致性**得到支持，但原始样本相对 DE440t 的
最大 residual 明显大于1 ps，不能宣称通过该目标，不能将小求积差解释为物理准确度。
证据：[地心报告](full_span_earth_report.md)、
[summary](../results/tables/full_span_earth_summary.json)、
[每日表](../results/tables/full_span_earth_daily.csv.gz)、
[独立累计表](../results/tables/full_span_earth_independent_cumulative.csv.gz)。

## 5. 月心最终结果复核

两套实际累计均已完成。审查侧从实际 gzip 每日表逐行重算 `raw_time_residual_s`
的统计，不采信 diagnostic JSON 的汇总值；重算与 JSON 逐项相符。

| 指标 | 全覆盖（401,794 时刻） | 现代 732 日点 |
|---|---:|---:|
| 原始时间 residual 最大绝对值 | 402.922349999 ps | 2.830296761 ps |
| 原始时间 residual RMS | 220.829010871 ps | 2.804314140 ps |
| 原始时间 residual 最小/最大 | −328.676638857 / +402.922349999 ps | +2.792562188 / +2.830296761 ps |
| 两算法完整累计最大差 | `7.28583859910259e−17 s` | `7.399679827213568e−17 s` |
| 锚点 `Q_TDB0` | 模型 `6.55e−05 s`，参考 `6.549999999999212e−05 s` | 同左 |

全覆盖每日点 401,793 个加精确 1977 原点，共 401,794 个时刻，两端为 TDB J2000 秒
`−14200747200`（1549-12-31）与 `20514081600`（2650-01-25）。最大原始残差出现在
TDB 2648-11-11；现代区间残差为同号约 +2.8 ps。

结论：月心全跨度的**数值求积一致性**得到支持，但原始样本相对固定 LTE440 产品的
最大残差远大于 1 ps，**不满足 1 ps 一致性目标**；现代 732 日点亦为 2.83 ps。小求积差
不能被解释为物理准确度。月心 full-span 的图像已产生
（`full_span_lunar_full_raw_residuals.png`、`full_span_lunar_modern_raw_residuals.png`），
审查侧未逐像素复核，以报告文字与表为准。
证据：[月心报告](full_span_lunar_report.md)、
[full diagnostic](../results/tables/full_span_lunar_full_diagnostic.json)、
[modern diagnostic](../results/tables/full_span_lunar_modern_diagnostic.json)。

后续验收前仍需说明的限制不变：

- 完整两套实际累计已结束并报告两端、全日点、1977锚点与现代732日点原始残差。
- 独立累计差已单列；两套求积相符不证明物理模型完整。
- 实际输入、源码与参考的哈希记录在 diagnostic JSON；结果对应冻结实现。
- <1 ps 阈值不能用于当前结果；SB441 版本混合、原始 DE440 小体轨道/环未获得、
  月球或地球时变重力遗漏等均为物理/输入限制。
