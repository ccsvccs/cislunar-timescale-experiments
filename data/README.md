# 真实 JPL 数据与读取约定

`python scripts/download_data.py` 下载完整 DE440（119,799,808 字节）、配套
GM 参数、真正 DE440 数值时间参考和审计元数据；`--core-only` 只下载核心。
`--verify-only` 离线核验已下载文件。下载先写 `.part`，大小和 SHA256 成功
后才改名。大文件不纳入 Git。可设置 `DE440_BSP`、`DE440_GM` 和
`JPL_TIME_BSP` 指向另一个位置。

文件 URL、字节数、SHA256 及实际对象清单见 [manifest.json](manifest.json)。
完整 DE440 SHA256 为
`a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5`；
其 MD5 `c9d581bfd84209dbeee8b1583939b148` 与
[NAIF 官方清单](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aa_checksums.txt)
一致。其他 SHA256 是本次由官方 HTTPS 取回后记录的复现指纹，并非独立签名。
实际共同覆盖为 TDB JD 2287184.5–2688976.5。

## 轨道状态

`de440.bsp` 实际为 14 个 Type 2 段。SSB→EMB 段为 `(0,3)`，
EMB→Earth 段为 `(3,399)`，两者相加得到地球质心位置和速度；
SSB→Moon 则为 `(0,3)+(3,301)`。不能将 EMB 当作地心。
jplephem 返回 km 与 km/day，本项目适配层统一转换成 m 和 m/s。
使用几何 J2000/ICRF 状态，不加光行时或恒星光行差。
输入是 TDB 秒（自 J2000）或双部分 TDB JD。

所有状态与 GM 为 **TDB-compatible** 数值；适配层不擅自缩放。
物理模型必须按常数审计说明其尺度转换。
`gm_de440.tpc` 原单位 km³/s²，转换乘 $10^9$。只读取 DE440 来源的
Sun、Moon、Mercury/Venus 和 Mars–Pluto 系统 GM；Earth 不进入外势，
EMB 不进入外势，也不再叠加已包含在行星系统 GM 内的卫星质量。
PCK 中额外卫星/小天体 GM 来源并不全部是 DE440，不能直接混入本模型。

## TT−TDB 数值参考：已取得真实 DE440 分量

普通 `de440.bsp` **没有时间分量**。普通原始 ASCII 的 `header.440`
GROUP 1050 最后两组系数数与子区间数均为 0，也不提供该分量。
但已逐项核实 [JPL 官方 BSP 目录](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/)
实际列有 **de440t.bsp**（152,747,008 字节）。实测它含前述轨道段和一个
`(1000000000,1000000001)` Type 2 时间段。
[官方 README](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/README.txt)
规定此特殊段 X 是 **TT−TDB，单位秒**，VX 是其相对 TDB 的一阶导数。
它虽然套用 SPK 容器，X 并不是 km，不能乘 1000。jplephem 的导数
按 day 给出，必须除 86400；模型符号 $Δ=\mathrm{TDB}-\mathrm{TT}$
因此取 **−X**。

`JPLTimeEphemeris.evaluate(t)` 返回 $(Δ,dΔ/d\mathrm{TDB})$。
独立 SPICE 原始读取验证了符号和单位。J2000 实测
$Δ=-9.9290115710493\times10^{-5}\,\mathrm{s}$；积分定义锚点
$t_0=T_0+\mathrm{TDB}_0$ 实测
$Δ=-6.550000000001775\times10^{-5}\,\mathrm{s}$，无需常数校准。
对照参考和核心积分保持独立。

## 明确缺失的外部引力源

[原始 header.440](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440)
列有 343 个小行星、30 个单独 KBO 以及 36 个 Kuiper 环质点的 GM；
这些天体轨道没有出现在本次 14 段行星 BSP 中。DE440 已将它们对行星
运动的摄动计入轨道，不能因此认为它们的**直接外势**也自动计入本计算。
[small_body_gm.json](small_body_gm.json) 从该头文件逐项解析，原单位
AU³/day²，取 AU=149597870700 m；记录所有质量项与分组总和：

| 缺失项 | 总 GM (m³/s²) |
|---|---:|
| 343 小行星 | $1.7053212133851767\times10^{11}$ |
| 30 单独 KBO | $2.315402864753561\times10^{12}$ |
| 36 环质点 | $8.916783186240352\times10^{12}$ |

Ceres、Pallas、Vesta 总 GM 分别为 $6.262888864440994\times10^{10}$、
$1.3665878145967424\times10^{10}$、$1.7288232879171516\times10^{10}$ m³/s²。
如仅作数量级估计，三者合计在地球距离 1 AU 时约为 0.6256 m²/s²，
即 $6.96\times10^{-18}$ 的速率、约 0.22 ns/year。
这是假设距离的说明性估计，不是缺失项的严格上界或实测误差。
报告需结合其实际缺失与 JPL 模型差异解释数值结果。

环项可建立条件明确的几何界：
[Park et al. (2021), Sec. 2.2](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf)
给定环的日心半径 44 AU。当地球日心距离不超过 1.02 AU 时，各环质点
到地球的距离处于 42.98–45.02 AU，于是环的直接牛顿势处于
1.3240–1.3868 m²/s²。仅该项的 $U/c^2$ 对时间积累约为
0.4649–0.4869 ns/year。这是给定距离条件下的上下界，尚不能替代
343 小行星和 30 KBO 的实际逐项计算，也不能直接当作与 JPL 的总差异。
