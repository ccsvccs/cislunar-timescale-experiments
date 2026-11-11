# 地心 TT–TDB 全覆盖独立积分报告

本轮模型采用 Park 主要天体与公开 ring leading 势、实际 SB441 小体 leading 势、
月球静态 J2+C22；不加入太阳直接 spin。太阳 J2 默认太阳赤道轴，黄道轴只作为
预先声明的副诊断，未依据残差选择轴。时间差统一为 Δ=TDB−TT，残差为 model−JPL。

## 覆盖与原始结果

真实共同覆盖为 JD 2287184.5–2688976.5，
实际对应 **TDB 1549-12-31 00:00:00 至 2650-01-25 00:00:00**，
即名义 1550–2650 年范围加星历边缘缓冲。包含两端共 401793 个每日 TDB 样本。所有时间值均由真实模型
从原始 1977 锚点完整累计，未减常数或拟合趋势。

| 指标 | 实测 |
|---|---:|
| 全覆盖原始时间 residual 最大绝对值 | 388.373581992 ps |
| 全覆盖原始时间 residual RMS | 215.883742323 ps |
| 全覆盖原始时间 residual 最小/最大 | -388.373581992 / 318.859723944 ps |
| 全覆盖原始速率 residual 最大绝对值 | 4.25114563952e-20 |
| 全覆盖原始速率 residual RMS | 2.54187979396e-20 |
| 现代 732 日点原始时间 residual 最大绝对值 | 2.0735609139 ps |
| 独立完整累计最大差 | 4.59701721134e-05 ps |

每日样本小于 1ps 的验收状态为 **False**。
这只是对指定 JPL 数值时间星历的采样一致性，不代表连续时间严格上界或绝对物理精度。
时间 residual 最大绝对值出现在 TDB 2650-01-25；速率 residual 最大绝对值
出现在 TDB 2157-10-08。两者是从完整每日 payload 独立复算得到的实际采样历元。

## 数值方法与独立验证

快通道（Phase1、Park 扩展、月球静态势与轴差）使用 GL12/day；
373 体 leading 势使用 GL4/day。独立方法为快通道 Romberg6/day、
小体 Romberg5/4day。独立方法也从同一物理锚点完整积到两个覆盖端点，
在全部 100449 个共同四日节点比较累计时间，
不是用稀疏 rate 比较代替长期 time 验证。

物理锚点保持 T0+TDB0；因日网格未恰好经过该时刻，先直接求积
anchor→邻近日网格点，再分别向过去和未来作 Neumaier 补偿累计。
该严格积分拆分不重新定义锚点、不读取参考值作为边界。
每个日/四日区间均实际计算模型；缓存按源码、全部数据 SHA、方法、网格与锚点指纹隔离。
三处分离的 32 日窗口先与高阶 GL24/GL16 比较；实测与全覆盖比较见 JSON。
求积差异是经验数值诊断，不能充当物理模型完整性证明。

现代结果另与已有 SB441 + 静态 lunar J2/C22 表逐点复现，最大差
1.23599047663e-05 ps。
所有旧表 epochs 均逐点核对；旧表哈希与当前完整模型指纹均记录在 JSON。
另在 `full_span_earth_archived_input_identity.json` 核对 Phase1、Park、SB441、
月球静态势与先前报告的完整物理 signature。主 session 将旧图文字改为中文后重跑
其 producer，旧报告 JSON 的 producer SHA/耗时发生变化；该身份记录同时保留原档案
与当前档案 SHA。排除 lunar 的 producer 实验脚本字段后，物理 signature 全部一致，
当前旧 modern CSV 的 SHA 也与本轮实际读取值一致；未为排版变化重跑全覆盖积分。

## 来源与限制

物理映射沿用 [Phase2 映射](phase2_model_mapping.md) 与
[月球静态势来源](../data/phase2_lunar_multipole_provenance.md)。旧物理源码未改。
SB441 是以 DE441 为参考单独拟合的真实小体轨道，不能称为原 DE440 小体状态。
公开 ring 几何、原 ring 动力学状态及小体/ring 完整 PN 项的限制继续保留；
月球仅保留静态 J2/C22，不含时变 tide/rotation distortion。
GMSDOT 时间势约定也未由本轮解决。JPL 时间数据仅在独立对照阶段读取。

副诊断黄道轴全覆盖最大时间 residual 为
376.440515901 ps；
该数值不用于替换预先指定的默认轴或宣称验收。

## 复现与产物

```sh
.venv/bin/python scripts/compare_full_span_earth.py
```

实际总耗时 765.13s；依赖版本、真实命令、输入/源码/产物 SHA
和方法统计保存在 `results/tables/full_span_earth_summary.json`。
原始每日序列为 `full_span_earth_daily.csv.gz`，无损 gzip 与 17 位有效数字；
独立累计序列为 `full_span_earth_independent_cumulative.csv.gz`。
图 `results/figures/full_span_earth_residuals.png` 使用全部每日数据，
上为 rate residual，下为 time residual，均保持原始零点及趋势。
图标题中的 JD 采用近似科学计数显示；精确覆盖端点以本报告首节及 JSON 为准。
