# Research 工作约定与新 session 交接入口

本文件适用于仓库内的主 session 和全部 subagent。用户明确授权主 session
调用、创建、指派和调度 subagent；主 session 对最终结论与交付负责。
本文件不能覆盖系统／开发者指令、用户当前指令或真实工具权限。
`AGENTS.md` 是自动加载入口，本文件保存完整约定，避免维护两套规则。

## 1. 主 session 与 subagent 分工

- **主 session** 是直接与用户交流的负责人，主要负责理解目标、划分任务、
  分配文件所有权、协调依赖、审核科学证据、简要汇报和整合交付。
- 对文献核查、物理实现、数据获取、数值实验和独立审查等实质性工作，
  优先划分为有清晰边界的 subagent 任务；简短、低开销的工作可直接完成。
  在用户已授权范围内调度，不为每次委派重复请求确认。
- 每次委派写清目标、必读来源、可修改文件、输入输出、验证标准和禁做事项。
  共享工作区内避免多个 agent 同时修改同一文件；并行只用于独立工作。
- 子 agent 承担详细阅读、计算与调试，向主 session 返回精简摘要：结论、
  修改路径、运行命令、实际数字、证据路径、限制和待决问题。
  原始日志、长公式推导、表格与图保存到仓库文件，避免整段回灌主 session。
- 科学关键项由不同 agent 或主 session 独立核查原式、符号、单位、
  时间尺度和参考数据独立性；主 session 不把子 agent 的结论直接当作证据。
- 主 session 用中文向用户汇报实测结果、变化、问题与下一步；持续工作期间
  按当前平台要求定期更新，不用大量原始工具输出来代替汇报。

## 2. 主 session 上下文预算

- 以**主 session 已用上下文不超过可用窗口的 30%**为工作预算，尽量保持
  主 session 精简，把深入实现和长实验交给 subagent。
- 先读本文件与短摘要，再按问题定位；使用 `rg`、定向读取、限制输出长度。
  不在主 session 一次性打印全部论文、JSON、CSV、源码或工具清单。
- 可用时给子 agent 最小必要背景，避免不必要的完整对话继承。
  检查 agent 结果时优先读取指标、差异与指向具体证据的路径。
- 平台提供上下文用量读数时按实际读数管理；没有读数时，不声称已精确
  测量或保证比例，采取保守的短摘要与文件交接方式。
- 接近预算时先将已完成内容、运行中任务、文件所有权、实际结果、阻碍和
  下一步写入 `docs/WORK_LOG.md` 或专项报告，再向用户提出新 session 交接。
  用户要求新 session 时，在提交后使用真实的新任务／新会话工具；
  不把继承全部历史的 fork 或 subagent 冒充新的主 session。

## 3. 研究范围与科学约束

- 仓库远程：<https://github.com/ccsvccs/cislunar-timescale-experiments>。
  本地目录名可能不同；以 `git remote -v` 和实际工作区为准。
- 当前只研究**地球质心 TT↔TDB**。未经用户扩展，不实现月球坐标时、
  地面钟或航天器钟。时间差统一为 `Δ = TDB − TT`，单位秒。
- 目标是独立复现 JPL DE440 时间星历的数值一致性，目标残差小于 1 ps；
  这不是绝对物理准确度小于 1 ps 的声明。
- Phase 1：保留 Turyshev arXiv:2507.22145v2 式 (21) 的地心截断模型。
  Phase 2：单独实现 Park 2021 式 (3)–(7)，逐项记录来源与模型限制。
- 原始边界固定为 IAU 定义的 `T0 = JD 2443144.5003725` 加
  `TDB0 = −65.5 μs`。使用定义的 `L_B`、`L_G` 和两部分 JD，不能用
  经验 `L_C` 替换精确定义常数而不恢复线性项。
- 使用真实几何星历、SI、TDB-compatible 状态与配套 GM。区分地心、EMB、
  行星系统质心，禁止重复计质量、混用尺度或默默退回 toy 轨道。
- 地球轨道已受小天体摄动，不代表时间 integrand 已包含其直接 `GM/r` 势。
  新物理项必须查原文，不猜系数、质量、距离或参数时间规律。
- JPL `de440t.bsp` 仅作独立对照；核心模型不能读取其时间差作为输入、
  修正表或校准点。SOFA/ERFA 同样仅作参考。
- 分别报告瞬时速率 residual、原始时间 residual、数值积分误差、物理模型
  缺项和输入星历不确定度。接受指标禁止拟合零点或线性趋势；诊断拟合
  必须单列，不能替代原始残差。独立求积相符不代表物理模型完整。

## 4. 代码与数据入口

| 路径 | 用途 |
|---|---|
| `src/tt_tdb/constants.py`、`timebase.py` | 定义常数、固定历元、两部分 JD |
| `src/tt_tdb/ephemeris.py`、`potential.py`、`pn_rate.py` | 原 Phase 1 DE440 reader、外部势和式 (21) |
| `src/tt_tdb/integration.py`、`transforms.py` | 独立 Romberg／Gauss、补偿累计、缓存、双向转换 |
| `src/tt_tdb/extended_ephemeris.py` | 同一 SPK Type 2 位置多项式的解析速度和加速度 |
| `src/tt_tdb/jpl_model.py` | Park 分项模型、开关、太阳 J2、PN 项及公开环几何候选 |
| `src/tt_tdb/small_body.py` | 真实 SB441 小天体直接势的独立、明确标记的诊断 |
| `src/tt_tdb/jpl_reference.py`、`reference.py` | 独立 JPL 时间段与 ERFA 对照 |
| `tests/`、`scripts/` | 自动测试、下载与复现实验 |
| `data/README.md`、`manifest.json`、`phase2_small_body_provenance.md` | 数据来源、哈希、覆盖与版本限制 |
| `results/tables/`、`results/figures/` | 实测 JSON／CSV／日志与图；Phase 2 用 `phase2_` 前缀 |

`de440.bsp` 提供轨道，`de440t.bsp` 的特殊段提供 TT−TDB 参考，两者不可混淆。
特殊时间段的 X 单位是秒，不是 km；本项目取负号并正确转换解析导数的
day／second 单位。大型 BSP、积分缓存和虚拟环境不提交 Git。

`sb441-n373s.bsp` 包含 343 小行星和 30 KBO 的真实日心状态，GM 与
DE440 header 逐项匹配，但轨道是以 DE441 行星为参考另行拟合的产品。
它是版本混合诊断，不能称作已取得原 DE440 小天体轨道。

## 5. Reference 与必读文档

新 session 先读 [README](README.md)、本文件和最新专项报告，再按任务读取：

- [原任务说明](guild.md)、[Phase 1 模型审计](docs/model_audit.md)、
  [Phase 1 报告](docs/tt_tdb_eq21_report.md)。
- [Phase 2 目标](docs/phase2_jpl_time_ephemeris.md)、
  [逐项数学／代码映射](docs/phase2_model_mapping.md)、
  [Phase 2 实测报告](docs/phase2_jpl_reproduction_report.md)。
- [Reference 索引](docs/references.md)、[来源记录](docs/source_notes.md)、
  [工作交接记录](docs/WORK_LOG.md)。
- Turyshev v2：[已归档 PDF](docs/sources/turyshev_2507.22145v2.pdf)，
  核查式 (13)、(21) 及 Appendix A；根目录 `paper.pdf` 不是固定 v2 的替代品。
- Park et al. 2021：[实际 15 页 PDF](docs/sources/park_2021_de440_de441.pdf)，
  Sec. 2.3 式 (3)–(7) 及相关动力学章节；
  [DOI](https://doi.org/10.3847/1538-3881/abd414)，保留 CC BY 4.0 归属信息。
- Lu, Yang & Xie：[LTE440 v1 PDF](docs/sources/lte440_2509.18511v1.pdf)，
  Section 3 和 Figure 1；区分作者的地心数值一致性与月球物理精度。
- Farnocchia：[SB441 原说明](docs/sources/sb441_perturbers.pdf)，
  核查 373 体标识、GM 和另行拟合／积分的来源，保留原版权与署名。

新增来源应记录原链接、版本、下载验证、SHA256、原式与使用范围；
不能将 HTML 错误页当 PDF，也不能将未读资料写成已核实。

## 6. 已有实测状态与下一步

以下为 2026-10-14、实现提交 `67b8fd0` 的历史结果，后续改动后需重新验证。
比较区间为 2024-01-01 至 2026-01-01，含端点共 732 个每日 TDB 历元：

| 模型／诊断 | 最大绝对原始时间 residual |
|---|---:|
| Phase 1 原模型 | 34.585021534 ns |
| 主天体 Park 扩展项与公开环直接势 | 11.282037595 ns |
| 再加入 SB441 实际小天体直接势，太阳赤道轴 | 2.818949569 ps |
| 另按 Park 字面黄道纬度定义的诊断 | 1.953519690 ps |

**1 ps 目标尚未达到**；最佳诊断 RMS 为 1.950451572 ps。
历史全套测试为 61 passed；模拟缺核时为 53 passed、8 skipped，不能当作
未来 session 已运行的测试。没有通过调零、趋势拟合或物理参数拟合减差。

下一步优先查证原 DE440 小天体状态、环实际动力学／几何、Park 黄道纬度与
Folkner 2014 太阳赤道纬度的实现约定，以及 `GMSDOT` 的单位、定义历元和
是否进入时间势。不能从字段名或 residual 反推参数规律。
剩余约 2 ps 尚未唯一归因，不以猜测修正追求过关；1100 年完整积分未验收。

## 7. 工程、验证与 Git 约定

- 修改前检查 `pwd`、`git status --short --branch`、`git remote -v`。
  不覆盖用户工作，不强制 reset；同步采用适当的 fast-forward-only 操作。
- Python-first，复用已有可靠模块；保持 Phase 1 接口和结果可复现。
  新物理项独立分解、提供开关，真实缺数据要明确报错或跳过，不伪造通过。
- 物理或输入改变时维护模型指纹与缓存失效，记录实际源码／数据哈希。
  不为整理格式随意修改已测物理源码，使结果与代码来源失去对应。
- 依据变更运行相关测试；提交科学实现前运行全套测试。纯文档变更检查
  路径、链接和 diff 即可，不要求重跑耗时的全部数值积分。
- 结果写清命令、模型配置、积分方法、数据版本、比较区间、原始指标和限制。
  不把少量历史采样或全覆盖速率采样称作长时段累计验收。
- 用户要求 commit 时只暂存本任务文件，检查 staged diff，保留其它改动。
  按会话授权提交／推送，不把先前任务的权限无条件扩展到无关工作。

常用命令（真实数据缺失时先查看 `data/README.md`）：

```sh
.venv/bin/python -m pytest -q
.venv/bin/python scripts/download_data.py --phase2-small-bodies --verify-only
.venv/bin/python scripts/compare_phase2_jpl.py
.venv/bin/python scripts/check_phase2_direct_integrand.py
.venv/bin/python scripts/compare_phase2_small_bodies.py
.venv/bin/python scripts/compare_phase2_solar_axis.py
.venv/bin/python scripts/phase2_history_profile.py
```

`--phase2-small-bodies` 的实际下载约 982 MB，按任务需要启用；不提交二进制核。
新 session 无需重做已验证的所有工作：先确认 Git 状态、阅读对应实测表和
未解决事项，再执行用户最新指派的下一项任务。
