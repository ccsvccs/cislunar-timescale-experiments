# 地心 TT↔TDB：Turyshev 式 (21) 的独立 DE440 积分

本项目使用真实 DE440 几何星历和配套 GM，直接积分地球轨道速度与外部势的 $c^{-2}$、$c^{-4}$ 贡献。核心不调用 SOFA/ERFA；ERFA 仅用于对照。研究范围是地球质心，非地面站或任意航天器。

## 安装和数据

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python scripts/download_data.py
```

本次工作环境为 Python 3.14；其他符合 `pyproject.toml` 的版本也可安装。当前环境的 `.venv` 已安装依赖，无需重新创建。本次系统未提供 pip/ensurepip，因此首次安装使用了 `venv --without-pip` 加官方 `https://bootstrap.pypa.io/get-pip.py` 引导 pip；已有 pip 的常规 Python 可直接使用上述命令。下载脚本从 NASA/JPL 官方源取得完整 DE440 和同套GM，校验信息见 [数据说明](data/README.md)。大型文件被 Git 忽略。可用 `DE440_BSP`、`DE440_GM` 指定本地路径。

## 复现实验

```sh
.venv/bin/python -m pytest -q
.venv/bin/python scripts/convergence_study.py
.venv/bin/python scripts/run_experiment.py
.venv/bin/python scripts/compare_sofa.py
.venv/bin/python scripts/compare_jpl.py
.venv/bin/python scripts/error_budget.py
```

默认先完成短区间验证，再从1977定义历元累计积分，输出2024–2026两年曲线。两种独立算法、积分参数、数据哈希、耗时和版本记录在 `results/tables/`。`results/cache/` 是可复用累计积分缓存，不提交Git；缓存签名不一致时拒绝使用。

没有星历时，单元测试照常运行，真实 DE440 集成测试明确跳过。需要数据的实验会报错，不会自动替换为简化轨道。

## 阅读顺序

1. [理论与常数审计](docs/model_audit.md)：原始公式、定义常数与版本差异。
2. [中文教学实验报告](docs/tt_tdb_eq21_report.md)：公式、真实代码、测试结果和图表。
3. [实验设计](docs/experiment_plan.md)：事先固定的比较方式和边界条件。

输出符号统一为 `TDB−TT`（秒）。对外接口接受 `TwoPartJD`，时间差单独保存；不能用两个完整 float64 JD 相减验证皮秒级往返。数值积分一致性与物理模型绝对准确度在报告中分别讨论。

## 核心代码

| 文件 | 内容 |
|---|---|
| [constants.py](src/tt_tdb/constants.py) | 定义常数、Decimal审计、1977锚点 |
| [timebase.py](src/tt_tdb/timebase.py) | 分量JD与小量时间算术 |
| [ephemeris.py](src/tt_tdb/ephemeris.py) | 真实DE440地心状态、SI适配 |
| [potential.py](src/tt_tdb/potential.py) / [pn_rate.py](src/tt_tdb/pn_rate.py) | 外部势、F2/F4、诊断 |
| [integration.py](src/tt_tdb/integration.py) | Romberg、Gauss-Legendre、补偿累计和缓存验签 |
| [transforms.py](src/tt_tdb/transforms.py) | TDB−TT、直接和隐式双向转换 |
| [reference.py](src/tt_tdb/reference.py) / [jpl_reference.py](src/tt_tdb/jpl_reference.py) | ERFA与真实JPL参考比较 |
| [validation.py](src/tt_tdb/validation.py) | 明确标记的解析椭圆轨道测试模型 |

冷启动结果见 `results/tables/experiment_summary_fresh.json`，后续缓存复现见 `experiment_summary.json`。测试日志为 `test_results.log` 和 `test_without_kernels.log`。六张图包括年度变化、ERFA/JPL差异、PN阶次贡献、收敛与解析轨道验证。
