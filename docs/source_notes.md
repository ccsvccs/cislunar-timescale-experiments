# 原始资料、版本与实际阅读证据

查阅日期：2026-10-13。除用户提供的 v3 PDF 外，下面都是本次实际打开或下载的原始资料。不会把搜索摘要、模型记忆或未下载的文件写成已读。数值公式审计先于积分器建立，见 [model_audit.md](model_audit.md)。

| 资料 | 原始链接 / 固定版本 | 本地证据与阅读范围 |
|---|---|---|
| Turyshev 主要复现论文 | [arXiv 2507.22145v2 PDF](https://arxiv.org/pdf/2507.22145v2) | `sources/turyshev_2507.22145v2.pdf`、`sources/turyshev_v2.txt`；Sec II.B.1–3、II.C、Appendix A.1/A.3；视觉读取 PDF 第7页的式21/22，截图 `sources/eq21_v2_page7.png` |
| 用户提供论文 | [arXiv v3](https://arxiv.org/pdf/2507.22145v3) | 根目录 `paper.pdf`，首页实际标记 v3 15 Nov 2025，文本 `sources/paper_local.txt`；SecII、AppendixA与v2核心公式一致，不能将它标成v2 |
| Turyshev HTML | [v2 HTML](https://arxiv.org/html/2507.22145v2) | 网页标记v2，但正文日期显示2026，附录编号161/162，官方PDF为A43/A44；以固定PDF核对原始符号 |
| Lu, Yang & Xie | [arXiv 2509.18511v1 PDF](https://arxiv.org/pdf/2509.18511v1)、[Section3 HTML](https://arxiv.org/html/2509.18511#S3) | `sources/lte440_2509.18511v1.pdf`、`sources/lte440_v1.txt`；已读Section3 Algorithm，额外检查式7含向量势与混合势；v2 URL实际404，查询历史只有v1，故保存v1 |
| Park, Folkner, Williams & Boggs (2021), DE440/DE441 | [JPL 官方 PDF](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf)、[DOI](https://doi.org/10.3847/1538-3881/abd414) | **参考链接已归档，PDF 尚未作为二进制文件提交**；具体获取和验证步骤见 [Park 2021 来源记录](sources/park_2021_de440_de441.reference.md)。本阶段应阅读 Sec. 2.3 Eq. (3)–(7) 及轨道模型相关章节；不得假称已经核验 PDF 哈希。 |
| LTE440仓库 | [GitHub](https://github.com/xlucn/LTE440) | 本地HEAD `76ea3a97fdbbd645352d3855b91dd2e878221351`，工作树只README，实际阅读完整README；没有LICENSE，没有原始积分源码 |
| LTE440 Release | [20250620](https://github.com/xlucn/LTE440/releases/tag/20250620)、[资产清单](https://github.com/xlucn/LTE440/releases/expanded_assets/20250620) | `sources/lte440_release_assets.html`；实际资产为 `lte440-tclmtdb.zip`、`lte440-zh.pdf` 和源树自动归档；未下载无关月球数据 |
| IAU 2006 B3 | [官方原链接](https://www.iau.org/static/resolutions/IAU2006_Resol3.pdf)、[ESO IAU档案镜像](https://iauarchive.eso.org/static/resolutions/IAU2006_Resol3.pdf) | `sources/iau2006_b3.pdf`、`sources/iau2006_b3.txt`；完整阅读定义式、注1/2/4/5；原网页工具失败后成功下载档案镜像 |
| IAU 定义常数 | [NSFA fixed constants](https://iau-a3.gitlab.io/NSFA/NSFA_fixed) | 已打开官方工作组表，核对LG/LB/TDB0性质；IAU的LG/LB为固定定义，而LC依赖星历 |
| TT原点与SI | [BIPM CCTF2017 Recommendation3](https://www.bipm.org/en/-/cctf-recommendation-2017-3)、[SI Brochure](https://www.bipm.org/en/publications/si-brochure) | 核对TAI原点与32.184s关系，光速采用SI精确定义 |
| ERFA dtdb | [原始源码](https://github.com/liberfa/erfa/blob/master/src/dtdb.c) | `sources/erfa_dtdb.c`；读参数/符号说明、Fairhead–Bretagnon、u=v=0、DE405约±3ns范围说明与最终求和代码；只作对照 |
| ERFA tttdb | [原始源码](https://github.com/liberfa/erfa/blob/master/src/tttdb.c) | `sources/erfa_tttdb.c`；读双部分JD加法的真实实现，未作为核心转换代码 |
| JPL DE440/441说明 | [JPL说明](https://ssd.jpl.nasa.gov/doc/de440_de441.html)、[JPL导出说明](https://ssd.jpl.nasa.gov/planets/eph_export.html) | 实际打开两页，记录SPICE、ASCII、binary不同产品，不假定数据分量一致 |
| DE440原始header | [header.440](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440) | `sources/de440_header.440`；读645个常数名、GROUP1050的3行15列指针；NCOEFF=1018，末尾第14/15列的阶数和分段数为0，标准ASCII该header没有TT−TDB系数 |
| JPL BSP README | [官方README](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/README.txt) | `../data/jpl_bsp_README.txt`；读TT−TDB伪目标1000000001/1000000000及X秒、Vx导数约定；该历史说明不代表每个BSP均含分量 |
| SPICE状态/时间 | [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)、[Time Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html) | 已打开官方资料，核对ET为TDB、J2000秒、6维状态及几何状态，禁止把SPICE简化时间函数当作高精度参考 |

LTE440 Section3 的可复用内容是积分组织思想：10阶Romberg、半天分段、减常数并最终恢复线性项，以及积分区间变化的误差测试。本文没有公开完整积分程序，本地也没有这些源码。其10fs与地球约1ps是作者原实验结果，与我们的式21截断模型及有限天体集合不能直接等同。数据/示例许可未明确，故本项目没有复制该仓库代码。

真实 JPL 时间参考的获取、段列表、哈希、SPK单位例外与实际核验由本项目数据模块记录在 [data/README.md](../data/README.md)；最终已经取得官方 `de440t.bsp`，并实际检查到时间段 `(center,target)=(1000000000,1000000001)`。不能用标准 `de440.bsp` 的无时间段读取失败替代对完整官方目录的调查，也不能靠一个猜测文件名宣称取得参考。

## 下载证据与哈希

| 文件 | SHA256 |
|---|---|
| 用户 `paper.pdf` | `c4d8ef3d9dc0a592585eaa985e84a434240fdb3632f13d47bcda38371c5db79c` |
| Turyshev v2 PDF | `dbbd9b74ee55f8bd1383d39d4b9aec505dca4a8fbb31531319ae5e4b110ece01` |
| LTE440 v1 PDF | `b6a620c0388f01b5349cbafdc2c7febc7ec8db43336b49dc65293ddaba2f72ac` |
| IAU B3 PDF | `1c8da5a20e4d67213f6008286cf133e0a495762719c4067b665f158ef83059da` |
| ERFA dtdb.c 快照 | `228bd7c7c1c88dda44c29a076f289c968f64d77a2dc389138b03ab5d1d442105` |
| ERFA tttdb.c 快照 | `452fa3604d8963b8c454d1b0c84aa4330ae4d08bc895d631b41217ecbb6863ef` |
| DE440 header | `0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359` |

ERFA源码下载自`master`，采用本地内容哈希固定证据，不伪造不存在的commit固定版本。源码自带原许可与版权说明。LTE440资产网页公布zip的SHA256为`baa74b47ef42633297617a1ed507daeb2157eb7a9a53bda26b5d5cc951b7b6dc`；本任务没有下载zip，因此这只是**服务器公布摘要**，不是本地验签结果。

壳层网络首次被沙箱拒绝，按权限规则使用批准的下载命令完成小型原始资料获取。GitHub API返回403，改为实际公开Release资产页；失败尝试不写成读取成功。CALCEPH网页工具未成功打开，故不声称已阅读其具体API；本项目实际采用SPK Type2独立读取路径。
