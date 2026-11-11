"""地心TT–TDB全覆盖独立积分：月球静态J2+C22，排除太阳直接spin。"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
import _common
import numpy as np
from tt_tdb.constants import A_SCALE,TDB0,T0_TDB
from tt_tdb.full_earth_model import FullEarthModel,anchored_prefix
from tt_tdb.integration import _quadrature_intervals
from tt_tdb.jpl_reference import JPLTimeEphemeris
from tt_tdb.small_body import file_sha256


PREFIX='results/tables/full_span_earth_'
ANCHOR=T0_TDB.seconds_since_j2000()


def metrics(values):
    d=np.asarray(values,float)
    return {'count':int(d.size),'min':float(d.min()),'max':float(d.max()),
            'max_abs':float(abs(d).max()),'mean':float(d.mean()),
            'rms':float(np.sqrt(np.mean(d*d))),'peak_to_peak':float(np.ptp(d))}


def integrate_grid(rate,grid,signature,label,method,order,channels,batch,no_cache):
    """Every cached increment is a direct model quadrature, never a fit/table."""
    config={'model':signature,'label':label,'method':method,'order':order,'channels':channels,
            'physical_anchor_seconds':ANCHOR,'grid_shape':list(grid.shape),
            'grid_sha256':hashlib.sha256(grid.tobytes()).hexdigest(),
            'accumulation':'Neumaier, separately toward each coverage boundary from anchor split'}
    digest=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    directory=Path('results/cache/full_span_earth')/digest;directory.mkdir(parents=True,exist_ok=True)
    metadata_path=directory/'metadata.json';data_path=directory/'increments.npy'
    metadata={'fingerprint':digest,'configuration':config,'completed_intervals':0,
              'evaluations':0,'quadrature_wall_seconds':0.,'sum_error_indicators_seconds':[0.]*channels}
    valid=False
    if metadata_path.exists() and data_path.exists() and not no_cache:
        old=json.loads(metadata_path.read_text())
        valid=old['fingerprint']==digest
        if valid:metadata=old
    increments=np.lib.format.open_memmap(data_path,mode='r+' if valid else 'w+',dtype='float64',shape=(len(grid)-1,channels))
    began=time.perf_counter();last_progress=began
    for start in range(metadata['completed_intervals'],len(grid)-1,batch):
        end=min(start+batch,len(grid)-1);t=time.perf_counter()
        result=_quadrature_intervals(rate,grid[start:end],grid[start+1:end+1],method,order,channels)
        increments[start:end]=result.value;increments.flush()
        metadata['completed_intervals']=end;metadata['evaluations']+=result.evaluations
        metadata['quadrature_wall_seconds']+=time.perf_counter()-t
        metadata['sum_error_indicators_seconds']=(np.asarray(metadata['sum_error_indicators_seconds'])+A_SCALE*result.error.sum(axis=0)).tolist()
        temporary=metadata_path.with_suffix('.tmp');temporary.write_text(json.dumps(metadata));temporary.replace(metadata_path)
        if time.perf_counter()-last_progress>30 or end==len(grid)-1:
            print(f'{label}: {end}/{len(grid)-1} intervals, this run {time.perf_counter()-began:.1f}s',flush=True)
            last_progress=time.perf_counter()
    origin=int(np.searchsorted(grid,ANCHOR,side='right')-1)
    if not 0<=origin<len(grid):raise ValueError('Physical anchor must lie within model coverage')
    boundary=_quadrature_intervals(rate,np.array([ANCHOR]),grid[origin:origin+1],method,order,channels)
    values=anchored_prefix(increments,origin,boundary.value[0])
    metadata.update(this_run_wall_seconds=time.perf_counter()-began,
                    numerical_split_grid_seconds=float(grid[origin]),
                    physical_anchor_to_split_integral_seconds=boundary.value[0].tolist(),
                    anchor_boundary_evaluations=boundary.evaluations,
                    increment_payload_sha256=file_sha256(data_path))
    return values,metadata


def convergence(model):
    """Three widely separated 32-day windows, independent of reference data."""
    outputs={}
    starts_jd=[model.start_jd+128,2460310.5,model.end_jd-160]
    for group,rate,n,configs in (
        ('fast',model.fast_channels,12,[('gauss_legendre',24),('gauss_legendre',12),('romberg',6)]),
        ('small',model.small_channels,2,[('gauss_legendre',16),('gauss_legendre',4),('romberg',5)])):
        windows=[]
        for jd in starts_jd:
            starts=(jd-2451545.)*86400.+np.arange(32)*86400.;reference=None;rows={}
            for method,order in configs:
                began=time.perf_counter()
                value=_quadrature_intervals(rate,starts,starts+86400.,method,order,n).value
                if reference is None:reference=value
                rows[f'{method}{order}_day']={'max_interval_difference_seconds':float(A_SCALE*abs(value-reference).max()),
                    'max_window_cumulative_difference_seconds':float(A_SCALE*abs(np.cumsum(value-reference,axis=0)).max()),
                    'wall_seconds':time.perf_counter()-began}
            windows.append({'start_jd':jd,'results':rows})
        outputs[group]=windows
    return outputs


def write_gzip(path,headers,arrays):
    # 17 significant digits preserve every binary64 input round-trip.
    with gzip.open(path,'wt',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(headers)
        for rows in zip(*arrays):writer.writerow([format(float(x),'.17g') for x in rows])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-cache',action='store_true');parser.add_argument('--benchmark-only',action='store_true')
    args=parser.parse_args();began=time.perf_counter()
    with FullEarthModel() as model:
        signature=model.fingerprint()
        signature['experiment_source_sha256']={str(p):file_sha256(p) for p in
            (Path(__file__),Path('tests/test_full_earth_model.py'),Path('src/tt_tdb/integration.py'))}
        signature['software']=_common.versions()
        short=convergence(model);_common.write_json(PREFIX+'short_convergence.json',short)
        if args.benchmark_only:
            print(json.dumps(short,indent=2));return
        # Reference identity/coverage is read here; no reference value enters
        # any integrand, quadrature, grid increment or anchor offset.
        with JPLTimeEphemeris() as reference:
            start_jd=max(model.start_jd,reference.start_jd);end_jd=min(model.end_jd,reference.end_jd)
            reference_identity={'path':str(reference.path),'sha256':file_sha256(reference.path),
                                'segment':[1000000000,1000000001]}
        count=int(round(end_jd-start_jd))+1
        if end_jd-start_jd != count-1:raise ValueError('Expected integer-day common coverage')
        times=(start_jd-2451545.)*86400.+np.arange(count)*86400.
        coarse_indices=np.unique(np.r_[np.arange(0,count,4),count-1]);coarse=times[coarse_indices]
        fast,fast_info=integrate_grid(model.fast_channels,times,signature,'fast_GL12_day','gauss_legendre',12,12,512,args.no_cache)
        small,small_info=integrate_grid(model.small_channels,times,signature,'SB441_GL4_day','gauss_legendre',4,2,512,args.no_cache)
        independent_fast,if_info=integrate_grid(model.fast_channels,times,signature,'fast_Romberg6_day','romberg',6,12,512,args.no_cache)
        independent_small,is_info=integrate_grid(model.small_channels,coarse,signature,'SB441_Romberg5_4day','romberg',5,2,256,args.no_cache)
        # Apply the original offset once. Extended precision here only limits
        # summation rounding; the integrands and stored outputs are binary64.
        park_delta=TDB0+A_SCALE*fast[:,:9].astype(np.longdouble).sum(axis=-1)
        small_delta=A_SCALE*small.astype(np.longdouble).sum(axis=-1)
        lunar_delta=A_SCALE*(fast[:,9].astype(np.longdouble)+fast[:,10])
        delta=np.asarray(park_delta+small_delta+lunar_delta,float)
        axis_delta=A_SCALE*fast[:,11]
        independent_delta=np.asarray(TDB0+A_SCALE*(independent_fast[coarse_indices,:11].astype(np.longdouble).sum(axis=-1)
                                                   +independent_small.astype(np.longdouble).sum(axis=-1)),float)
        quadrature_time=delta[coarse_indices]-independent_delta
        fast_difference=A_SCALE*(fast-independent_fast)
        small_difference=A_SCALE*(small[coarse_indices]-independent_small)
        rate=np.empty(count);axis_rate=np.empty(count)
        for start in range(0,count,1024):
            end=min(start+1024,count);value=model(times[start:end])
            rate[start:end]=A_SCALE*value[:,model.primary_indices].astype(np.longdouble).sum(axis=-1)
            axis_rate[start:end]=A_SCALE*value[:,model.axis_index]
        # Only after all actual model integrals exist, evaluate the reference.
        with JPLTimeEphemeris() as reference:
            ref_delta,ref_rate=reference.evaluate(times)
            anchor_ref,_=reference.evaluate(ANCHOR)
        residual=delta-ref_delta;rate_residual=rate-ref_rate
        secondary=residual+axis_delta;secondary_rate=rate_residual+axis_rate
        modern_times=757339200.+np.arange(732)*86400.
        modern_indices=np.searchsorted(times,modern_times)
        np.testing.assert_array_equal(times[modern_indices],modern_times)
        old_sb=np.genfromtxt('results/tables/phase2_sb441_time_series.csv',delimiter=',',names=True)
        old_lunar=np.genfromtxt('results/tables/phase2_lunar_multipole_modern.csv',delimiter=',',names=True)
        np.testing.assert_array_equal(old_sb['tdb_j2000_seconds'],modern_times)
        np.testing.assert_array_equal(old_lunar['tdb_j2000_seconds'],modern_times)
        old_expected=old_sb['sb441_delta_s']+old_lunar['primary_J2_C22_delta_s']
        reproduction=metrics(delta[modern_indices]-old_expected)
        result={'command':' '.join([sys.executable,*sys.argv]),'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'configuration_chinese':'Park主要体+公开ring leading势+真实SB441小体leading势+月球静态J2+C22；不加太阳直接spin；太阳J2默认太阳赤道轴',
            'model_signature':signature,'coverage_jd':[start_jd,end_jd],'daily_sample_count':count,
            'independent_cumulative_sample_count':len(coarse),'daily_sampling':'All endpoints and every TDB midnight; no decimation',
            'anchor_tdb_seconds':ANCHOR,'anchor_jd_parts':T0_TDB.jd_parts(),'anchor_delta_seconds':TDB0,
            'anchor_model_minus_JPL_seconds':TDB0-anchor_ref,'numerical_anchor_handling':'Exact integral split at nearby grid node; both directions Neumaier accumulation; no reanchoring',
            'sign':'Delta=TDB-TT; residual=model-JPL; rate=dDelta/dTDB',
            'reference':reference_identity,'short_convergence':short,
            'methods':{'primary_fast':fast_info,'primary_small_body':small_info,'independent_fast':if_info,'independent_small_body':is_info},
            'primary_raw_time_residual_seconds':metrics(residual),'primary_raw_rate_residual_dimensionless':metrics(rate_residual),
            'secondary_ecliptic_axis_time_residual_seconds':metrics(secondary),'secondary_ecliptic_axis_rate_residual_dimensionless':metrics(secondary_rate),
            'independent_full_span_cumulative_time_difference_seconds':metrics(quadrature_time),
            'independent_fast_max_channel_difference_seconds':abs(fast_difference).max(axis=0).tolist(),
            'independent_small_body_max_channel_difference_seconds':abs(small_difference).max(axis=0).tolist(),
            'modern_732_raw_time_residual_seconds':metrics(residual[modern_indices]),
            'modern_732_raw_rate_residual_dimensionless':metrics(rate_residual[modern_indices]),
            'modern_previous_SB441_plus_lunar_reproduction_seconds':reproduction,
            'modern_archived_input_sha256':{n:file_sha256('results/tables/'+n) for n in
                        ('phase2_sb441_time_series.csv','phase2_lunar_multipole_modern.csv')},
            'lunar_static_endpoint_contribution_seconds':[float(lunar_delta[0]),float(lunar_delta[-1])],
            'acceptance':{'threshold_seconds':1e-12,'passed_raw_daily_samples':bool(abs(residual).max()<1e-12),
                          'zero_point_fitted':False,'linear_trend_fitted':False,
                          'absolute_physical_accuracy_claimed':False},
            'wall_seconds_before_output':time.perf_counter()-began}
        write_gzip(PREFIX+'daily.csv.gz',
             ['tdb_j2000_seconds','model_TDB_minus_TT_s','JPL_TDB_minus_TT_s','raw_time_residual_s','model_rate','JPL_rate',
              'raw_rate_residual','secondary_ecliptic_time_residual_s','secondary_ecliptic_rate_residual',
              'lunar_static_J2_C22_delta_s','SB441_leading_delta_s','ecliptic_minus_equator_delta_s'],
             (times,delta,ref_delta,residual,rate,ref_rate,rate_residual,secondary,secondary_rate,lunar_delta,small_delta,axis_delta))
        write_gzip(PREFIX+'independent_cumulative.csv.gz',
             ['tdb_j2000_seconds','primary_delta_s','independent_delta_s','primary_minus_independent_s',
              'asteroid_integral_difference_s','KBO_integral_difference_s'],
             (coarse,delta[coarse_indices],independent_delta,quadrature_time,small_difference[:,0],small_difference[:,1]))
        _common.csv_table(PREFIX+'modern.csv',['tdb_j2000_seconds','raw_time_residual_s','raw_rate_residual','previous_reproduction_difference_s'],
                           np.column_stack((modern_times,residual[modern_indices],rate_residual[modern_indices],delta[modern_indices]-old_expected)))
        make_plot(times,residual,rate_residual,secondary,secondary_rate,start_jd,end_jd)
        result['wall_seconds_total']=time.perf_counter()-began
        result['output_sha256']={str(p):file_sha256(p) for p in
             (Path(PREFIX+'daily.csv.gz'),Path(PREFIX+'independent_cumulative.csv.gz'),Path(PREFIX+'modern.csv'),
              Path('results/figures/full_span_earth_residuals.png'))}
        _common.write_json(PREFIX+'summary.json',result)
        write_report(result)
        print(json.dumps({k:result[k] for k in ('primary_raw_time_residual_seconds','primary_raw_rate_residual_dimensionless',
                     'independent_full_span_cumulative_time_difference_seconds','modern_previous_SB441_plus_lunar_reproduction_seconds','wall_seconds_total')},indent=2))


def make_plot(times,residual,rate,secondary,secondary_rate,start_jd,end_jd):
    directory=_common.make_plots_directory()
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font_path=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if not font_path.exists():raise FileNotFoundError('Chinese font unavailable: '+str(font_path))
    font=FontProperties(fname=str(font_path))
    year=2000.+times/(365.25*86400.)
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True,constrained_layout=True)
    axs[0].plot(year,rate/1e-22,linewidth=.65,color='#165a9f',label='主模型：太阳赤道 J2 轴')
    axs[1].plot(year,residual*1e12,linewidth=.65,color='#165a9f',label='主模型：太阳赤道 J2 轴')
    axs[0].plot(year,secondary_rate/1e-22,linewidth=.55,color='#b75a27',alpha=.65,label='副诊断：黄道 J2 轴')
    axs[1].plot(year,secondary*1e12,linewidth=.55,color='#b75a27',alpha=.7,label='副诊断：黄道 J2 轴')
    axs[0].set_ylabel('原始速率残差 / $10^{-22}$',fontproperties=font)
    axs[1].set_ylabel('原始时间残差 / ps',fontproperties=font)
    axs[1].set_xlabel('TDB 年份（近似公历年）',fontproperties=font)
    fig.suptitle('地心 TT–TDB 全覆盖残差：月球静态 J2+C22，不加太阳直接 spin\n'
                 f'JD {start_jd:g}–{end_jd:g}；原始 1977 锚点；每日采样；未拟合',fontproperties=font,fontsize=13)
    for ax in axs:
        ax.grid(alpha=.25);ax.axhline(0,color='black',linewidth=.4);ax.legend(prop=font,loc='best')
    fig.savefig(directory/'full_span_earth_residuals.png',dpi=190);plt.close(fig)


def write_report(r):
    t=r['primary_raw_time_residual_seconds'];v=r['primary_raw_rate_residual_dimensionless']
    q=r['independent_full_span_cumulative_time_difference_seconds'];modern=r['modern_732_raw_time_residual_seconds']
    text=f'''# 地心 TT–TDB 全覆盖独立积分报告

本轮模型采用 Park 主要天体与公开 ring leading 势、实际 SB441 小体 leading 势、
月球静态 J2+C22；不加入太阳直接 spin。太阳 J2 默认太阳赤道轴，黄道轴只作为
预先声明的副诊断，未依据残差选择轴。时间差统一为 Δ=TDB−TT，残差为 model−JPL。

## 覆盖与原始结果

真实共同覆盖为 JD {r['coverage_jd'][0]}–{r['coverage_jd'][1]}，
包含两端共 {r['daily_sample_count']} 个每日 TDB 样本。所有时间值均由真实模型
从原始 1977 锚点完整累计，未减常数或拟合趋势。

| 指标 | 实测 |
|---|---:|
| 全覆盖原始时间 residual 最大绝对值 | {t['max_abs']*1e12:.12g} ps |
| 全覆盖原始时间 residual RMS | {t['rms']*1e12:.12g} ps |
| 全覆盖原始时间 residual 最小/最大 | {t['min']*1e12:.12g} / {t['max']*1e12:.12g} ps |
| 全覆盖原始速率 residual 最大绝对值 | {v['max_abs']:.12g} |
| 全覆盖原始速率 residual RMS | {v['rms']:.12g} |
| 现代 732 日点原始时间 residual 最大绝对值 | {modern['max_abs']*1e12:.12g} ps |
| 独立完整累计最大差 | {q['max_abs']*1e12:.12g} ps |

每日样本小于 1ps 的验收状态为 **{r['acceptance']['passed_raw_daily_samples']}**。
这只是对指定 JPL 数值时间星历的采样一致性，不代表连续时间严格上界或绝对物理精度。

## 数值方法与独立验证

快通道（Phase1、Park 扩展、月球静态势与轴差）使用 GL12/day；
373 体 leading 势使用 GL4/day。独立方法为快通道 Romberg6/day、
小体 Romberg5/4day。独立方法也从同一物理锚点完整积到两个覆盖端点，
在全部 {r['independent_cumulative_sample_count']} 个共同四日节点比较累计时间，
不是用稀疏 rate 比较代替长期 time 验证。

物理锚点保持 T0+TDB0；因日网格未恰好经过该时刻，先直接求积
anchor→邻近日网格点，再分别向过去和未来作 Neumaier 补偿累计。
该严格积分拆分不重新定义锚点、不读取参考值作为边界。
每个日/四日区间均实际计算模型；缓存按源码、全部数据 SHA、方法、网格与锚点指纹隔离。
三处分离的 32 日窗口先与高阶 GL24/GL16 比较；实测与全覆盖比较见 JSON。
求积差异是经验数值诊断，不能充当物理模型完整性证明。

现代结果另与已有 SB441 + 静态 lunar J2/C22 表逐点复现，最大差
{r['modern_previous_SB441_plus_lunar_reproduction_seconds']['max_abs']*1e12:.12g} ps。
所有旧表 epochs 均逐点核对；旧表哈希与当前完整模型指纹均记录在 JSON。

## 来源与限制

物理映射沿用 [Phase2 映射](phase2_model_mapping.md) 与
[月球静态势来源](../data/phase2_lunar_multipole_provenance.md)。旧物理源码未改。
SB441 是以 DE441 为参考单独拟合的真实小体轨道，不能称为原 DE440 小体状态。
公开 ring 几何、原 ring 动力学状态及小体/ring 完整 PN 项的限制继续保留；
月球仅保留静态 J2/C22，不含时变 tide/rotation distortion。
GMSDOT 时间势约定也未由本轮解决。JPL 时间数据仅在独立对照阶段读取。

副诊断黄道轴全覆盖最大时间 residual 为
{r['secondary_ecliptic_axis_time_residual_seconds']['max_abs']*1e12:.12g} ps；
该数值不用于替换预先指定的默认轴或宣称验收。

## 复现与产物

```sh
.venv/bin/python scripts/compare_full_span_earth.py
```

实际总耗时 {r['wall_seconds_total']:.2f}s；依赖版本、真实命令、输入/源码/产物 SHA
和方法统计保存在 `results/tables/full_span_earth_summary.json`。
原始每日序列为 `full_span_earth_daily.csv.gz`，无损 gzip 与 17 位有效数字；
独立累计序列为 `full_span_earth_independent_cumulative.csv.gz`。
图 `results/figures/full_span_earth_residuals.png` 使用全部每日数据，
上为 rate residual，下为 time residual，均保持原始零点及趋势。
'''
    Path('docs/full_span_earth_report.md').write_text(text)


if __name__=='__main__':main()
