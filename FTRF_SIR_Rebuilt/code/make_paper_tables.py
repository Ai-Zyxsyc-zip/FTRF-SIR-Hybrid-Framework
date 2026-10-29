"""Generate all manuscript numbers from saved CSV/JSON. No fitting imports."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; T=ROOT/'results/tables'; O=ROOT/'results/tables'
ORDER=['Persistence','Seasonal','SIR','RF_target','RF_pool','FTRF','SIR_RF','FTRF-SIR','Guarded']
def esc(x): return str(x).replace('_',r'\_').replace('%',r'\%')
def dual(en,cn): return r'\ifdefined\CN '+cn+r'\else '+en+r'\fi'
def table(caption,label,headers,rows,long=False):
    env='longtable' if long else 'tabular'; spec='l'+'r'*(len(headers)-1)
    prefix=(r'\begin{center}\small\begin{longtable}{'+spec+r'}\caption{'+caption+r'}\label{'+label+r'}\\' if long else r'\begin{table}[tbp]\centering\small\caption{'+caption+r'}\label{'+label+r'}\begin{tabular}{'+spec+'}')
    return prefix+'\n'+r'\toprule '+' & '.join(headers)+r'\\\midrule'+'\n'+('\n'.join(' & '.join(map(esc,row))+r'\\' for row in rows))+'\n'+r'\bottomrule\end{'+env+'}'+(r'\end{center}' if long else r'\end{table}')+'\n'

def main():
    O.mkdir(exist_ok=True); a=pd.read_csv(T/'forecast_summary.csv'); p=pd.read_csv(T/'predictions.csv')
    counts=p.groupby('dataset').size()/9; n=len(p)//9
    english=[f'The completed benchmark contains {n} admissible forecast tasks and {len(p):,} model forecasts. ILI and CLS each contribute {int(counts.ILI)} tasks; CPT contributes {int(counts.CPT)}, with {144-int(counts.CPT)} prespecified tasks excluded for data availability. All nine models are scored on the same tasks within each dataset.']
    chinese=[f'完整基准包含 {n} 个有效预测任务、{len(p):,} 条模型预测。ILI 与 CLS 各有 {int(counts.ILI)} 个任务，CPT 有 {int(counts.CPT)} 个任务，因数据可用性排除 {144-int(counts.CPT)} 个设计任务。每组数据内九个模型的评价任务完全一致。']
    for dataset in ['ILI','CLS','CPT']:
        d=a[a.dataset==dataset].set_index('method'); x=d.loc['FTRF-SIR']; b=d.loc['RF_target']; best=d.wis.idxmin()
        rm=100*(1-x.rmse/b.rmse); wi=100*(1-x.wis/b.wis)
        english.append(f'For {dataset}, the FTRF-SIR model RMSE is {x.rmse:.3f} percentage points versus {b.rmse:.3f} for target RF, a {rm:.1f}\\% reduction. Its WIS is {x.wis:.3f} versus {b.wis:.3f}, a {wi:.1f}\\% reduction. The lowest WIS among all nine models is achieved by {esc(best)} ({d.loc[best,"wis"]:.3f}). FTRF-SIR 95\\% interval coverage is {100*x.coverage95:.1f}\\%, below nominal coverage.')
        chinese.append(f'{dataset} 上，混合模型 RMSE 为 {x.rmse:.3f} 个百分点，目标森林为 {b.rmse:.3f}，相对下降 {rm:.1f}\\%；WIS 分别为 {x.wis:.3f} 与 {b.wis:.3f}，相对下降 {wi:.1f}\\%。九种模型中 WIS 最小者为 {esc(best)}（{d.loc[best,"wis"]:.3f}）。混合模型 95\\% 区间覆盖为 {100*x.coverage95:.1f}\\%，低于名义水平。')
        rows=[[m,f'{d.loc[m,"rmse"]:.3f}',f'{d.loc[m,"mae"]:.3f}',f'{d.loc[m,"wis"]:.3f}',f'{100*d.loc[m,"coverage95"]:.1f}',f'{d.loc[m,"width95"]:.3f}'] for m in ORDER]
        tab=table(dual(f'{dataset}: accuracy and interval quality; errors and widths are percentage points.',f'{dataset}：精度与区间质量；误差、宽度单位为百分点。'),'tab:'+dataset,['Model','RMSE','MAE','WIS',r'Cov. (\%)','Width'],rows)
        english.append(tab); chinese.append(tab)
    english.append('The FTRF-SIR model is therefore better than target RF in this benchmark, but not uniformly best. Pooled RF leads on ILI WIS and RMSE; target-only SIR residual RF has slightly lower CLS RMSE than the full FTRF-SIR model. The guarded variant does not improve the full FTRF-SIR model\'s aggregate WIS on any of the three datasets. These comparisons suggest that the mechanism offset contributes substantially, while the added value of source residual transfer depends on the outcome and scoring rule. Undercoverage is a remaining problem, not evidence of calibrated 95\\% uncertainty.')
    chinese.append('因此，混合模型在本基准中优于目标森林，但并非普遍最优。合并森林在 ILI 的 WIS 和 RMSE 上领先；目标 SIR 残差森林的 CLS RMSE 略低于完整混合。防护版本在三个数据集上的总体 WIS 均未优于完整混合。这提示机制偏置贡献较大，来源残差迁移的额外价值依赖观测类型与评分规则。覆盖不足仍是问题，不能据此声称得到校准的 95\\% 不确定性。')
    pairs=pd.read_csv(T/'paired_differences.csv'); ds=pairs[(pairs.method=='FTRF-SIR')&(pairs.metric=='wis')]
    for _,r in ds.iterrows():
        english.append(f'The equal-cluster FTRF-SIR model-minus-target-RF WIS difference for {r.dataset} is {r.mean_difference:.3f}, with descriptive 95\\% interval [{r.lo95:.3f}, {r.hi95:.3f}].')
        chinese.append(f'{r.dataset} 中，按簇等权计算的“混合减目标森林”WIS 差为 {r.mean_difference:.3f}，描述性 95\\% 区间为 [{r.lo95:.3f}, {r.hi95:.3f}]。')
    pm=pd.read_csv(T/'posterior_metrics.csv'); s=pm.groupby(['scenario','prior','method']).mean(numeric_only=True)
    for scenario,cn in [('matched','同源'),('shifted','失配')]:
        local=s.loc[(scenario,'Local','Grid')]; trans=s.loc[(scenario,'Transfer','Grid')]
        english.append(f'In the {scenario}-source simulation, grid-reference forecast RMSE averages {local.rmse:.3f} under the local prior and {trans.rmse:.3f} under the transfer prior; corresponding NLL values are {local.nll:.3f} and {trans.nll:.3f}. Each average uses ten independently generated target/source datasets. These figures do not establish a universal benefit from prior transfer.')
        chinese.append(f'{cn}来源模拟中，网格参考预测 RMSE 在局部先验下平均为 {local.rmse:.3f}，迁移先验下为 {trans.rmse:.3f}，相应 NLL 为 {local.nll:.3f} 和 {trans.nll:.3f}。每项平均来自十组独立目标与来源数据，不能据此建立先验迁移的普遍收益。')
    mh=pm[pm.method=='MH']; la=pm[pm.method=='Laplace']
    english.append(f'Across the 40 posterior targets, MH has maximum conventional split-$\\widehat R$ {mh.rhat.max():.3f} and minimum ESS {mh.ess.min():.0f}. Mean standardized marginal Wasserstein discrepancy is {mh.wasserstein_std.mean():.3f} for MH and {la.wasserstein_std.mean():.3f} for Laplace. These are fixed-target approximation diagnostics, separate from forecast skill and from the choice of prior.')
    chinese.append(f'40 个后验目标中，MH 最大传统分链 $\\widehat R$ 为 {mh.rhat.max():.3f}，最小 ESS 为 {mh.ess.min():.0f}。平均标准化边缘 Wasserstein 差异为 MH {mh.wasserstein_std.mean():.3f}、Laplace {la.wasserstein_std.mean():.3f}。这些是固定统计目标下的近似诊断，与预测能力和先验选择分开解释。')
    for name,text in [('results_en',english),('results_cn',chinese)]: (O/f'{name}.tex').write_text('\n\n'.join(text),encoding='utf-8')
    dq=pd.read_csv(T/'data_quality.csv')
    rows=[[r.dataset,r.region,r.rows,r.missing,f'{r.first_year}--{r.last_year}'] for r in dq.itertuples()]
    (O/'quality.tex').write_text(table(dual('Raw snapshot audit for the six selected states.','六州原始快照审计。'),'tab:quality',['Dataset','Region','Rows','Missing','Years'],rows).replace('[tbp]','[H]'),encoding='utf-8')
    full=[]
    for dataset in ['ILI','CLS','CPT']:
        d=a[a.dataset==dataset].set_index('method')
        rows=[[m,f'{d.loc[m,"nll"]:.3f}',f'{d.loc[m,"crps"]:.3f}',f'{d.loc[m,"wis"]:.3f}',f'{d.loc[m,"coverage_error95"]:.3f}',f'{d.loc[m,"fit_seconds"]:.3f}',f'{d.loc[m,"mechanism_seconds"]:.3f}'] for m in ORDER]
        full.append(table(dual(f'{dataset}: probabilistic scores and accounting components. Forest and mechanism columns are seconds.',f'{dataset}：概率得分与成本。Forest、Mechanism 列单位为秒。'),'tab:full'+dataset,['Model','NLL','CRPS','WIS','Cov. error','Forest','Mechanism'],rows))
    (O/'full_metrics.tex').write_text('\n'.join(full),encoding='utf-8')
    rows=[]
    for (scenario,prior,method),r in s.iterrows(): rows.append([scenario,prior,method,f'{r.rmse:.3f}',f'{r.nll:.3f}',f'{r.wasserstein_std:.3f}',f'{r.mmd2:.4f}',f'{r.seconds:.3f}'])
    (O/'posterior.tex').write_text(table(dual('Mean results over ten independent datasets per source condition. W1 is standardized marginal distance.','每种来源条件十组独立数据的平均结果；W1 为标准化边缘距离。'),'tab:posterior',['Source','Prior','Method','RMSE','NLL','W1','MMD2','Seconds'],rows),encoding='utf-8')
    v=json.loads((ROOT/'audit/validation.json').read_text())
    err=v['solver_max_abs_error_8_16_32']; conv=pd.read_csv(T/'grid_convergence.csv')
    numerical=table(dual('Discrete solver refinement audit: maximum absolute compartment error.','离散求解细化审计：各仓室最大绝对误差。'),'tab:numerical',['Substeps per week','Maximum error'],[[k,f'{e:.6f}'] for k,e in zip([8,16,32],err)])
    numerical+=dual(f'Maximum coarse-to-fine posterior mean difference in the four checked grid cases is {conv.absolute_difference.max():.6f}; maximum boundary mass over all posterior targets is {pm.edge_mass.max():.2e}.',f'四个检查网格中，粗细网格后验均值的最大差为 {conv.absolute_difference.max():.6f}；全部目标最大边缘质量为 {pm.edge_mass.max():.2e}。')
    (O/'numerical.tex').write_text(numerical,encoding='utf-8')
    print('Generated synchronized manuscript tables and result narratives.')

if __name__=='__main__': main()
