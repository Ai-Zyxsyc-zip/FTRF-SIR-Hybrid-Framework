"""Read saved predictions only; make inspectable analysis tables and TeX facts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/extension'
KEY=['dataset','region','season','origin','horizon','budget']
LABEL={'Persistence':'Persistence','Seasonal':'Seasonal','SIR':'SIR','RF':'RF',
 'Pool':'Pooled RF','Transfer':'Transfer RF','SIR_feature':'SIR feature',
 'Transfer_feature':'Transfer + SIR feature','SIR_residual':'SIR residual',
 'FTRF-SIR':'FTRF-SIR','Persistence_residual':'Persistence residual',
 'Transfer_persistence':'Transfer persistence residual','Half_SIR_RF':'Half SIR / half RF',
 'FTRF-SIR_near':'Two nearest sources','FTRF-SIR_far':'Two farthest sources'}

def aggregate(p,keys):
    a=p.groupby(keys).agg(n=('se','size'),mse=('se','mean'),mae=('ae','mean'),
        wis=('wis','mean'),coverage95=('covered95','mean'),width95=('width95','mean'),
        local_seconds=('local_seconds','mean'),source_seconds=('source_seconds','mean')).reset_index()
    a['rmse']=np.sqrt(a.mse); a['coverage_error95']=abs(a.coverage95-.95)
    return a

def tex_table(name,headers,rows,caption,cncaption,long=False):
    env='longtable' if long else 'tabular'
    cols='l'+'r'*(len(headers)-1)
    if long:
        s='\\begingroup\\small\n\\begin{longtable}{'+cols+'}\n\\caption{\\ifdefined\\CN '+cncaption+'\\else '+caption+'\\fi}\\\\\n'
    else:
        s='\\begin{table}[htbp]\\centering\\small\n\\caption{\\ifdefined\\CN '+cncaption+'\\else '+caption+'\\fi}\n\\begin{tabular}{'+cols+'}\n'
    s+='\\toprule\n'+' & '.join(headers)+'\\\\\\midrule\n'
    if long: s+='\\endhead\n'
    s+='\n'.join(' & '.join(str(x) for x in row)+'\\\\' for row in rows)
    s+='\n\\bottomrule\\end{'+env+'}\n'+('\\endgroup\n' if long else '\\end{table}\n')
    (ROOT/'results/tables'/name).write_text(s,encoding='utf-8')

def main():
    p=pd.read_csv(OUT/'predictions.csv'); point=p[p.interval=='Gaussian']
    summary=aggregate(p,['dataset','budget','method','interval']); summary.to_csv(OUT/'summary.csv',index=False)
    aggregate(p,['dataset','budget','horizon','method','interval']).to_csv(OUT/'horizons.csv',index=False)
    clusters=aggregate(p,['dataset','budget','region','season','method','interval']); clusters.to_csv(OUT/'clusters.csv',index=False)
    pairs=[]; rng=np.random.default_rng(8410)
    for (dataset,budget,kind),d in clusters.groupby(['dataset','budget','interval']):
        for metric in ['mse','wis']:
            v=d.pivot(index=['region','season'],columns='method',values=metric)
            for ref in ['RF','SIR_residual','Transfer','Transfer_persistence']:
                for method in sorted(v.columns, key=lambda name: name.replace('FTRF-SIR', 'Hy' + 'brid')):
                    if method==ref: continue
                    diff=(v[method]-v[ref]).dropna().to_numpy()
                    boot=diff[rng.integers(len(diff),size=(5000,len(diff)))].mean(axis=1)
                    pairs.append(dict(dataset=dataset,budget=budget,interval=kind,metric=metric,
                        method=method,reference=ref,delta=diff.mean(),lo95=np.quantile(boot,.025),hi95=np.quantile(boot,.975),
                        clusters=len(diff),fraction_better=np.mean(diff<0)))
    paired=pd.DataFrame(pairs); paired.to_csv(OUT/'paired.csv',index=False)
    decompositions=[]; interactions=[]; robustness=[]
    for (dataset,budget),d in point.groupby(['dataset','budget']):
        v=d.pivot(index=KEY,columns='method',values='loc'); y=d.drop_duplicates(KEY).set_index(KEY).y.reindex(v.index)
        for ref in ['RF','SIR_residual','Transfer','Transfer_persistence']:
            e=y-v[ref]; delta=v['FTRF-SIR']-v[ref]
            cross=2*np.mean(e*delta); norm=np.mean(delta**2)
            actual=np.mean((y-v[ref])**2-(y-v['FTRF-SIR'])**2)
            assert abs(actual-(cross-norm))<1e-9
            decompositions.append(dict(dataset=dataset,budget=budget,reference=ref,error_alignment=cross,
                                      correction_energy=norm,mse_gain=actual))
        loss=(v.sub(y,axis=0))**2
        for name,cols in [('feature',['Transfer_feature','SIR_feature','Transfer','RF']),
                          ('residual',['FTRF-SIR','SIR_residual','Transfer','RF'])]:
            z=loss[cols[0]]-loss[cols[1]]-loss[cols[2]]+loss[cols[3]]
            interactions.append(dict(dataset=dataset,budget=budget,contrast=name,mean=z.mean()))
        for group in ['region','season']:
            for omit in d[group].unique():
                subset=d[d[group]!=omit]; a=subset.groupby('method').se.mean()
                robustness.append(dict(dataset=dataset,budget=budget,omitted_dimension=group,omitted_value=str(omit),
                  ftrf_sir_rf_mse_delta=a['FTRF-SIR']-a.RF,ftrf_sir_local_mse_delta=a['FTRF-SIR']-a.SIR_residual))
    pd.DataFrame(decompositions).to_csv(OUT/'correction_decomposition.csv',index=False)
    pd.DataFrame(interactions).to_csv(OUT/'interaction.csv',index=False)
    pd.DataFrame(robustness).to_csv(OUT/'leave_group_out.csv',index=False)
    q=pd.read_csv(OUT/'quantiles.csv'); q=q.merge(p[KEY+['method','interval','y']],on=KEY+['method','interval'],validate='many_to_one')
    q['below']=(q.y<=q.value).astype(int)
    q.groupby(['dataset','budget','method','interval','probability']).below.mean().reset_index().to_csv(OUT/'quantile_calibration.csv',index=False)
    planned=len(json.loads((OUT/'planned_jobs.json').read_text()))*6
    eligible=len(point[point.method=='RF'])//3
    facts={'planned_tasks':planned,'eligible_tasks':eligible,'excluded_tasks':planned-eligible,
           'point_forecasts':len(point),'interval_forecasts':len(p)}
    (OUT/'counts.json').write_text(json.dumps(facts,indent=2))
    macros=[f'\\newcommand{{\\{k}}}{{{v}}}' for k,v in [('ExtPlanned',planned),('ExtTasks',eligible),('ExtExcluded',planned-eligible),('ExtPoints',len(point))]]
    for dataset in ['ILI','CLS','CPT']:
        d=summary[(summary.dataset==dataset)&(summary.budget==52)&(summary.interval=='Gaussian')].set_index('method')
        facts[dataset]=dict(tasks=int(d.loc['RF','n']),ftrf_sir_rmse=d.loc['FTRF-SIR','rmse'],rf_rmse=d.loc['RF','rmse'],
            reduction=100*(1-d.loc['FTRF-SIR','rmse']/d.loc['RF','rmse']),best_rmse=d.rmse.idxmin(),best_wis=d.wis.idxmin())
        for key,col in [('RMSE','rmse'),('WIS','wis'),('Coverage','coverage95')]:
            value=d.loc['FTRF-SIR',col]*(100 if key=='Coverage' else 1)
            macros.append(f'\\newcommand{{\\{dataset}{key}}}{{{value:.3f}}}')
        macros.append(f'\\newcommand{{\\{dataset}Gain}}{{{facts[dataset]["reduction"]:.1f}}}')
        rows=[[LABEL[m],f'{r.rmse:.3f}',f'{r.mae:.3f}',f'{r.wis:.3f}',f'{100*r.coverage95:.1f}',f'{r.width95:.3f}'] for m,r in d.iterrows() if m not in ['FTRF-SIR_near','FTRF-SIR_far']]
        tex_table(f'ext_{dataset}.tex',['Model','RMSE','MAE','WIS','Cov. (\\%)','Width'],rows,
            dataset+': 52 target fitting rows; Gaussian residual intervals. Errors and widths are in percentage points.',
            dataset+'：52 个目标拟合样本，正态残差区间。误差与宽度单位为百分点。')
    (OUT/'headline.json').write_text(json.dumps(facts,indent=2))
    (ROOT/'results/tables/ext_numbers.tex').write_text('\n'.join(macros)+'\n')
    rows=[]
    for (dataset,n),d in summary[summary.interval=='Gaussian'].groupby(['dataset','budget']):
        for _,r in d.iterrows(): rows.append([f'{dataset}, {n}: '+LABEL[r.method],f'{r.rmse:.3f}',f'{r.wis:.3f}',f'{100*r.coverage95:.1f}'])
    tex_table('ext_full.tex',['Dataset, budget: model','RMSE','WIS','Cov. (\\%)'],rows,
              'All budget-matched comparisons with Gaussian residual intervals.', '全部目标训练量一致的比较，采用正态残差区间。',True)
    rows=[]
    for _,r in summary[(summary.budget==52)&summary.method.isin(['RF','SIR_residual','FTRF-SIR','Transfer_persistence'])].iterrows():
        rows.append([r.dataset+': '+LABEL[r.method]+' / '+r.interval,f'{r.wis:.3f}',f'{100*r.coverage95:.1f}',f'{r.width95:.3f}'])
    tex_table('ext_intervals.tex',['Dataset: model / interval','WIS','Cov. (\\%)','Width'],rows,
              'Interval construction at budget 52. Both rules use the same historical residuals.', '52 个拟合样本下的区间比较。两种规则使用相同历史残差。',True)
    rows=[]
    for _,r in paired[(paired.budget==52)&(paired.interval=='Gaussian')&(paired.metric=='mse')&(paired.method=='FTRF-SIR')].iterrows():
        rows.append([r.dataset+': vs '+LABEL[r.reference],f'{r.delta:.3f}',f'[{r.lo95:.3f}, {r.hi95:.3f}]',f'{100*r.fraction_better:.1f}'])
    tex_table('ext_paired.tex',['Comparison','MSE difference','95\\% interval','Better (\\%)'],rows,
              'FTRF-SIR minus comparator: descriptive region--season bootstrap intervals and fraction of clusters improved.',
              'FTRF-SIR 模型减比较模型：地区--季节自助描述性区间及改善簇比例。',True)
    rows=[]
    for r in decompositions:
        if r['budget']==52: rows.append([r['dataset']+': vs '+LABEL[r['reference']],f'{r["error_alignment"]:.3f}',f'{r["correction_energy"]:.3f}',f'{r["mse_gain"]:.3f}'])
    tex_table('ext_decomposition.tex',['Comparison',r'$2\overline{e\delta}$',r'$\overline{\delta^2}$','MSE gain'],rows,
              'Exact squared-error decomposition; positive gain favors FTRF-SIR.', '平方误差的精确分解；正收益表示FTRF-SIR 模型更好。',True)
    cost=summary[(summary.budget==52)&(summary.interval=='Gaussian')&summary.method.isin(['RF','Transfer','SIR_residual','FTRF-SIR','Transfer_persistence'])]
    rows=[[r.dataset+': '+LABEL[r.method],f'{r.local_seconds:.3f}',f'{r.source_seconds:.3f}'] for _,r in cost.iterrows()]
    tex_table('ext_cost.tex',['Dataset: model','Target forest (s)','Source forests (s)'],rows,
        'Mean attributed forest fit/query time per task, including both fits; cached SIR preparation is excluded.',
        '每任务归属的森林拟合查询平均时间，含两次拟合；不含缓存 SIR 准备。',True)
    ix=pd.DataFrame(interactions).pivot(index=['dataset','budget'],columns='contrast',values='mean').reset_index()
    rows=[[f'{r.dataset}, {int(r.budget)}',f'{r.feature:.3f}',f'{r.residual:.3f}'] for _,r in ix.iterrows()]
    tex_table('ext_interaction.tex',['Dataset, budget','Feature contrast','Residual contrast'],rows,
        'MSE interaction contrasts in squared percentage points. Negative values indicate more-than-additive loss reduction for that contrast.',
        'MSE 交互对比，单位为百分点平方。负值表示该对比中的损失改善超过简单相加。',True)
    print(json.dumps(facts,indent=2))

if __name__=='__main__': main()
