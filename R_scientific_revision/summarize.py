from pathlib import Path
import json, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; D=ROOT/'R_scientific_revision/results'
F=ROOT/'R_scientific_revision/figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
groups={'conditional':['task','model'],'forecast':['task','model'],'simulation':['condition','model']}
frames={};summaries={};counts={}
for name,keys in groups.items():
 d=pd.read_csv(D/f'{name}_predictions.csv');assert len(d)>0 and np.isfinite(d[['observed','predicted']]).all().all()
 if name!='simulation':assert (d.test_time>d['cutoff' if name=='conditional' else 'origin']).all()
 rowkeys=([*keys[:1],'target','cutoff','test_time'] if name=='conditional' else [*keys[:1],'target','origin','test_time','h'] if name=='forecast' else ['condition','replicate','origin','test_time','h'])
 assert not d.duplicated(rowkeys+['model']).any()
 assert d.groupby(rowkeys).model.nunique().eq(7 if name=='conditional' else 9).all()
 d['SE']=(d.predicted-d.observed)**2;d['AE']=abs(d.predicted-d.observed)
 s=d.groupby(keys).agg(MSE=('SE','mean'),MAE=('AE','mean'),N=('SE','size')).reset_index();s['RMSE']=np.sqrt(s.MSE)
 r2=d.groupby(keys).apply(lambda v:1-v.SE.sum()/((v.observed-v.observed.mean())**2).sum() if ((v.observed-v.observed.mean())**2).sum()>0 else np.nan,include_groups=False).rename('R2').reset_index()
 s=s.merge(r2,on=keys)
 s.to_csv(D/f'{name}_summary.csv',index=False);frames[name]=d;summaries[name]=s
 counts[name]={'prediction_rows':len(d),'evaluation_rows':len(d[rowkeys].drop_duplicates()),'models':d.model.nunique()}
 if name=='conditional':
  blocks=d.groupby(['task','target','cutoff','model']).SE.mean().reset_index();blocks.to_csv(D/'conditional_block_errors.csv',index=False)
  fig,axs=plt.subplots(2,3,figsize=(12,7));models=list(d.model.unique())
  for ax,task in zip(axs.flat,sorted(d.task.unique())):
   v=blocks[blocks.task==task];ax.boxplot([v[v.model==m].SE for m in models],tick_labels=models,showfliers=False)
   ax.tick_params(axis='x',rotation=60);ax.set_title(task);ax.set_ylabel('MSE (squared percentage points)' if task in ['ILI','CLS_B'] else 'MSE (squared counts)')
  axs.flat[-1].axis('off');fig.tight_layout();fig.savefig(F/'original_five_tasks.pdf');plt.close(fig)
 else:
  split=keys[0];fig,axs=plt.subplots(1,d[split].nunique(),figsize=(12,3.5),squeeze=False)
  for ax,k in zip(axs.flat,sorted(d[split].unique())):
   v=s[s[split]==k].sort_values('RMSE');ax.barh(v.model,v.RMSE);ax.set_title(k);ax.set_xlabel('RMSE (percentage points)');ax.invert_yaxis()
  fig.tight_layout();fig.savefig(F/f'{name}_comparison.pdf');plt.close(fig)
# Plot the saved SIR, RF, and FTRF-SIR trajectory comparisons.
d=frames['forecast'];z=d[(d.target=='Alabama')&(d.h==1)&d.model.isin(['SIR','RF_tar','FTRF-SIR'])]
fig,axs=plt.subplots(3,3,figsize=(12,8),sharex=False)
for i,task in enumerate(sorted(z.task.unique())):
 for j,m in enumerate(['SIR','RF_tar','FTRF-SIR']):
  v=z[(z.task==task)&(z.model==m)].sort_values('test_time');ax=axs[i,j];ax.plot(range(len(v)),v.observed,'o-',label='Observed',ms=3);ax.plot(range(len(v)),v.predicted,'s--',label='Predicted',ms=3);ax.set_title(f'{task}: {m}');ax.set_ylabel('Percentage points');ax.set_xlabel('Ordered forecast origin');ax.legend(fontsize=7)
fig.tight_layout();fig.savefig(F/'original_forecast_trajectories.pdf');plt.close(fig)
# The simulation uses independent replicate-level paired resampling.
rng=np.random.default_rng(20260921);paired=[]
d=frames['simulation'];by=d.groupby(['condition','replicate','model']).SE.mean().unstack('model')
for condition in ['matched','shifted']:
 for model in ['RF_tar','SIR_residual','Persistence_residual','Auxiliary_SIR_mean','All_SIR_mean']:
  v=(by.loc[condition,'FTRF-SIR']-by.loc[condition,model]).to_numpy();boot=rng.choice(v,(5000,len(v)),replace=True).mean(1)
  paired.append({'condition':condition,'comparator':model,'MSE_difference':v.mean(),'lower':np.quantile(boot,.025),'upper':np.quantile(boot,.975),'replicates':len(v),'improved_fraction':np.mean(v<0)})
pd.DataFrame(paired).to_csv(D/'simulation_paired_comparisons.csv',index=False)
def tex(s,cols):
 return s[cols].to_latex(index=False,float_format=lambda x:f'{x:.3f}',escape=True,longtable=True)
for lang in ['en','cn']:
 s=summaries['forecast'];h=s[s.model=='FTRF-SIR'].set_index('task');r=s[s.model=='RF_tar'].set_index('task');g=100*(1-h.RMSE/r.RMSE)
 sim=summaries['simulation'];hh=sim[sim.model=='FTRF-SIR'].set_index('condition');rr=sim[sim.model=='RF_tar'].set_index('condition');sg=100*(1-hh.RMSE/rr.RMSE)
 title='Repaired original tasks and additional R validation' if lang=='en' else '原任务修订与新增 R 验证'
 t='\\subsection{'+title+'}\n'
 if lang=='en':
  t+=f'The repaired five-task regression experiment produced {counts["conditional"]["prediction_rows"]:,} saved predictions. These conditional regression results retain the original targets without equating contemporaneous reconstruction with forecasting. Figure~\\ref{{fig:original_five_tasks}} displays the block-loss distributions in their native units. The additional R forecasting experiment contains {counts["forecast"]["evaluation_rows"]} common target--origin--horizon rows. Relative to target-only random forests, the FTRF-SIR model changed RMSE by {g["ILI"]:.1f}\\%, {g["CLS_B"]:.1f}\\%, and {g["CPT_B_rate"]:.1f}\\% reductions on ILI, clinical B positivity and historical B positivity, respectively; a negative reduction denotes deterioration. Auxiliary-SIR averages are retained as direct controls (Figure~\\ref{{fig:forecast_comparison}}). These are different outcomes and fitting rules from the principal total-positivity benchmark.\n\n'
  t+=f'The full-predictor simulation produced {counts["simulation"]["evaluation_rows"]} evaluation rows across 20 independent replicates and two source conditions. FTRF-SIR RMSE reductions relative to target forests were {sg["matched"]:.1f}\\% with matched sources and {sg["shifted"]:.1f}\\% with shifted sources. All competitors use the same generated observations. The paired replicate summaries in the supplement assess variability; they are exploratory and not adjusted for the multiple comparisons. The experiment validates executable prediction under controlled conditions, not the biological adequacy of the real-data observation model.\n'
 else:
  t+=f'修订后的五任务回归共保存 {counts["conditional"]["prediction_rows"]:,} 条预测，保留原目标但不把同期重建等同于未来预测。图~\\ref{{fig:original_five_tasks}} 按原单位展示测试块误差。新增 R 周预测含 {counts["forecast"]["evaluation_rows"]} 个共同地区--起点--步长任务。相对目标森林，FTRF-SIR 在 ILI、临床 B 型阳性率和历史 B 型阳性率上的 RMSE 降幅分别为 {g["ILI"]:.1f}\\%、{g["CLS_B"]:.1f}\\%、{g["CPT_B_rate"]:.1f}\\%；负降幅表示变差。辅助 SIR 平均作为直接对照保留（图~\\ref{{fig:forecast_comparison}}）。这些结局与拟合规则不同于主实验的总阳性率，不能合并计算收益。\n\n'
  t+=f'完整预测器模拟在 20 个独立重复和两种来源条件下产生 {counts["simulation"]["evaluation_rows"]} 个评价任务。相对目标森林，匹配与失配来源下 FTRF-SIR 的 RMSE 降幅分别为 {sg["matched"]:.1f}\\% 和 {sg["shifted"]:.1f}\\%。所有方法使用相同生成观测。附件按独立重复给出配对变异；这些比较是探索性的，未作多重比较调整。该实验验证可控条件下的预测程序，不证明真实监测观测模型的生物学充分性。\n'
 for stem,caption in [('original_five_tasks','Chronological conditional-regression block losses for the five original targets. Units differ across panels; overlapping blocks are not independent replicates.' if lang=='en' else '原五任务的时间留出条件回归误差；各面板单位不同，重叠测试块并非独立重复。'),('forecast_comparison','Additional R forecasts, with auxiliary-region SIR averages retained as controls.' if lang=='en' else '新增 R 周预测，包含辅助地区 SIR 平均对照。'),('simulation_comparison','Full-predictor synthetic comparison under matched and shifted sources.' if lang=='en' else '完整预测器在来源匹配与失配条件下的模拟比较。')]:
  t+='\\fig{'+stem+'}{1}{'+caption+'}\n'
 (D/f'paper_results_{lang}.tex').write_text(t,encoding='utf8')
 t='\\section{'+('Repaired R experiment tables' if lang=='en' else '修订 R 实验完整结果')+'}\n'
 t+=('RMSE and MAE use percentage points for ILI and CLS B, and counts for conditional CPT B and PHL tasks; all forecasting and synthetic outcomes use percentage points. N counts evaluation rows, including overlapping conditional blocks, not independent replicates. R2 is the coefficient of determination and can be negative. The R forest uses 80 trees, minimum node size 3, and floor(p/3) candidate variables (at least one). Its annual phase uses the ordered row index divided by 52, not exact calendar weeks.\n' if lang=='en' else 'ILI 与临床 B 型回归的 RMSE、MAE 单位为百分点，CPT B 与 PHL 条件回归为计数；周预测和模拟均为百分点。N 是评价行数，包含重叠条件测试块，不代表独立重复数。R2 为决定系数，可为负。R 森林使用 80 棵树、最小叶节点 3、候选变量数 floor(p/3)（至少 1 个）。年度相位按排序后的行号除以 52 构造，并非精确日历周。\n')
 for name,s in summaries.items():
  cols=groups[name]+['RMSE','MAE','R2','N'];tab=tex(s,cols)
  t+='\\subsection{'+name.replace('_',' ')+'}\n{\\small\n'+tab+'}\n'
 t+='\\subsection{'+('Paired simulation differences' if lang=='en' else '模拟配对差异')+'}\n'
 t+=('The difference is FTRF-SIR minus comparator MSE (squared percentage points); negative favors FTRF-SIR. Bounds are descriptive 95\\% percentile intervals from 5,000 paired resamples of 20 independent replicates, without multiplicity adjustment.\n' if lang=='en' else '差值为 FTRF-SIR 减对照的 MSE（平方百分点），负值有利于 FTRF-SIR。上下限由 20 个独立重复的 5,000 次配对自助抽样产生，为未作多重比较调整的描述性 95\\% 百分位区间。\n')
 t+=tex(pd.DataFrame(paired),['condition','comparator','MSE_difference','lower','upper'])
 t+='\\clearpage\n'
 t+='\\fig{original_forecast_trajectories}{1}{'+('One-week-ahead Alabama trajectories at the saved origins; all three methods use the same observations.' if lang=='en' else 'Alabama 保存起点的一周预测轨迹；三种方法使用相同观测。')+'}\n'
 (D/f'paper_appendix_{lang}.tex').write_text(t,encoding='utf8')
for folder in ['English','Chinese']:
 dest=ROOT/'latex'/folder/'figures'
 import shutil
 for p in F.glob('*.pdf'):shutil.copy2(p,dest/p.name)
(D/'validation_summary.json').write_text(json.dumps(counts,indent=2),encoding='utf8')
print(json.dumps(counts,indent=2));print(summaries['forecast'].to_string(index=False));print(summaries['simulation'].to_string(index=False))

# Generate publication figures from the saved numerical results.
import subprocess, sys
subprocess.run([sys.executable, str(ROOT/"R_scientific_revision/plot_paper.py")], check=True)
