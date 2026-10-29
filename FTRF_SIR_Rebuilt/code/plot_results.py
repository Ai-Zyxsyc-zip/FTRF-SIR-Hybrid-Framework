"""Rebuild EVERY figure from saved tables; does not import or execute models."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; T=ROOT/'results/tables'; F=ROOT/'results/figures'
ORDER=['Persistence','Seasonal','SIR','RF_target','RF_pool','FTRF','SIR_RF','FTRF-SIR','Guarded']
COLORS=dict(zip(ORDER,['#9a9a9a','#bababa','#d98e32','#4477aa','#66aacc','#aa4499','#669944','#228833','#cc6677']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})

def save(fig,name):
    F.mkdir(exist_ok=True); fig.savefig(F/f'{name}.pdf'); fig.savefig(F/f'{name}.png',dpi=180); plt.close(fig)

def main(posterior_only=False):
    if not posterior_only:
        a=pd.read_csv(T/'forecast_summary.csv'); p=pd.read_csv(T/'predictions.csv')
        fig,axes=plt.subplots(1,3,figsize=(12,4.8),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            d=a[a.dataset==dataset].set_index('method').reindex(ORDER)
            ax.barh(ORDER,d.wis,color=[COLORS[x] for x in ORDER]); ax.invert_yaxis(); ax.set_title(dataset); ax.set_xlabel('WIS (pp; lower is better)'); ax.grid(axis='x',alpha=.2)
        save(fig,'forecast_scores')
        h=pd.read_csv(T/'horizon_summary.csv')
        fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            for method in ['Persistence','RF_target','FTRF','FTRF-SIR','Guarded']:
                d=h[(h.dataset==dataset)&(h.method==method)]
                ax.plot(d.horizon,d.rmse,'o-',label=method,color=COLORS[method])
            ax.set(title=dataset,xlabel='Forecast horizon (weeks)',ylabel='RMSE (percentage points)',xticks=[1,2,4]); ax.grid(alpha=.2)
        axes[2].legend(fontsize=10); save(fig,'horizons')
        fig,axes=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            d=a[a.dataset==dataset].set_index('method').reindex(ORDER)
            ax.barh(ORDER,d.coverage95,color=[COLORS[x] for x in ORDER]); ax.axvline(.95,color='black',ls='--',lw=1); ax.set(xlim=(0,1),title=dataset,xlabel='95% interval coverage'); ax.invert_yaxis()
        save(fig,'calibration')
        delta=pd.read_csv(T/'paired_differences.csv')
        fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            d=delta[(delta.dataset==dataset)&(delta.metric=='wis')].set_index('method').reindex([m for m in ORDER if m!='RF_target'])
            ax.errorbar(d.mean_difference,np.arange(len(d)),xerr=np.vstack([d.mean_difference-d.lo95,d.hi95-d.mean_difference]),fmt='o',color='#4477aa',capsize=3)
            ax.set_yticks(np.arange(len(d)),d.index); ax.axvline(0,color='black',lw=1); ax.invert_yaxis(); ax.set(title=dataset,xlabel='WIS difference from RF_target')
        save(fig,'paired_differences')
        fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            # Fixed illustrative target and horizon, never chosen by model performance.
            d=p[(p.dataset==dataset)&(p.region=='California')&(p.horizon==4)]
            season=int(d.season.min()); d=d[d.season==season]
            for method in ['RF_target','FTRF-SIR','Guarded']:
                x=d[d.method==method].sort_values('origin')
                ax.plot(x.origin-x.origin.min(),x['mean'],'o-',color=COLORS[method],label=method)
                if method=='Guarded': ax.fill_between(x.origin-x.origin.min(),x.lo95,x.hi95,color=COLORS[method],alpha=.15)
            obs=d[d.method=='RF_target'].sort_values('origin'); ax.scatter(obs.origin-obs.origin.min(),obs.y,c='black',s=18,label='Observed',zorder=5)
            ax.set(title=f'{dataset}, California, {season}-{str(season+1)[-2:]}',xlabel='Weeks since first selected origin',ylabel='Outcome (%)'); ax.grid(alpha=.2)
        axes[-1].legend(fontsize=10); save(fig,'trajectories')
        fig,axes=plt.subplots(1,3,figsize=(11,4.5),layout='constrained')
        for ax,dataset in zip(axes,['ILI','CLS','CPT']):
            d=a[a.dataset==dataset].set_index('method').reindex(ORDER)
            ax.barh(ORDER,d.fit_seconds,label='Forest fitting / querying',color='#4477aa')
            ax.barh(ORDER,d.mechanism_seconds,left=d.fit_seconds,label='Attributed cached SIR work',color='#ddaa33')
            ax.set(title=dataset,xlabel='Accounted seconds / task'); ax.invert_yaxis()
        axes[-1].legend(fontsize=9); save(fig,'costs')
    s=pd.read_csv(T/'posterior_metrics.csv')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
    for ax,metric,label in zip(axes,['rmse','nll','beta_error'],['Forecast RMSE (pp)','Forecast NLL','Absolute beta error (per week)']):
        data=[]; labels=[]
        for scenario in ['matched','shifted']:
            for prior in ['Local','Transfer']:
                data.append(s[(s.scenario==scenario)&(s.prior==prior)&(s.method=='Grid')][metric]); labels.append(f'{scenario}\n{prior}')
        ax.boxplot(data,tick_labels=labels,showfliers=True); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    save(fig,'synthetic_transfer')
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
    for ax,metric,label in zip(axes,['wasserstein_std','mmd2','seconds'],['Standardized marginal W1','Biased MMD squared','Inference seconds']):
        d=s[s.method!='Grid']; data=[d[d.method==m][metric] for m in ['Laplace','MH']]
        ax.boxplot(data,tick_labels=['Laplace','MH']); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    save(fig,'posterior_fidelity')
    q=pd.read_csv(T/'posterior_samples.csv'); q=q[(q.scenario=='matched')&(q.rep==0)&(q.prior=='Transfer')]
    fig,axes=plt.subplots(1,2,figsize=(8,3.2),layout='constrained')
    for ax,parameter in zip(axes,['beta','gamma']):
        for method in ['Grid','Laplace','MH']:
            ax.hist(q[q.method==method][parameter],bins=45,density=True,histtype='step',label=method,lw=1.5)
        ax.set(xlabel=f'{parameter} (per week)',ylabel='Density')
    axes[-1].legend(); save(fig,'posterior_marginals')

if __name__=='__main__':
    import sys
    main(posterior_only='--posterior-only' in sys.argv)
