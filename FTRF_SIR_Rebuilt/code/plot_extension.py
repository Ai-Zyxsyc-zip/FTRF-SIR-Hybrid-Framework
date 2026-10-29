"""Publication figures from saved CSV only; no training modules imported."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/extension'; FIG=ROOT.parent/'latex/English/figures'
STYLE={'RF':('#555555','o','-'),'Transfer':('#3876a8','s','--'),
       'SIR_residual':('#b4872c','^','-.'),'FTRF-SIR':('#964875','D','-'),
       'Transfer_persistence':('#718043','v',':')}
NAMES={'RF':'Target RF','Transfer':'Transfer RF','SIR_residual':'Local SIR residual',
       'FTRF-SIR':'FTRF-SIR','Transfer_persistence':'Transfer persistence residual'}
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})

def save(fig,name):
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight'); plt.close(fig)

def main():
    s=pd.read_csv(OUT/'summary.csv'); p=pd.read_csv(OUT/'paired.csv')
    fig,axs=plt.subplots(2,3,figsize=(11,6),sharex=True)
    for j,ds in enumerate(['ILI','CLS','CPT']):
        d=s[(s.dataset==ds)&(s.interval=='Gaussian')]
        for method,(color,marker,line) in STYLE.items():
            z=d[d.method==method].sort_values('budget')
            for i,metric in enumerate(['rmse','wis']):
                axs[i,j].plot(z.budget,z[metric],color=color,marker=marker,ls=line,label=NAMES[method],lw=1.4)
                axs[i,j].grid(axis='y',alpha=.2)
        axs[0,j].set_title(ds); axs[1,j].set_xlabel('Target training sample size')
        axs[1,j].set_xticks([24,52,104])
    axs[0,0].set_ylabel('RMSE (percentage points)'); axs[1,0].set_ylabel('WIS (percentage points)')
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,1.055))
    fig.tight_layout(); save(fig,'ext_learning')
    methods=['Pool','Transfer','SIR_feature','Transfer_feature','SIR_residual','FTRF-SIR','Persistence_residual','Transfer_persistence','Half_SIR_RF']
    labels=['Pooled RF','Transfer RF','SIR feature','Transfer + SIR feature','Local SIR residual','FTRF-SIR','Local persistence residual','Transfer persistence residual','Half SIR / half RF']
    fig,axs=plt.subplots(1,3,figsize=(11,4.6),sharey=True)
    for j,ds in enumerate(['ILI','CLS','CPT']):
        d=p[(p.dataset==ds)&(p.budget==52)&(p.interval=='Gaussian')&(p.metric=='wis')&(p.reference=='RF')].set_index('method').loc[methods]
        axs[j].errorbar(d.delta,np.arange(len(d)),xerr=np.vstack([d.delta-d.lo95,d.hi95-d.delta]),fmt='o',color='#3876a8',capsize=2,ms=4)
        axs[j].axvline(0,color='#444444',lw=.8); axs[j].set_title(ds); axs[j].set_xlabel('WIS difference from target RF')
        axs[j].set_yticks(np.arange(len(d)),labels); axs[j].grid(axis='x',alpha=.2)
    axs[0].invert_yaxis(); fig.tight_layout(); save(fig,'ext_ablation')
    fig,axs=plt.subplots(1,3,figsize=(11,3.8))
    for j,ds in enumerate(['ILI','CLS','CPT']):
        for method,(color,marker,line) in STYLE.items():
            if method not in ['RF','FTRF-SIR','SIR_residual']: continue
            d=s[(s.dataset==ds)&(s.budget==52)&(s.method==method)].set_index('interval').loc[['Gaussian','Order']]
            axs[j].plot(d.width95,100*d.coverage95,color=color,lw=1,label=NAMES[method])
            for kind,m in [('Gaussian','o'),('Order','s')]:
                axs[j].scatter(d.loc[kind,'width95'],100*d.loc[kind,'coverage95'],marker=m,color=color,s=32,facecolors=color if kind=='Gaussian' else 'white',zorder=3)
        axs[j].axhline(95,color='#444444',ls='--',lw=.8); axs[j].set_title(ds)
        axs[j].set_xlabel('Mean 95% width (percentage points)'); axs[j].grid(alpha=.15)
    axs[0].set_ylabel('Empirical 95% coverage (%)')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,1.06))
    fig.tight_layout(); save(fig,'ext_intervals')
    fig,axs=plt.subplots(1,3,figsize=(10,3.4))
    for j,ds in enumerate(['ILI','CLS','CPT']):
        d=s[(s.dataset==ds)&(s.budget==52)&(s.interval=='Gaussian')].set_index('method').loc[['FTRF-SIR_near','FTRF-SIR','FTRF-SIR_far']]
        axs[j].bar(np.arange(3),d.wis,color=['#3876a8','#964875','#b4872c'],width=.6)
        axs[j].set_xticks(np.arange(3),['Nearest 2','All 5','Farthest 2']); axs[j].set_title(ds); axs[j].grid(axis='y',alpha=.2)
    axs[0].set_ylabel('WIS (percentage points)'); fig.tight_layout(); save(fig,'ext_sources')
    h=pd.read_csv(OUT/'horizons.csv'); fig,axs=plt.subplots(1,3,figsize=(10,3.6))
    for j,ds in enumerate(['ILI','CLS','CPT']):
        for method,(color,marker,line) in STYLE.items():
            d=h[(h.dataset==ds)&(h.budget==52)&(h.interval=='Gaussian')&(h.method==method)].sort_values('horizon')
            axs[j].plot(d.horizon,d.rmse,color=color,marker=marker,ls=line,label=NAMES[method]); axs[j].grid(axis='y',alpha=.2)
        axs[j].set_xticks([1,2,4]); axs[j].set_title(ds); axs[j].set_xlabel('Forecast horizon (weeks)')
    axs[0].set_ylabel('RMSE (percentage points)')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.5,1.12))
    fig.tight_layout(); save(fig,'ext_horizons')
    print('Five figure families exported from CSV')

if __name__=='__main__': main()
