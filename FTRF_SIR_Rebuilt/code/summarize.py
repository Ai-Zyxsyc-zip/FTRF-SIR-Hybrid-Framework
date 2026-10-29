"""Summaries, clustered paired differences, ranks and workbook from saved CSV only."""
import numpy as np
import pandas as pd
from prepare_data import ROOT

def summary(d,keys):
    a=d.groupby(keys,observed=True).agg(n=('se','size'),mse=('se','mean'),mae=('ae','mean'),nll=('nll','mean'),crps=('crps','mean'),wis=('wis','mean'),coverage95=('covered95','mean'),width95=('width95','mean'),fit_seconds=('fit_seconds','mean'),mechanism_seconds=('mechanism_seconds','mean'),accounted_seconds=('accounted_seconds','mean')).reset_index()
    a['rmse']=np.sqrt(a.mse); a['coverage_error95']=abs(a.coverage95-.95)
    return a

def main():
    t=ROOT/'results/tables'; p=pd.read_csv(t/'predictions.csv')
    a=summary(p,['dataset','method']); a.to_csv(t/'forecast_summary.csv',index=False)
    summary(p,['dataset','horizon','method']).to_csv(t/'horizon_summary.csv',index=False)
    summary(p,['dataset','region','season','method']).to_csv(t/'cluster_summary.csv',index=False)
    rows=[]; rng=np.random.default_rng(987)
    cluster=p.groupby(['dataset','region','season','method'])[['se','wis','nll']].mean().reset_index()
    for dataset,d in cluster.groupby('dataset'):
        for metric in ['se','wis','nll']:
            pivot=d.pivot(index=['region','season'],columns='method',values=metric)
            for method in pivot.columns:
                if method=='RF_target': continue
                delta=(pivot[method]-pivot.RF_target).to_numpy()
                boot=delta[rng.integers(0,len(delta),(5000,len(delta)))].mean(axis=1)
                rows.append(dict(dataset=dataset,method=method,metric=metric,mean_difference=delta.mean(),lo95=np.quantile(boot,.025),hi95=np.quantile(boot,.975),clusters=len(delta),fraction_better=np.mean(delta<0)))
    pd.DataFrame(rows).to_csv(t/'paired_differences.csv',index=False)
    ranks=[]
    for dataset,d in a.groupby('dataset'):
        for metric in ['rmse','mae','nll','crps','wis','coverage_error95','accounted_seconds']:
            for method,value,rank in zip(d.method,d[metric],d[metric].rank(method='average')):
                ranks.append(dict(dataset=dataset,method=method,metric=metric,value=value,rank=rank,direction='lower'))
    pd.DataFrame(ranks).to_csv(t/'metric_ranks.csv',index=False)
    pm=pd.read_csv(t/'posterior_metrics.csv')
    pm.groupby(['scenario','prior','method']).mean(numeric_only=True).reset_index().to_csv(t/'posterior_summary.csv',index=False)
    # Keep large posterior draws in CSV rather than embedding them in the workbook.
    with pd.ExcelWriter(ROOT/'results/results.xlsx',engine='openpyxl') as writer:
        for name in ['forecast_summary','horizon_summary','cluster_summary','paired_differences','metric_ranks','data_quality','posterior_summary','posterior_metrics','grid_convergence','predictions','synthetic_observations','source_parameters']:
            pd.read_csv(t/f'{name}.csv').to_excel(writer,sheet_name=name[:31],index=False)
    print(a[['dataset','method','n','rmse','wis','coverage95']].to_string(index=False))

if __name__=='__main__': main()
