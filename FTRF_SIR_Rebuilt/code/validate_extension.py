"""Independent score checks, common-case reconciliation, and future-data test."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/extension'
KEY=['dataset','region','season','origin','horizon','budget','method','interval']

def main():
    p=pd.read_csv(OUT/'predictions.csv'); q=pd.read_csv(OUT/'quantiles.csv')
    assert not p.duplicated(KEY).any() and not q.duplicated(KEY+['probability']).any()
    v=q.pivot(index=KEY,columns='probability',values='value').sort_index()
    pp=p.set_index(KEY).loc[v.index]; yy=pp.y.to_numpy()[:,None]; probs=v.columns.to_numpy()
    assert v.shape[1]==21 and np.isfinite(v.to_numpy()).all()
    assert np.min(np.diff(v.to_numpy(),axis=1))>=-1e-10
    assert v.min().min()>=0 and v.max().max()<=100
    residual=yy-v.to_numpy()
    wis=2*np.mean(np.maximum(probs*residual,(probs-1)*residual),axis=1)
    wis_error=float(np.max(abs(wis-pp.wis.to_numpy())))
    assert wis_error<1e-8
    assert np.max(abs(p.se-(p.y-p['loc'])**2))<1e-8
    assert np.max(abs(p.ae-abs(p.y-p['loc'])))<1e-8
    assert (p.n_cal>=26).all()
    for keys,g in p.groupby(['dataset','region','origin','horizon']):
        assert set(g.budget)=={24,52,104}
        assert g.groupby('budget').size().to_dict()=={24:26,52:30,104:26}
        assert g.y.nunique()==1
    jobs=json.loads((OUT/'planned_jobs.json').read_text())
    assert len(list((OUT/'folds').glob('*.done')))==len(jobs)
    exclusions=pd.read_csv(OUT/'excluded.csv')
    tasks=p[['dataset','region','origin','horizon']].drop_duplicates()
    assert len(tasks)+len(exclusions)==len(jobs)*6
    assert tasks.merge(exclusions,on=['dataset','region','origin','horizon']).empty
    point=p[p.interval=='Gaussian']; other=p[p.interval=='Order']
    assert np.max(abs(point.set_index(KEY[:-1])['loc']-other.set_index(KEY[:-1])['loc']))==0
    # Independently interpolate interval endpoints from the saved quantiles.
    i0=np.argmin(abs(probs-.025)); i1=np.argmin(abs(probs-.975))
    cover=((yy[:,0]>=v.iloc[:,i0])&(yy[:,0]<=v.iloc[:,i1])).astype(int)
    assert np.array_equal(cover,pp.covered95)
    # Perturb information strictly unavailable at the final cutoff.
    from extension import panel,source_bank,fit_predict
    row=point[(point.budget==52)&(point.method=='FTRF-SIR')].iloc[0]
    ds,origin,h,target=row.dataset,int(row.origin),int(row.horizon),row.region
    d=panel(ds,h); query=d[(d.region==target)&(d.t==origin)].copy()
    b,bc=source_bank(d,origin,h); expected,_,_=fit_predict(d,target,query,origin,h,52,b,bc)
    changed=d.copy(); changed.loc[changed.t+h>origin,'label']=99.123
    from prepare_data import FEATURES
    changed.loc[changed.t>origin,FEATURES]=77.777
    b2,bc2=source_bank(changed,origin,h); result,_,_=fit_predict(changed,target,query,origin,h,52,b2,bc2)
    future_error=max(float(np.max(abs(result[m]-expected[m]))) for m in expected)
    assert future_error<1e-12
    saved=point[(point.dataset==ds)&(point.origin==origin)&(point.horizon==h)&(point.region==target)&(point.budget==52)].set_index('method')
    replay_error=max(abs(float(expected[m][0])-saved.loc[m,'loc']) for m in expected)
    assert replay_error<1e-9
    report=dict(eligible_tasks=len(tasks),excluded_tasks=len(exclusions),planned_tasks=len(jobs)*6,
                forecast_rows=len(p),quantile_rows=len(q),max_wis_identity_error=wis_error,
                future_perturbation_max_difference=future_error,replayed_task_max_difference=replay_error,
                replay_scope=f'{ds}/{target}/{origin}/h={h}/budget=52; all 15 models',
                raw_scores_recomputed=True,common_tasks_verified=True)
    (ROOT/'audit/extension_validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
