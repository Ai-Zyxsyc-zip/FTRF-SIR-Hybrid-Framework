"""Scientific checks: solver accuracy/invariants, scores, time keys and saved results."""
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from mechanism import sir,fit_forecast
from scoring import score_one
from prepare_data import ROOT

def main():
    checks={}
    z=sir(1.05,.65,.97,.005,24)
    assert np.min(z)>=-1e-12 and np.max(abs(z.sum(axis=1)-1))<1e-12
    checks['mass_and_positivity']=True
    def rhs(t,y):
        return [-1.05*y[0]*y[1],1.05*y[0]*y[1]-.65*y[1],.65*y[1]]
    ref=solve_ivp(rhs,[0,24],[.97,.005,.025],t_eval=np.arange(25),rtol=1e-10,atol=1e-12).y.T
    errors=[float(np.max(abs(sir(1.05,.65,.97,.005,24,k)-ref))) for k in [8,16,32]]
    assert errors[2]<errors[1]<errors[0]
    checks['solver_max_abs_error_8_16_32']=errors
    # The simulation uses this discrete model rather than the exact ODE flow.
    pp,_=score_one(3.,3.,.5); bad,_=score_one(3.,6.,.5)
    assert pp['wis']<bad['wis'] and pp['nll']<bad['nll']
    checks['score_sanity']=True
    tab=ROOT/'results/tables'
    if (tab/'predictions.csv').exists():
        p=pd.read_csv(tab/'predictions.csv'); q=pd.read_csv(tab/'predictive_quantiles.csv')
        keys=['dataset','region','origin','horizon','method']
        assert not p.duplicated(keys).any()
        assert p.groupby(keys[:-1]).method.nunique().eq(9).all()
        assert np.isfinite(p[['se','ae','nll','wis','crps','scale']]).all().all()
        assert p.groupby(['dataset','method']).size().groupby(level=0).nunique().eq(1).all()
        obs=pd.read_csv(tab/'observations.csv')
        completed=set(map(tuple,p[['dataset','region','origin','horizon']].drop_duplicates().to_numpy()))
        excluded=set()
        for file in (ROOT/'results/folds').glob('*_excluded.json'):
            r=json.loads(file.read_text()); excluded.add((r['dataset'],r['region'],r['origin'],r['horizon']))
        planned=set()
        for dataset,d in obs.groupby('dataset'):
            seasons=[2013,2014] if dataset=='CPT' else [2018,2022]
            origins=d[(d.region=='Alabama')&d.season.isin(seasons)&d.week.isin([48,52,4,8,12,16])].t
            for origin in origins:
                for region in ['Alabama','California','Florida','Texas']:
                    for h in [1,2,4]: planned.add((dataset,region,int(origin),h))
        assert not completed.intersection(excluded)
        assert completed.union(excluded)==planned
        checks['planned_tasks']=len(planned); checks['excluded_tasks']=len(excluded)
        checks['every_planned_task_accounted_for']=True
        for _,d in q.groupby(keys): assert (np.diff(d.sort_values('probability').value)>=0).all()
        for _,r in p.sample(min(30,len(p)),random_state=7).iterrows():
            scores,_=score_one(r.y,r['loc'],r.scale)
            for metric in ['se','ae','nll','wis','crps']: assert np.isclose(scores[metric],r[metric],rtol=1e-10)
        checks['forecast_records']=len(p); checks['unique_complete_and_scores_recomputed']=True
    if (tab/'posterior_metrics.csv').exists():
        p=pd.read_csv(tab/'posterior_metrics.csv')
        checks['mh_max_conventional_split_rhat']=float(p.rhat.max())
        checks['mh_min_initial_positive_ess']=float(p.ess.min())
        checks['max_grid_boundary_mass']=float(p.edge_mass.max())
        assert p.rhat.max()<1.1
        assert p.edge_mass.max()<.01
    (ROOT/'audit/validation.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))

if __name__=='__main__': main()
