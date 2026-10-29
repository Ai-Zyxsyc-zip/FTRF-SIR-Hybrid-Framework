"""Chronological multicentre forecast benchmark. Writes checkpoint tables per fold."""
from pathlib import Path
import argparse,json,platform,time
import numpy as np
import pandas as pd
from prepare_data import ROOT,STATES
from forecast_models import fit_all,METHODS
from scoring import score_one

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--smoke',action='store_true'); ap.add_argument('--dataset',choices=['ILI','CLS','CPT']); ap.add_argument('--collect',action='store_true'); ap.add_argument('--resume',action='store_true'); args=ap.parse_args()
    tab=ROOT/'results/tables'; features=pd.read_csv(tab/'causal_features.csv'); obs=pd.read_csv(tab/'observations.csv')
    out=ROOT/'results'/('smoke' if args.smoke else 'folds'); out.mkdir(exist_ok=True)
    for dataset in ([] if args.collect else ([args.dataset] if args.dataset else ['ILI','CLS','CPT'])):
        f=features[features.dataset==dataset]; o=obs[obs.dataset==dataset]
        seasons=[2013,2014] if dataset=='CPT' else [2018,2022]
        origins=o[(o.region==STATES[0]) & o.season.isin(seasons) & o.week.isin([48,52,4,8,12,16])].t.unique()
        for origin in origins[:1] if args.smoke else origins:
            for target in STATES[:1] if args.smoke else STATES[:4]:
                for h in [1] if args.smoke else [1,2,4]:
                    stem=f'{dataset}_{target.replace(" ","")}_{origin}_{h}'
                    if args.resume and (out/f'{stem}.csv').exists(): continue
                    start=time.perf_counter()
                    truth=o[['region','t','y']].rename(columns={'y':'label'}).copy(); truth.t-=h
                    dh=f.merge(truth,on=['region','t'])
                    seasonal=o[['region','t','y']].rename(columns={'y':'seasonal_h'}).copy(); seasonal.t+=52-h
                    dh=dh.merge(seasonal,on=['region','t'],how='left')
                    cutoff=int(origin)-20
                    cal=dh[(dh.region==target)&(dh.t>cutoff)&(dh.t+h<=origin)].dropna(subset=['label','seasonal_h'])
                    query=dh[(dh.region==target)&(dh.t==origin)].dropna(subset=['label','seasonal_h'])
                    if len(cal)<10 or len(query)!=1:
                        (out/f'{stem}_excluded.json').write_text(json.dumps({'dataset':dataset,'region':target,'origin':int(origin),'horizon':h,'n_cal':len(cal),'n_query':len(query),'reason':'Insufficient observed calibration labels or observed seasonal comparator / test truth. All methods excluded.'},indent=2))
                        print('excluded',stem,flush=True)
                        continue
                    tr=dh[(dh.t+h<=cutoff)&(dh.t>=cutoff-156)]
                    cp,cc,cs,_=fit_all(tr,target,cal,h,20260907)
                    # Gate uses only earlier half of the retrospective calibration block.
                    split=len(cal)//2; gate_y=cal.label.to_numpy()[:split]
                    weights=np.array([0,.25,.5,.75,1.])
                    w=weights[np.argmin([np.mean((v*cp['FTRF-SIR'][:split]+(1-v)*cp['RF_target'][:split]-gate_y)**2) for v in weights])]
                    cp['Guarded']=w*cp['FTRF-SIR']+(1-w)*cp['RF_target']
                    final=dh[(dh.t+h<=origin)&(dh.t>=origin-156)]
                    fp,fc,fs,n=fit_all(final,target,query,h,20260907)
                    fp['Guarded']=w*fp['FTRF-SIR']+(1-w)*fp['RF_target']
                    rows=[]; quant=[]
                    for method in METHODS:
                        # Same last calibration half for every method, including gate.
                        residual=cal.label.to_numpy()[split:]-cp[method][split:]
                        sigma=max(.1,float(np.sqrt(np.mean(residual**2))))
                        y=float(query.label.iloc[0]); loc=float(fp[method][0])
                        metrics,q=score_one(y,loc,sigma)
                        src_key='FTRF-SIR' if method=='Guarded' else method
                        source_cost=cs.get(src_key,0)+fs.get(src_key,0)
                        local=cc.get(method,0)+fc.get(method,0)
                        if method=='Guarded': local=cc['RF_target']+fc['RF_target']+cc['FTRF-SIR']+fc['FTRF-SIR']
                        mechanism_cost=0.
                        if method=='SIR':
                            mechanism_cost=pd.concat([cal,query]).drop_duplicates(['region','t']).mechanism_seconds.sum()
                        elif method in ['SIR_RF','FTRF-SIR','Guarded']:
                            relevant=pd.concat([tr,final,cal,query]).drop_duplicates(['region','t'])
                            if method=='SIR_RF': relevant=relevant[relevant.region==target]
                            mechanism_cost=relevant.mechanism_seconds.sum()
                        row=dict(dataset=dataset,region=target,origin=int(origin),season=int(query.season.iloc[0]),
                                 origin_year=int(query.year.iloc[0]),origin_week=int(query.week.iloc[0]),horizon=h,
                                 method=method,y=y,loc=loc,scale=sigma,gate_weight=w,n_target=n,
                                 n_calibration=len(residual),source_seconds=source_cost,local_seconds=local,
                                 fit_seconds=source_cost+local,mechanism_seconds=mechanism_cost,
                                 accounted_seconds=source_cost+local+mechanism_cost,**metrics)
                        rows.append(row)
                        quant.extend(dict(dataset=dataset,region=target,origin=int(origin),horizon=h,method=method,probability=float(p),value=float(v)) for p,v in q.items())
                    pd.DataFrame(rows).to_csv(out/f'{stem}.csv',index=False)
                    pd.DataFrame(quant).to_csv(out/f'{stem}_quantiles.csv',index=False)
                    print(stem,'done',round(time.perf_counter()-start,2),'s',flush=True)
        if args.smoke: break
    if not args.smoke and not args.dataset:
        pd.concat([pd.read_csv(p) for p in sorted(out.glob('*.csv')) if not p.stem.endswith('_quantiles')]).to_csv(tab/'predictions.csv',index=False)
        pd.concat([pd.read_csv(p) for p in sorted(out.glob('*_quantiles.csv'))]).to_csv(tab/'predictive_quantiles.csv',index=False)
        (ROOT/'audit/forecast_run.json').write_text(json.dumps({'seed':20260907,'python':platform.python_version(),
          'target_window':52,'source_window':156,'calibration_block':20,'trees':80,'states':STATES,
          'timing':'Local and source forest fit+query prediction; historical SIR feature construction separately tabulated. No networking included.'},indent=2))

if __name__=='__main__': main()
