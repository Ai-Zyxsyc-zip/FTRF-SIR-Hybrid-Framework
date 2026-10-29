"""Budget-matched retrospective experiment. See audit/EXTENSION_PROTOCOL.md.

Run --workers 4; resume requires identical code/input fingerprint. No networking.
Each job shares source fits across targets/budgets, but never uses future labels.
"""
from pathlib import Path
import argparse, hashlib, json, os, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.ensemble import RandomForestRegressor
from prepare_data import ROOT, FEATURES, STATES

OUT=ROOT/'results/extension'
SEED=20260916
BUDGETS=[24,52,104]
ALPHA=np.array([.05,.1,.2,.3,.4,.5,.6,.7,.8,.9])
PROBS=np.unique(np.r_[.5,ALPHA/2,1-ALPHA/2])
METHODS=['Persistence','Seasonal','SIR','RF','Pool','Transfer','SIR_feature',
         'Transfer_feature','SIR_residual','FTRF-SIR','Persistence_residual',
         'Transfer_persistence','Half_SIR_RF']

def rf(seed=SEED):
    return RandomForestRegressor(n_estimators=80,min_samples_leaf=3,
                                 max_features=1.,random_state=seed,n_jobs=1)

def panel(dataset,h):
    f=pd.read_csv(ROOT/'results/tables/causal_features.csv')
    o=pd.read_csv(ROOT/'results/tables/observations.csv')
    f=f[f.dataset==dataset].copy(); o=o[o.dataset==dataset].copy()
    y=o[['region','t','y']].rename(columns={'y':'label'}); y.t-=h
    s=o[['region','t','y']].rename(columns={'y':'seasonal_h'}); s.t+=52-h
    return f.merge(y,on=['region','t'],how='left').merge(s,on=['region','t'],how='left')

def history(d,cutoff,h,n):
    return d[(d.t+h<=cutoff)&d.label.notna()].sort_values('t').tail(n)

def source_bank(d,cutoff,h):
    bank={}; durations={}
    for region in STATES:
        tr=history(d[d.region==region],cutoff,h,156)
        if len(tr)<24: continue
        x=tr[FEATURES].to_numpy(); y=tr.label.to_numpy()
        bank[region]={}; durations[region]={}
        for key,offset in [('raw',0),('sir',tr[f'sir{h}'].to_numpy()),('persist',tr.lag0.to_numpy())]:
            start=time.perf_counter()
            bank[region][key]=rf(SEED+100+STATES.index(region)).fit(x,y-offset)
            durations[region][key]=time.perf_counter()-start
    return bank,durations

def source_order(d,target,cutoff):
    target_rows=d[(d.region==target)&(d.t<=cutoff)].sort_values('t').tail(26)
    scores=[]
    for source in STATES:
        if source==target: continue
        z=target_rows[['t','lag0']].merge(d[d.region==source][['t','lag0']],on='t')
        distance=float(np.sqrt(np.mean((z.lag0_x-z.lag0_y)**2))/max(np.std(z.lag0_x),.5)) if len(z)>=20 else float('inf')
        scores.append((source,distance,len(z)))
    return sorted(scores,key=lambda x:(x[1],x[0]))

def fit_predict(d,target,query,cutoff,h,budget,bank,bankcost):
    tr=history(d[d.region==target],cutoff,h,budget)
    assert len(tr)==budget and (tr.t+h<=cutoff).all()
    x=tr[FEATURES].to_numpy(); q=query[FEATURES].to_numpy(); y=tr.label.to_numpy()
    m=tr[f'sir{h}'].to_numpy(); qm=query[f'sir{h}'].to_numpy()
    p=tr.lag0.to_numpy(); qp=query.lag0.to_numpy()
    sources=[s for s in STATES if s!=target]
    aug={}; costs={}; used={}
    for key in ['raw','sir','persist']:
        a=[]; b=[]; querytime=0.
        for s in sources:
            tic=time.perf_counter(); a.append(bank[s][key].predict(x)); b.append(bank[s][key].predict(q))
            querytime+=time.perf_counter()-tic
        aug[key]=(np.column_stack(a),np.column_stack(b))
        used[key]=sum(bankcost[s][key] for s in sources)+querytime
    raw,rq=aug['raw']; res,sq=aug['sir']; prs,pq=aug['persist']
    pool=pd.concat([tr]+[history(d[d.region==s],cutoff,h,156) for s in sources])
    tasks=[('RF',x,y,q,0,None),('Pool',pool[FEATURES].to_numpy(),pool.label.to_numpy(),q,0,None),
      ('Transfer',np.c_[x,raw],y,np.c_[q,rq],0,'raw'),
      ('SIR_feature',np.c_[x,m],y,np.c_[q,qm],0,None),
      ('Transfer_feature',np.c_[x,m,raw],y,np.c_[q,qm,rq],0,'raw'),
      ('SIR_residual',np.c_[x,m],y-m,np.c_[q,qm],qm,None),
      ('FTRF-SIR',np.c_[x,m,res],y-m,np.c_[q,qm,sq],qm,'sir'),
      ('Persistence_residual',np.c_[x,p],y-p,np.c_[q,qp],qp,None),
      ('Transfer_persistence',np.c_[x,p,prs],y-p,np.c_[q,qp,pq],qp,'persist')]
    ranks=source_order(d,target,cutoff)
    if budget==52:
        for name,selected in [('FTRF-SIR_near',[r[0] for r in ranks[:2]]),('FTRF-SIR_far',[r[0] for r in ranks[-2:]])]:
            ix=[sources.index(s) for s in selected]
            tasks.append((name,np.c_[x,m,res[:,ix]],y-m,np.c_[q,qm,sq[:,ix]],qm,'sir'))
    preds={'Persistence':qp,'Seasonal':query.seasonal_h.to_numpy(),'SIR':qm}
    for name,xt,yt,xq,offset,key in tasks:
        tic=time.perf_counter(); preds[name]=np.clip(rf().fit(xt,yt).predict(xq)+offset,0,100)
        costs[name]=(time.perf_counter()-tic,used[key] if key else 0.)
    preds['Half_SIR_RF']=(preds['RF']+qm)/2
    costs['Half_SIR_RF']=costs['RF']
    for name in ['Persistence','Seasonal','SIR']: costs[name]=(0.,0.)
    return preds,costs,ranks

def quantiles(loc,errors,kind):
    if kind=='Gaussian':
        sigma=max(.1,float(np.sqrt(np.mean(errors**2))))
        return np.clip(loc+sigma*norm.ppf(PROBS),0,100)
    e=np.sort(np.abs(errors)); n=len(e); values=[]
    for prob in PROBS:
        if prob==.5: values.append(loc); continue
        alpha=2*min(prob,1-prob)
        k=int(np.ceil((n+1)*(1-alpha)-1e-10))
        radius=e[k-1] if k<=n else 100.
        values.append(np.clip(loc+np.sign(prob-.5)*radius,0,100))
    return np.asarray(values)

def metrics(y,loc,q):
    lower=np.interp(ALPHA/2,PROBS,q); upper=np.interp(1-ALPHA/2,PROBS,q)
    iscore=upper-lower+2/ALPHA*np.maximum(lower-y,0)+2/ALPHA*np.maximum(y-upper,0)
    wis=(.5*abs(y-np.interp(.5,PROBS,q))+np.sum(ALPHA/2*iscore))/(len(ALPHA)+.5)
    lo,hi=np.interp([.025,.975],PROBS,q)
    # Quantile forecasts do not define a unique density, so NLL is not reported.
    return dict(se=(loc-y)**2,ae=abs(loc-y),wis=wis,covered95=int(lo<=y<=hi),width95=hi-lo,lo95=lo,hi95=hi)

def job(spec):
    dataset,origin,h=spec; start=time.perf_counter(); d=panel(dataset,h)
    cutoff=origin-36; banks=[]; bc=[]
    for c in [cutoff,origin]:
        b,cost=source_bank(d,c,h); banks.append(b); bc.append(cost)
    records=[]; qs=[]; residuals=[]; exclusions=[]; ranking=[]
    for target in STATES:
        t=d[d.region==target]
        cal=t[(t.t>cutoff)&(t.t+h<=origin)].dropna(subset=['label','seasonal_h']).sort_values('t')
        query=t[t.t==origin].dropna(subset=['label','seasonal_h'])
        eligible=(len(cal)>=26 and len(query)==1 and all(len(b)==6 for b in banks)
                  and all(len(history(t,c,h,104))==104 for c in [cutoff,origin]))
        if not eligible:
            exclusions.append(dict(dataset=dataset,origin=origin,horizon=h,region=target,n_cal=len(cal),
                n_history=len(history(t,cutoff,h,104)),reason='Common calibration/history/source/query eligibility'))
            continue
        for n in BUDGETS:
            cp,cc,_=fit_predict(d,target,cal,cutoff,h,n,banks[0],bc[0])
            fp,fc,ranks=fit_predict(d,target,query,origin,h,n,banks[1],bc[1])
            common=dict(dataset=dataset,region=target,origin=origin,season=int(query.season.iloc[0]),
                        year=int(query.year.iloc[0]),week=int(query.week.iloc[0]),horizon=h,budget=n)
            if n==52:
                ranking.extend(dict(**common,source=s,distance=dist,rank=idx+1,n_overlap=nn) for idx,(s,dist,nn) in enumerate(ranks))
            for method in fp:
                errors=cal.label.to_numpy()-cp[method]
                loc=float(fp[method][0]); y=float(query.label.iloc[0])
                residuals.extend(dict(**common,method=method,cal_origin=int(ct),error=float(e)) for ct,e in zip(cal.t,errors))
                for kind in ['Gaussian','Order']:
                    q=quantiles(loc,errors,kind)
                    records.append(dict(**common,method=method,interval=kind,y=y,loc=loc,n_cal=len(errors),
                        local_seconds=cc[method][0]+fc[method][0],source_seconds=cc[method][1]+fc[method][1],
                        **metrics(y,loc,q)))
                    qs.extend(dict(**common,method=method,interval=kind,probability=float(p),value=float(v)) for p,v in zip(PROBS,q))
    stem=f'{dataset}_{origin}_{h}'
    for suffix,rows in [('predictions',records),('quantiles',qs),('errors',residuals),('excluded',exclusions),('sources',ranking)]:
        pd.DataFrame(rows).to_csv(OUT/'folds'/f'{stem}_{suffix}.csv',index=False)
    (OUT/'folds'/f'{stem}.done').write_text(json.dumps({'seconds':time.perf_counter()-start}))
    return stem,len(records),round(time.perf_counter()-start,1)

def collect():
    for name in ['predictions','quantiles','errors','excluded','sources']:
        frames=[]
        for p in sorted((OUT/'folds').glob(f'*_{name}.csv')):
            try: frames.append(pd.read_csv(p))
            except pd.errors.EmptyDataError: pass
        pd.concat(frames,ignore_index=True).to_csv(OUT/f'{name}.csv',index=False) if frames else (OUT/f'{name}.csv').write_text('')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--smoke',action='store_true'); ap.add_argument('--collect',action='store_true')
    args=ap.parse_args(); (OUT/'folds').mkdir(parents=True,exist_ok=True)
    if args.collect: collect(); return
    paths=[Path(__file__),ROOT/'code/prepare_data.py',ROOT/'results/tables/causal_features.csv',ROOT/'results/tables/observations.csv']
    fp={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    manifest=OUT/'run_manifest.json'
    if manifest.exists(): assert json.loads(manifest.read_text())['hashes']==fp,'Inputs/code changed; use a new output directory'
    else: manifest.write_text(json.dumps({'hashes':fp,'seed':SEED,'budgets':BUDGETS,'workers':args.workers},indent=2))
    obs=pd.read_csv(ROOT/'results/tables/observations.csv'); specs=[]
    for dataset in ['ILI','CLS','CPT']:
        seasons=[2012,2013,2014] if dataset=='CPT' else [2018,2019,2022,2023]
        origins=obs[(obs.dataset==dataset)&(obs.region==STATES[0])&obs.season.isin(seasons)&obs.week.isin([48,4,12,16])].t.unique()
        specs.extend((dataset,int(t),h) for t in origins for h in [1,2,4])
    (OUT/'planned_jobs.json').write_text(json.dumps(specs))
    pending=[s for s in specs if not (OUT/'folds'/f'{s[0]}_{s[1]}_{s[2]}.done').exists()]
    if args.smoke: print(job(pending[0]),flush=True); return
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures=[ex.submit(job,s) for s in pending]
        for i,f in enumerate(as_completed(futures)): print(i+1,'/',len(futures),f.result(),flush=True)
    collect()

if __name__=='__main__': main()
