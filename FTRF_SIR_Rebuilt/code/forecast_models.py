"""Nine baselines/ablations; source models only see labels before the cutoff."""
import time
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from prepare_data import FEATURES

METHODS=['Persistence','Seasonal','SIR','RF_target','RF_pool','FTRF','SIR_RF','FTRF-SIR','Guarded']

def forest(seed):
    return RandomForestRegressor(n_estimators=80,min_samples_leaf=3,max_features=1.,random_state=seed,n_jobs=1)

def fit_all(train,target,query,h,seed):
    """Query covariates are causal; train is already filtered by label availability."""
    train=train.dropna(subset=FEATURES+['label',f'sir{h}'])
    tar=train[train.region==target].sort_values('t').tail(52)
    if len(tar)<20: raise ValueError('Insufficient target training history')
    x=tar[FEATURES].to_numpy(); y=tar.label.to_numpy(); q=query[FEATURES].to_numpy()
    base=tar[f'sir{h}'].to_numpy(); qb=query[f'sir{h}'].to_numpy()
    pred={'Persistence':query.lag0.to_numpy(),'Seasonal':query.seasonal_h.to_numpy(),'SIR':qb}
    costs={m:0. for m in METHODS}; srcpred=[]; srcquery=[]; srcr=[]; srcqr=[]
    source_seconds={'FTRF':0.,'FTRF-SIR':0.}
    for idx,(region,r) in enumerate(train[train.region!=target].groupby('region')):
        r=r.sort_values('t').tail(156)
        if len(r)<20: continue
        start=time.perf_counter()
        f=forest(seed+idx+100).fit(r[FEATURES].to_numpy(),r.label.to_numpy())
        srcpred.append(f.predict(x)); srcquery.append(f.predict(q))
        source_seconds['FTRF']+=time.perf_counter()-start
        start=time.perf_counter()
        f=forest(seed+idx+100).fit(r[FEATURES].to_numpy(),r.label.to_numpy()-r[f'sir{h}'].to_numpy())
        srcr.append(f.predict(x)); srcqr.append(f.predict(q))
        source_seconds['FTRF-SIR']+=time.perf_counter()-start
    if not srcpred: raise ValueError('No usable source models')
    xx=np.column_stack([x]+srcpred); qq=np.column_stack([q]+srcquery)
    xr=np.column_stack([x,base]+srcr); qr=np.column_stack([q,qb]+srcqr)
    tasks=[('RF_target',x,y,q,0),('RF_pool',train[FEATURES].to_numpy(),train.label.to_numpy(),q,0),
           ('FTRF',xx,y,qq,0),('SIR_RF',np.column_stack([x,base]),y-base,np.column_stack([q,qb]),qb),
           ('FTRF-SIR',xr,y-base,qr,qb)]
    for name,xt,yt,xq,offset in tasks:
        start=time.perf_counter()
        pred[name]=forest(seed).fit(xt,yt).predict(xq)+offset
        costs[name]=time.perf_counter()-start
    return {k:np.clip(v,0,100) for k,v in pred.items()},costs,source_seconds,len(tar)
