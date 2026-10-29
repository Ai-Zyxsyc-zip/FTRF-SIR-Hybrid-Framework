"""Rank-normalized and folded split R-hat from the archived four MH chains."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import norm,rankdata
ROOT=Path(__file__).resolve().parents[1]

def split_rhat(x):
    n=x.shape[1]//2; x=np.concatenate([x[:,:n],x[:,-n:]])
    w=np.var(x,axis=1,ddof=1).mean(); b=n*np.var(x.mean(axis=1),ddof=1)
    return np.sqrt(((n-1)/n*w+b/n)/w)

def ranked(x):
    z=norm.ppf((rankdata(x.ravel())-.375)/(x.size+.25)).reshape(x.shape)
    return split_rhat(z)

def main():
    d=pd.read_csv(ROOT/'results/tables/posterior_samples.csv'); d=d[d.method=='MH']; rows=[]
    for keys,g in d.groupby(['scenario','rep','prior']):
        r=dict(zip(['scenario','rep','prior'],keys))
        for param in ['beta','gamma']:
            x=np.stack([v.sort_values('draw')[param].to_numpy() for _,v in g.groupby('chain')])
            assert x.shape==(4,1600)
            r[param+'_rank_rhat']=ranked(x)
            r[param+'_folded_rhat']=ranked(abs(x-np.median(x)))
        r['max_rank_folded_rhat']=max(v for k,v in r.items() if k.endswith('_rhat'))
        rows.append(r)
    a=pd.DataFrame(rows); a.to_csv(ROOT/'results/tables/posterior_rank_diagnostics.csv',index=False)
    value=float(a.max_rank_folded_rhat.max())
    (ROOT/'results/tables/rank_rhat.tex').write_text(f'\\newcommand{{\\MaxRankRhat}}{{{value:.3f}}}\n')
    # Independent sanity cases: matched stationary chains and separated locations.
    rng=np.random.default_rng(77); x=rng.normal(size=(4,2000))
    assert ranked(x)<1.02 and ranked(x+np.arange(4)[:,None]*3)>1.2
    print('Targets:',len(a),'max rank/folded split Rhat:',value)

if __name__=='__main__': main()
