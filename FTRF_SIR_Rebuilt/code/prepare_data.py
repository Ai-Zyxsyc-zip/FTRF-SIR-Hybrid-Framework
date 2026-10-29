"""Audit raw snapshots, align CDC year/week keys, and create causal features."""
from pathlib import Path
import hashlib, json, time
import numpy as np
import pandas as pd
from mechanism import fit_forecast

ROOT=Path(__file__).resolve().parents[1]
STATES=['Alabama','California','Florida','Texas','Illinois','New York']
SPECS={'ILI':('ILINet.csv','%UNWEIGHTED ILI'),
       'CLS':('WHO_NREVSS_Clinical_Labs.csv','PERCENT POSITIVE'),
       'CPT':('WHO_NREVSS_Combined_prior_to_2015_16.csv','PERCENT POSITIVE')}
FEATURES=['lag0','lag1','lag2','lag3','mean4','slope','sin_week','cos_week']

def main():
    tables=ROOT/'results/tables'; audit=[]; manifests=[]; features=[]; panels=[]
    for file in sorted((ROOT/'data/raw').glob('*.csv')):
        manifests.append(dict(file=file.name,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),bytes=file.stat().st_size))
    for dataset,(filename,column) in SPECS.items():
        d=pd.read_csv(ROOT/'data/raw'/filename,skiprows=1,na_values=['X','NA',''])
        assert not d.duplicated(['REGION','YEAR','WEEK']).any()
        keys=d[['YEAR','WEEK']].drop_duplicates().sort_values(['YEAR','WEEK']).reset_index(drop=True)
        keys['t']=np.arange(len(keys))
        keys['season']=keys.YEAR-(keys.WEEK<40).astype(int)
        for region in STATES:
            r=keys.merge(d[d.REGION==region][['YEAR','WEEK',column]],how='left')
            y=pd.to_numeric(r[column],errors='coerce')
            y=y.where(y.between(0,100))
            audit.append(dict(dataset=dataset,region=region,rows=len(r),missing=int(y.isna().sum()),
                              first_year=int(r.YEAR.min()),last_year=int(r.YEAR.max()),outcome=column))
            # At most two missing past weeks carried forward; never fill truth.
            yf=y.ffill(limit=2).to_numpy(); truth=y.to_numpy()
            for j,row in r.iterrows():
                panels.append(dict(dataset=dataset,region=region,t=j,year=int(row.YEAR),week=int(row.WEEK),season=int(row.season),y=truth[j]))
                if j<11 or not np.isfinite(yf[j-11:j+1]).all(): continue
                started=time.perf_counter()
                pred,beta,s,success=fit_forecast(yf[j-11:j+1])
                cost=time.perf_counter()-started
                x={'lag0':yf[j],'lag1':yf[j-1],'lag2':yf[j-2],'lag3':yf[j-3],
                   'mean4':np.mean(yf[j-3:j+1]),'slope':yf[j]-yf[j-3],
                   'sin_week':np.sin(2*np.pi*row.WEEK/52.1775),
                   'cos_week':np.cos(2*np.pi*row.WEEK/52.1775)}
                features.append(dict(dataset=dataset,region=region,t=j,year=int(row.YEAR),week=int(row.WEEK),
                                     season=int(row.season),**x,**{f'sir{h}':pred[h-1] for h in [1,2,3,4]},
                                     beta=beta,susceptible=s,fit_success=success,mechanism_seconds=cost,
                                     seasonal=yf[j-52] if j>=52 else np.nan))
            print('prepared',dataset,region,flush=True)
    pd.DataFrame(panels).to_csv(tables/'observations.csv',index=False)
    pd.DataFrame(features).to_csv(tables/'causal_features.csv',index=False)
    pd.DataFrame(audit).to_csv(tables/'data_quality.csv',index=False)
    (ROOT/'audit/raw_manifest.json').write_text(json.dumps(manifests,indent=2),encoding='utf-8')
    # The short PHL snapshot is unsuitable for the same historical protocol.
    p=pd.read_csv(ROOT/'data/raw/WHO_NREVSS_Public_Health_Labs.csv',skiprows=1)
    (ROOT/'audit/phl_exclusion.json').write_text(json.dumps({'rows':len(p),'seasons':int(p.SEASON_DESCRIPTION.nunique()),'reason':'Season-level totals; no YEAR/WEEK fields. Not a weekly forecasting dataset.'},indent=2))

if __name__=='__main__': main()
