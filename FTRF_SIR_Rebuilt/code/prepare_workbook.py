"""Typed JSON transport for workbook creation; all scientific sources are CSV."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; T=ROOT/'results/extension'
blocks=[]
for name,file in [('Summary','summary'),('Horizons','horizons'),('Paired comparisons','paired'),
                  ('Correction','correction_decomposition'),('Group sensitivity','leave_group_out')]:
    d=pd.read_csv(T/f'{file}.csv')
    if file=='paired': d=d[(d.budget==52)&(d.interval=='Gaussian')&(d.method=='FTRF-SIR')]
    note='RMSE, MAE, WIS, width: percentage points. MSE: squared points. Coverage: fraction. Time: seconds.'
    if file in ['correction_decomposition','leave_group_out']: note='All error decomposition and MSE difference columns: squared percentage points.'
    if file=='paired': note='Paired cluster differences: MSE in squared points; WIS in points. Intervals are descriptive.'
    blocks.append(dict(name=name,title=name+': retrospective influenza forecasting',note=note,source=f'results/extension/{file}.csv',columns=list(d.columns),rows=json.loads(d.to_json(orient='values'))))
d=pd.read_csv(ROOT/'results/tables/posterior_summary.csv')
d=d[['scenario','prior','method','rmse','nll','coverage95','wasserstein_std','mmd2','seconds','ess','rhat']]
blocks.append(dict(name='Posterior',title='Independent SIR parameter-inference simulation',note='Ten independent datasets per source condition. Conventional ESS and split Rhat; time in seconds.',source='results/tables/posterior_summary.csv',columns=list(d.columns),rows=json.loads(d.to_json(orient='values'))))
(T/'workbook_data.json').write_text(json.dumps(blocks,ensure_ascii=False),encoding='utf-8')
