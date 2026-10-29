"""Collect existing, verified reporting assets without fitting any models."""
from pathlib import Path
import shutil,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'advisor_report_20260917'; OUT.mkdir(exist_ok=True)
(OUT/'figures').mkdir(exist_ok=True); (OUT/'data').mkdir(exist_ok=True)
names=['ext_ablation','ext_learning','ext_intervals','ext_sources','ext_horizons','posterior_fidelity','posterior_marginals']
manifest=[]
for name in names:
    for ext in ['png','pdf']:
        src=ROOT/'latex/figures'/f'{name}.{ext}'
        dest=OUT/'figures'/src.name; shutil.copy2(src,dest)
        manifest.append({'source':str(src),'copy':str(dest)})
for name in ['summary','paired','correction_decomposition','leave_group_out','interaction','horizons']:
    shutil.copy2(ROOT/'results/extension'/f'{name}.csv',OUT/'data'/f'{name}.csv')
p=pd.read_csv(ROOT/'results/extension/summary.csv')
d=p[(p.budget==52)&(p.interval=='Gaussian')&p.method.isin(['RF','FTRF-SIR','SIR_residual','Transfer_persistence'])]
print(d[['dataset','method','rmse','wis']].to_string(index=False))
(OUT/'asset_sources.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
print('Report assets:',OUT)
