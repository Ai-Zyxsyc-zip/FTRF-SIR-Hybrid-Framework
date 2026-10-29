"""Rerun one actual forecast task and compare its saved numerical outputs."""
from pathlib import Path
import subprocess,sys,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,str(ROOT/'code/run_forecasts.py'),'--smoke'],cwd=ROOT,check=True)
files=[p for p in (ROOT/'results/smoke').glob('*.csv') if not p.stem.endswith('_quantiles')]
fresh=pd.read_csv(files[0]); full=pd.read_csv(ROOT/'results/tables/predictions.csv')
keys=['dataset','region','origin','horizon','method']
joined=fresh.merge(full,on=keys,suffixes=('_new','_saved'),validate='one_to_one')
assert len(joined)==9
metrics=['y','loc','scale','gate_weight','mean','median','se','ae','nll','wis','crps','lo95','hi95']
errors={m:float(np.max(abs(joined[m+'_new']-joined[m+'_saved']))) for m in metrics}
assert max(errors.values())<1e-10
report={'rerun':'ILI / Alabama / first selected origin / horizon 1','models':9,
        'maximum_absolute_differences':errors,'timing_excluded':True,'passed':True}
(ROOT/'audit/reproduction.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
