"""Archive obsolete generated assets and remove only known build intermediates."""
from pathlib import Path
import json,shutil
ROOT=Path(__file__).resolve().parents[1]
archive=ROOT/'audit/legacy_assets'; archive.mkdir(exist_ok=True)
figures={'ext_ablation','ext_horizons','ext_intervals','ext_learning','ext_sources','posterior_fidelity','posterior_marginals'}
tables={'full_metrics.tex','results_cn.tex','results_en.tex'}
moved=[]
for folder in ['figures','tables']:
    for p in (ROOT/'latex'/folder).iterdir():
        obsolete=p.stem not in figures if folder=='figures' else p.name in tables
        if p.is_file() and obsolete:
            assert p.resolve().is_relative_to(ROOT.resolve())
            dest=archive/folder/p.name; dest.parent.mkdir(exist_ok=True)
            if dest.exists(): assert dest.read_bytes()==p.read_bytes(); p.unlink()
            else: shutil.move(str(p),str(dest))
            moved.append(str(p.relative_to(ROOT)))
old=ROOT/'audit/DELIVERY_AUDIT.md'
if not (ROOT/'audit/DELIVERY_AUDIT_PRELIMINARY.md').exists():
    shutil.copy2(old,ROOT/'audit/DELIVERY_AUDIT_PRELIMINARY.md')
shutil.copy2(ROOT/'audit/EXTENSION_DELIVERY.md',old)
logs=ROOT/'audit/compilation'; logs.mkdir(exist_ok=True)
for p in (ROOT/'latex').iterdir():
    if p.is_file() and (p.suffix in {'.aux','.blg','.log','.out','.spl'} or p.name.endswith('.build.txt')):
        assert p.resolve().is_relative_to((ROOT/'latex').resolve())
        if p.suffix=='.log': shutil.copy2(p,logs/(p.stem+'.txt'))
        p.unlink()
helper=ROOT/'audit/refresh_readme.py'
if helper.exists(): helper.unlink()
(ROOT/'audit/cleanup.json').write_text(json.dumps({'archived_obsolete_latex_assets':moved},indent=2),encoding='utf8')
print('Archived obsolete LaTeX assets:',len(moved))
