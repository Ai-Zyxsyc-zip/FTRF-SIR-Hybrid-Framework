"""Synchronize plotted figures without replacing manuscript text."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parents[2]/'latex'
def sync():
    for source in (ROOT/'English/figures').glob('*.pdf'):
        shutil.copy2(source,ROOT/'Chinese/figures'/source.name)
if __name__=='__main__': sync()
