"""Build review sheets from already-rendered PDF page images (no model imports)."""
from pathlib import Path
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]
for stem in ['main_en','main_cn','supplement_en','supplement_cn']:
    pages=sorted((ROOT/'audit/pages').glob(stem+'-*.png'))
    for start in range(0,len(pages),4):
        canvas=Image.new('RGB',(1100,1620),'#d0d0d0'); pen=ImageDraw.Draw(canvas)
        for j,path in enumerate(pages[start:start+4]):
            im=Image.open(path).convert('RGB'); im.thumbnail((530,770))
            x=(j%2)*550+10; y=(j//2)*810+30
            canvas.paste(im,(x,y)); pen.text((x,y-20),path.name,fill='black')
        canvas.save(ROOT/'audit/pages'/f'contact_{stem}_{start//4+1}.png')
