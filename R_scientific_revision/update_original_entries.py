from pathlib import Path
import zipfile,hashlib,json
root=Path(__file__).resolve().parents[1]
backup=root/'R_scientific_revision/original_R_sources_20260921.zip'
files=[p for folder in ['CLs','CPT','ILI','PHLs'] for p in (root/folder).iterdir() if p.suffix.lower()=='.r']
assert not backup.exists(),'Original backup already exists; do not overwrite.'
with zipfile.ZipFile(backup,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,p.relative_to(root))
manifest=[]
for p in files:
 manifest.append({'file':str(p.relative_to(root)),'original_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 mode='helpers' if p.stem.lower()=='funcs' else 'forecast' if any(x in p.stem.upper() for x in ['SIR','RFTAR']) else 'conditional'
 s='''# Scientifically revised entry point. Original source is preserved in
# R_scientific_revision/original_R_sources_20260921.zip.
# This shared runner evaluates all datasets with the same documented protocol.
locate_root <- function() {
  of <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  a <- commandArgs(trailingOnly=FALSE)
  f <- sub('^--file=', '', a[grep('^--file=', a)])
  candidates <- c(if(!is.null(of))dirname(of), if(length(f))dirname(f[1]), getwd())
  for(c in candidates) for(k in 0:3) {
    if(dir.exists(file.path(c,'R_scientific_revision'))) return(normalizePath(c,winslash='/'))
    c <- dirname(c)
  }
  stop('Run from this project or provide the script path to Rscript.')
}
project_root <- locate_root()
'''
 if mode=='helpers':s+="options(scientific.helpers.only=TRUE)\n"
 s+=f"source(file.path(project_root,'R_scientific_revision/{'forecast.R' if mode=='forecast' else 'run.R'}'))\n"
 p.write_text(s,encoding='utf8')
 manifest[-1]['revised_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'R_scientific_revision/source_revision_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(len(files),'original R files backed up and redirected to repaired shared implementations.')
