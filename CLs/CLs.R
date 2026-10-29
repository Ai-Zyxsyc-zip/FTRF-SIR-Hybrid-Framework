# Run the shared analysis for this dataset using the documented protocol.
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
options(scientific.root=project_root)
options(scientific.helpers.only=FALSE, scientific.definitions.only=FALSE)
source(file.path(project_root,'R_scientific_revision/run.R'))
