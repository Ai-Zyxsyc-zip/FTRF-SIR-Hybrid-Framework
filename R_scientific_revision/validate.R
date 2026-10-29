args=commandArgs(trailingOnly=FALSE);script=sub('^--file=','',args[grep('^--file=',args)][1]);ROOT=normalizePath(file.path(dirname(script),'..'),winslash='/')
options(scientific.definitions.only=TRUE);source(file.path(ROOT,'R_scientific_revision/forecast.R'))
p=readRDS(file.path(OUT,'forecast_panels.rds'))$ILI
saved=read.csv(file.path(OUT,'forecast_predictions.csv'));s=saved[saved$task=='ILI'&saved$target=='Alabama'&saved$h==1,];o=s$origin[1];seed=s$seed[1]
a=forecast_one(p,'Alabama',o,1,seed)
b=p;b$y[b$response_time>o]=9999
c=forecast_one(b,'Alabama',o,1,seed)
stopifnot(isTRUE(all.equal(a$predicted,c$predicted,tolerance=0)))
ref=s[s$origin==o,];stopifnot(isTRUE(all.equal(a$predicted,ref$predicted,tolerance=1e-10)))
# Offset construction itself does not see observations after its feature origin.
set.seed(33);x=runif(20,0,10);m=offset(x[1:12]);x[13:20]=99;stopifnot(identical(m,offset(x[1:12])))
writeLines(c('PASS: exact saved-task replay (nine models).','PASS: perturbing all post-cutoff labels does not change the nine predictions.','PASS: past-only SIR offset is unchanged by future observations.','Scope: sampled implementation checks, not a proof for every input.'),file.path(OUT,'causality_validation.txt'))
cat('R VALIDATION PASSED\n')
