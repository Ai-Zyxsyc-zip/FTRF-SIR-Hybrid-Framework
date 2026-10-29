args=commandArgs(trailingOnly=FALSE);script=sub('^--file=','',args[grep('^--file=',args)][1]);ROOT=normalizePath(file.path(dirname(script),'..'),winslash='/')
options(scientific.definitions.only=TRUE);source(file.path(ROOT,'R_scientific_revision/forecast.R'))
all=list()
for(rep in 1:20)for(condition in c('matched','shifted'))for(k in 1:3) {
 b=if(k==1||condition=='matched')1.15+.02*(rep%%5) else 1.65;r=if(k==1||condition=='matched')1 else .4
 y=simulate_series(20260921+rep*10+k,b,r)
 all[[length(all)+1]]=data.frame(replicate=rep,condition=condition,region=c('target','source1','source2')[k],time=1:160,latent_percent=attr(y,'latent'),expected_reported_percent=r*attr(y,'latent'),observed=as.numeric(y),beta=b,gamma=.7,reporting_multiplier=r,noise_sd=.2)
}
d=do.call(rbind,all);savecsv(d,'simulation_truth')
p=read.csv(file.path(OUT,'simulation_predictions.csv'));check=merge(p,d[d$region=='target',],by.x=c('replicate','condition','test_time'),by.y=c('replicate','condition','time'))
stopifnot(nrow(check)==nrow(p),max(abs(check$observed.x-check$observed.y))<1e-10)
cat('Simulation truth archived and every saved test observation matched.\n')
