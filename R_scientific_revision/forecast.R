args=commandArgs(trailingOnly=FALSE);script=sub('^--file=','',args[grep('^--file=',args)][1]);ROOT=getOption('scientific.root',normalizePath(file.path(dirname(script),'..'),winslash='/'))
options(scientific.helpers.only=TRUE);source(file.path(ROOT,'R_scientific_revision/run.R'))
# SIR is an effective local percentage offset, not an observed infection process.
sir <- function(beta,a,y,h=4) {
 i=max(1e-6,min(.7,y[1]/100));s=a*(1-i);out=numeric(length(y)+h);out[1]=i*100
 for(t in 2:length(out)){for(j in 1:4){A=s*(-expm1(-beta*i/4));D=i*(-expm1(-.7/4));s=s-A;i=i+A-D};out[t]=100*i;if(t==length(y)){i=max(1e-6,min(.99,tail(y,1)/100));s=min(s,1-i)}}
 out
}
offset <- function(y) {fit=optim(c(.5,.8),function(p)mean((head(sir(p[1],p[2],y),length(y))-y)^2),method='L-BFGS-B',lower=c(.03,.03),upper=c(3.5,.999),control=list(maxit=25));tail(sir(fit$par[1],fit$par[2],y),4)}
panel <- function(d) {
 out=list()
 for(r in unique(d$region)) {
  z=d[d$region==r,];z=z[order(z$time),];n=nrow(z);z$t=seq_len(n)
  # Do not fill response values. Past input values may be carried at most 2 rows.
  yy=z$y;age=0
  for(j in seq_len(n)){if(is.na(yy[j])){age=age+1;if(j>1&&age<=2)yy[j]=yy[j-1]}else age=0}
  for(j in 12:(n-4)) {
   hist=yy[(j-11):j];if(anyNA(hist))next
   m=offset(hist)
   for(h in c(1,2,4)) {
    if(is.na(z$y[j+h]))next
    out[[length(out)+1]]=data.frame(region=r,time=z$time[j],response_time=z$time[j+h],t=j,h=h,y=z$y[j+h],lag0=yy[j],lag1=yy[j-1],lag2=yy[j-2],lag3=yy[j-3],avg=mean(yy[(j-3):j]),sine=sin(2*pi*j/52),cosine=cos(2*pi*j/52),sir=m[h])
   }
  }
 }
 do.call(rbind,out)
}
forecast_one <- function(d,target,origin,h,seed) {
 features=c('lag0','lag1','lag2','lag3','avg','sine','cosine');d=d[d$h==h,]
 te=d[d$region==target&d$time==origin,];if(nrow(te)!=1)return(NULL)
 tr=tail(d[d$region==target&d$response_time<=origin,],80);if(nrow(tr)<60)return(NULL)
 sources=split(d[d$region!=target&d$response_time<=origin,],d$region[d$region!=target&d$response_time<=origin]);sources=lapply(sources,tail,120)
 if(length(sources)!=2||any(vapply(sources,nrow,1L)<60))return(NULL)
 x=tr[,features];q=te[,features];y=tr$y
 bank=lapply(seq_along(sources),function(k)rf(sources[[k]][,features],sources[[k]]$y-sources[[k]]$sir,seed+k))
 P=sapply(bank,pr,x=x);Q=matrix(sapply(bank,pr,x=q),nrow=1);colnames(P)=colnames(Q)=c('source1','source2')
 rf0=pr(rf(x,y,seed+10),q)
 local=te$sir+pr(rf(cbind(x,sir=tr$sir),y-tr$sir,seed+11),cbind(q,sir=te$sir))
 ftrf_sir=te$sir+pr(rf(cbind(x,sir=tr$sir,P),y-tr$sir,seed+12),cbind(q,sir=te$sir,Q))
 persist=te$lag0+pr(rf(x,y-tr$lag0,seed+13),q)
 srcsir=d[d$region!=target&d$time==origin,'sir'];if(length(srcsir)!=2)return(NULL)
 pred=c(Persistence=te$lag0,SIR=te$sir,RF_tar=rf0,SIR_residual=local,`FTRF-SIR`=ftrf_sir,Persistence_residual=persist,Auxiliary_SIR_mean=mean(srcsir),All_SIR_mean=mean(c(te$sir,srcsir)),Half_SIR_RF=(te$sir+rf0)/2)
 data.frame(target=target,origin=origin,test_time=te$response_time,h=h,model=names(pred),observed=te$y,predicted=pmin(100,pmax(0,pred)),train_n=nrow(tr),source_n=sum(vapply(sources,nrow,1L)),seed=seed)
}
if(!isTRUE(getOption('scientific.definitions.only',FALSE))) {
if(!('--simulation-only'%in%commandArgs(TRUE))) {
 states=c('Alabama','Alaska','Arizona','Arkansas');ili=readraw('ILINet.csv');cls=readraw('WHO_NREVSS_Clinical_Labs.csv');cpt=readraw('WHO_NREVSS_Combined_prior_to_2015_16.csv')
 raws=list(ILI=data.frame(region=ili$REGION,time=ili$YEAR*100+ili$WEEK,y=ili$X.UNWEIGHTED.ILI),CLS_B=data.frame(region=cls$REGION,time=cls$YEAR*100+cls$WEEK,y=cls$PERCENT.B),CPT_B_rate=data.frame(region=cpt$REGION,time=cpt$YEAR*100+cpt$WEEK,y=100*cpt$B/cpt$TOTAL.SPECIMENS))
 all=list();panels=list()
 for(task in names(raws)) {
  d=raws[[task]];d$y=suppressWarnings(as.numeric(d$y));d$y[!is.finite(d$y)|d$y<0|d$y>100]=NA;d=d[d$region%in%states,]
  d=do.call(rbind,lapply(split(d,d$region),tail,220));selected=head(names(which(table(d$region[!is.na(d$y)])>=180)),3);d=d[d$region%in%selected,];stopifnot(length(selected)==3);p=panel(d);panels[[task]]=p
  for(target in selected)for(h in c(1,2,4)) {
   orig=tail(p$time[p$region==target&p$h==h],52);if(!length(orig))next;orig=orig[seq(1,length(orig),by=4)]
   for(o in orig){v=forecast_one(p,target,o,h,20260921+o+h);if(!is.null(v)){v$task=task;all[[length(all)+1]]=v}}
  }
  savecsv(do.call(rbind,all),'forecast_predictions');cat('Forecast',task,'complete\n');flush.console()
 }
 saveRDS(panels,file.path(OUT,'forecast_panels.rds'))
}
}
# Independent recurrent epidemic episodes. Later origins may use earlier test-episode observations.
simulate_series <- function(seed,beta,reporting) {
 set.seed(seed);Y=c();latent=c()
 for(season in 1:4){s=.92;i=.008;z=numeric(40);for(t in 1:40){for(j in 1:4){A=s*(-expm1(-beta*i/4));D=i*(-expm1(-.7/4));s=s-A;i=i+A-D};z[t]=100*i};Y=c(Y,pmin(100,pmax(0,reporting*z+rnorm(40,0,.2))));latent=c(latent,z)}
 structure(Y,latent=latent)
}
if(!isTRUE(getOption('scientific.definitions.only',FALSE))) {
sim=list()
for(rep in 1:20)for(condition in c('matched','shifted')) {
 beta=1.15+.02*(rep%%5);raw=do.call(rbind,lapply(1:3,function(k)data.frame(region=c('target','source1','source2')[k],time=1:160,y=simulate_series(20260921+rep*10+k,if(k==1||condition=='matched')beta else 1.65,if(k==1||condition=='matched')1 else .4))))
 p=panel(raw)
 for(h in c(1,4))for(o in c(128,136,144,152)) {v=forecast_one(p,'target',o,h,20260921+rep*100+o+h);if(!is.null(v)){v$replicate=rep;v$condition=condition;sim[[length(sim)+1]]=v}}
 savecsv(do.call(rbind,sim),'simulation_predictions');cat('Simulation',rep,condition,'complete\n');flush.console()
}
}
