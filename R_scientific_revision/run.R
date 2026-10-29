# All tuning and stacking use observations available at the historical cutoff.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grep('^--file=',args)][1])
ROOT <- getOption('scientific.root',normalizePath(file.path(dirname(script),'..'), winslash='/'))
.libPaths(c(file.path(ROOT,'R_scientific_revision/library'),.libPaths()))
suppressPackageStartupMessages(library(randomForest))
OUT <- file.path(ROOT,'R_scientific_revision/results'); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
set.seed(20260921)
savecsv <- function(x,n) write.csv(x,file.path(OUT,paste0(n,'.csv')),row.names=FALSE,na='')
rf <- function(x,y,seed) {set.seed(seed); randomForest(x=as.data.frame(x),y=y,ntree=80,nodesize=3,mtry=max(1,floor(ncol(x)/3)))}
pr <- function(m,x) as.numeric(predict(m,as.data.frame(x)))
kernel <- function(a,b) exp(-pmax(outer(rowSums(a*a),rowSums(b*b),'+')-2*a%*%t(b),0)/max(1,ncol(a)))
krr <- function(x,y,q) {mu=colMeans(x); sc=apply(x,2,sd);sc[!is.finite(sc)|sc<1e-8]=1;x=sweep(sweep(x,2,mu),2,sc,'/');q=sweep(sweep(q,2,mu),2,sc,'/'); ym=mean(y); as.numeric(ym+kernel(q,x)%*%solve(kernel(x,x)+diag(.1,nrow(x)),y-ym))}
# Greedy convex aggregation uses the specified zero prior penalty.
gma <- function(P,y,Q) {f=rep(0,length(y));w=rep(0,ncol(P));for(l in 1:100){a=(l-1)/l;mu=if(l==1)0 else if(l==2).05 else .01*(l-1)/l^2;loss=colSums((a*f+(1-a)*P-y)^2)+mu*colSums((f-P)^2);j=which.min(loss);w=a*w;w[j]=w[j]+1-a;f=as.numeric(P%*%w)};as.numeric(Q%*%w)}
compare <- function(tr,src,te,features,seed) {
  x=as.matrix(tr[,features,drop=FALSE]);y=tr$y;q=as.matrix(te[,features,drop=FALSE]);n=nrow(x)
  core=seq_len(max(2,floor(.65*n)));cal=setdiff(seq_len(n),core)
  sources=split(src,src$region);sources=sources[vapply(sources,nrow,1L)>=4]
  if(!length(sources)||length(cal)<2) stop('Insufficient stacking records')
  bank=c(list(rf(x[core,,drop=FALSE],y[core],seed)),lapply(seq_along(sources),function(k)rf(sources[[k]][,features,drop=FALSE],sources[[k]]$y,seed+k)))
  P=sapply(bank,pr,x=x[cal,,drop=FALSE]);Q=sapply(bank,pr,x=q)
  if(is.null(dim(Q)))Q=matrix(Q,nrow=nrow(q))
  # Target meta learners see only the held-out tail; base target fit excludes it.
  meta=rf(cbind(x[cal,,drop=FALSE],P[,-1,drop=FALSE]),y[cal],seed+10)
  ftrf=pr(meta,cbind(q,Q[,-1,drop=FALSE]))
  pm=colMeans(P);ps=apply(P,2,sd);ps[!is.finite(ps)|ps<1e-8]=1
  Z=cbind(1,sweep(sweep(P,2,pm),2,ps,'/'));W=cbind(1,sweep(sweep(Q,2,pm),2,ps,'/'))
  sr=as.numeric(W%*%solve(crossprod(Z)+diag(c(0,rep(.1,ncol(P)))),crossprod(Z,y[cal])))
  all=rbind(tr,src);all=all[order(all$time),];if(nrow(all)>600)all=tail(all,600)
  sp=src[order(src$time),];if(nrow(sp)>312)sp=tail(sp,312)
  sm=rf(sp[,features,drop=FALSE],sp$y,seed+20);res=rf(x,y-pr(sm,x),seed+21)
  kp=krr(as.matrix(sp[,features,drop=FALSE]),sp$y,rbind(x,q));kr=kp[(n+1):(n+nrow(q))]+krr(x,y-kp[seq_len(n)],q)
  ans=cbind(RF_tar=pr(rf(x,y,seed+30),q),FTRF=ftrf,SR=sr,GMA=gma(P,y[cal],Q),RF_agg=pr(rf(all[,features,drop=FALSE],all$y,seed+40),q),TL_RF=pr(sm,q)+pr(res,q),TL_KRR=kr)
  pmax(ans,0)
}
readraw <- function(name) {p=file.path(ROOT,'FTRF_SIR_Rebuilt/data/raw',name);d=read.csv(p,skip=1,check.names=TRUE,na.strings=c('X','NA',''));d=d[d$REGION.TYPE=='States',];d}
if (!isTRUE(getOption('scientific.helpers.only',FALSE))) {
states=c('Alabama','Alaska','Arizona','Arkansas')
ili=readraw('ILINet.csv');cls=readraw('WHO_NREVSS_Clinical_Labs.csv');cpt=readraw('WHO_NREVSS_Combined_prior_to_2015_16.csv');phl=readraw('WHO_NREVSS_Public_Health_Labs.csv')
make <- function(d,y,features,season=FALSE) {
  tm=if(season)as.integer(sub('Season ([0-9]+).*','\\1',d$SEASON_DESCRIPTION)) else d$YEAR*100+d$WEEK
  out=data.frame(region=d$REGION,time=tm,y=suppressWarnings(as.numeric(d[[y]])))
  for(k in seq_along(features))out[[paste0('x',k)]]=suppressWarnings(as.numeric(d[[features[k]]]))
  out=out[out$region%in%states,];out=out[order(out$region,out$time),];stopifnot(!anyDuplicated(out[c('region','time')]))
  out
}
tasks=list(ILI=make(ili,'X.UNWEIGHTED.ILI',c('YEAR','WEEK','NUM..OF.PROVIDERS')),
 CLS_B=make(cls,'PERCENT.B',c('YEAR','WEEK','TOTAL.SPECIMENS','PERCENT.A')),
 CPT_B=make(cpt,'B',c('YEAR','WEEK','TOTAL.SPECIMENS','A..2009.H1N1.','A..H1.','A..H3.','A..Subtyping.not.Performed.','A..Unable.to.Subtype.','H3N2v')),
 PHL_A=make(phl,'A..Subtyping.not.Performed.',c('TOTAL.SPECIMENS','A..2009.H1N1.','A..H3.'),TRUE),
 PHL_B=make(phl,'B',c('TOTAL.SPECIMENS','BYam','BVic'),TRUE))
records=list();exclusions=list();selection=list();counter=0
for(task in names(tasks)) {
 d=tasks[[task]];features=grep('^x',names(d),value=TRUE);isphl=grepl('PHL',task)
 eligible=names(which(table(d$region[complete.cases(d)])>=if(isphl)7 else 130)); selected=head(eligible,3)
 selection[[task]]=data.frame(task,region=selected);d=d[d$region%in%selected,]
 for(target in selected) {
  local=d[d$region==target,];local=local[complete.cases(local),];n=nrow(local)
  origins=if(isphl)seq.int(6,n-1) else unique(round(seq(max(104,floor(.7*n)),n-4,length.out=100)))
  for(o in origins) {
   cutoff=local$time[o];tr=tail(local[seq_len(o),],104);te=local[(o+1):min(n,o+if(isphl)1 else 4),]
   src=d[d$region!=target&d$time<=cutoff,];src=src[complete.cases(src),];src=do.call(rbind,lapply(split(src,src$region),tail,156))
   counter=counter+1;seed=20260921+counter
   pred=tryCatch(compare(tr,src,te,features,seed),error=function(e)e)
   if(inherits(pred,'error')) {exclusions[[length(exclusions)+1]]=data.frame(task,target,cutoff,reason=conditionMessage(pred));next}
   for(m in colnames(pred))records[[length(records)+1]]=data.frame(task,target,cutoff,test_time=te$time,model=m,observed=te$y,predicted=pred[,m],seed=seed,train_n=nrow(tr),source_n=nrow(src))
  }
  cat(task,target,'complete\n');flush.console()
 }
 savecsv(do.call(rbind,records),'conditional_predictions')
}
savecsv(if(length(exclusions))do.call(rbind,exclusions)else data.frame(task=character(),target=character(),cutoff=numeric(),reason=character()),'exclusions')
savecsv(do.call(rbind,selection),'region_selection')
saveRDS(tasks,file.path(OUT,'prepared_conditional.rds'))
capture.output(sessionInfo(),file=file.path(OUT,'R_session.txt'))
cat('CONDITIONAL COMPLETE\n')
}
