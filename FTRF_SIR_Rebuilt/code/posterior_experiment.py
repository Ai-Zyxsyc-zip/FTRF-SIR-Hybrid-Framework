"""Synthetic SIR: NLL/RMSE/MAE/coverage, Wasserstein/MMD, ESS/Rhat and cost.

This is a distinct parameter-transfer diagnostic, not a Bayesian interpretation
of a random forest. Quadrature, Laplace and MH target the same prior/likelihood.
"""
import time,json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import wasserstein_distance
from scipy.spatial.distance import cdist
from mechanism import sir,njit
from prepare_data import ROOT

@njit(cache=True)
def log_target(z,y,mu,sd):
    b,g=np.exp(z[0]),np.exp(z[1])
    if b<.2 or b>2.2 or g<.12 or g>1.3: return -1e100
    traj=sir(b,g,.97,.005,len(y)-1)
    error=(np.log(100*traj[:,1])-np.log(y))/.15
    return -.5*np.sum(error**2)-.5*np.sum(((z-mu)/sd)**2)

@njit(cache=True)
def grid_log(points,y,mu,sd):
    out=np.empty(len(points))
    for j in range(len(points)):
        z=np.log(points[j]); out[j]=log_target(z,y,mu,sd)-np.sum(z)
    return out

@njit(cache=True)
def mh_chain(y,mu,sd,start,chol,seed,burn=800,draws=1600):
    np.random.seed(seed); z=start.copy(); lp=log_target(z,y,mu,sd)
    out=np.empty((draws,2)); accepted=0
    for j in range(burn+draws):
        proposal=z+chol@np.random.normal(0.,1.,2)*1.7
        value=log_target(proposal,y,mu,sd)
        if np.log(np.random.random())<value-lp:
            z=proposal; lp=value
            if j>=burn: accepted+=1
        if j>=burn: out[j-burn]=np.exp(z)
    return out,accepted/draws

def hessian(fun,z,step=1e-4):
    eye=np.eye(2)*step; out=np.empty((2,2))
    for i in range(2):
        for j in range(2):
            out[i,j]=(fun(z+eye[i]+eye[j])-fun(z+eye[i]-eye[j])-fun(z-eye[i]+eye[j])+fun(z-eye[i]-eye[j]))/(4*step**2)
    vals,vec=np.linalg.eigh(out)
    return (vec*(1/np.maximum(vals,1e-5)))@vec.T

def diagnostic(chains):
    """Conventional split-Rhat and initial-positive-sequence ESS; not rank ESS."""
    m,n,d=chains.shape; half=n//2
    x=np.concatenate([chains[:,:half],chains[:,half:2*half]],axis=0)
    within=np.mean(np.var(x,axis=1,ddof=1),axis=0)
    between=half*np.var(np.mean(x,axis=1),axis=0,ddof=1)
    rhat=np.sqrt(((half-1)/half*within+between/half)/within)
    ess=[]
    for j in range(d):
        ac=[]
        for c in chains[:,:,j]:
            c=c-c.mean(); spec=np.fft.rfft(c,n=2*n)
            cov=np.fft.irfft(spec*np.conj(spec))[:n]/np.arange(n,0,-1)
            ac.append(cov/cov[0])
        rho=np.mean(ac,axis=0); total=0.
        for k in range(1,n-1,2):
            pair=rho[k]+rho[k+1]
            if pair<=0: break
            total+=pair
        ess.append(min(m*n,m*n/max(1.,1+2*total)))
    return min(ess),max(rhat)

def distances(x,ref,rng):
    scale=np.std(ref,axis=0); scale=np.maximum(scale,1e-5)
    wd=np.mean([wasserstein_distance(x[:,j],ref[:,j])/scale[j] for j in range(2)])
    a=x[rng.choice(len(x),400,replace=False)]/scale
    b=ref[rng.choice(len(ref),400,replace=False)]/scale
    both=np.vstack([a,b]); bandwidth=np.median(cdist(both,both,'sqeuclidean'))
    # Use the biased nonnegative MMD estimator, including diagonal terms.
    kernel=lambda u,v:np.exp(-cdist(u,v,'sqeuclidean')/max(2*bandwidth,1e-8))
    return wd,max(0.,float(kernel(a,a).mean()+kernel(b,b).mean()-2*kernel(a,b).mean()))

def main():
    tab=ROOT/'results/tables'; rows=[]; samples=[]; predictions=[]; observations=[]; sources=[]; convergence=[]
    # Compile the shared solver before recording method runtimes.
    sir(1.,.7,.97,.005,24)
    log_target(np.log([1.,.7]),np.ones(17),np.log([1.,.7]),np.ones(2))
    grid_log(np.array([[1.,.7]]),np.ones(17),np.log([1.,.7]),np.ones(2))
    mh_chain(np.ones(17),np.log([1.,.7]),np.ones(2),np.log([1.,.7]),np.eye(2)*.01,1,1,1)
    for scenario in ['matched','shifted']:
        for rep in range(10):
            rng=np.random.default_rng(7000+rep+(scenario=='shifted')*100)
            target=np.exp(rng.normal(np.log([1.05,.65]),[.09,.08]))
            noiseless=sir(*target,.97,.005,24)[:,1]*100
            y=noiseless*np.exp(rng.normal(0,.15,25))
            observations.extend(dict(scenario=scenario,rep=rep,t=t,y=value,latent=noiseless[t],train=t<=16,beta=target[0],gamma=target[1]) for t,value in enumerate(y))
            src=[]; source_start=time.perf_counter()
            for k in range(6):
                center=[1.05,.65] if scenario=='matched' else [.55,.9]
                truth=np.exp(rng.normal(np.log(center),[.09,.08]))
                sy=sir(*truth,.97,.005,32)[:,1]*100*np.exp(rng.normal(0,.15,33))
                fit=minimize(lambda z:-log_target(z,sy,np.log([1.,.65]),np.array([.65,.65])),np.log(center),method='Nelder-Mead',options={'maxiter':300,'xatol':1e-7})
                src.append(fit.x)
                sources.append(dict(scenario=scenario,rep=rep,source=k,true_beta=truth[0],true_gamma=truth[1],fit_beta=np.exp(fit.x[0]),fit_gamma=np.exp(fit.x[1])))
            source_seconds=time.perf_counter()-source_start
            src=np.array(src); tmu=src.mean(axis=0); tsd=np.maximum(src.std(axis=0,ddof=1),.18)
            for prior,mu,sd in [('Local',np.log([1.,.65]),np.array([.65,.65])),('Transfer',tmu,tsd)]:
                started=time.perf_counter(); grid_n=181
                bb,gg=np.meshgrid(np.linspace(.2,2.2,grid_n),np.linspace(.12,1.3,grid_n),indexing='ij')
                points=np.column_stack([bb.ravel(),gg.ravel()]); lp=grid_log(points,y[:17],mu,sd)
                weights=np.exp(lp-logsumexp(lp)); ref=points[rng.choice(len(points),6000,p=weights)]
                grid_seconds=time.perf_counter()-started
                edge=np.zeros((grid_n,grid_n),bool); edge[[0,-1],:]=True; edge[:,[0,-1]]=True
                edge_mass=float(weights[edge.ravel()].sum())
                if rep==0:
                    b2,g2=np.meshgrid(np.linspace(.2,2.2,271),np.linspace(.12,1.3,271),indexing='ij')
                    p2=np.column_stack([b2.ravel(),g2.ravel()]); l2=grid_log(p2,y[:17],mu,sd); w2=np.exp(l2-logsumexp(l2))
                    for j,name in enumerate(['beta','gamma']):
                        convergence.append(dict(scenario=scenario,prior=prior,parameter=name,coarse_mean=np.sum(weights*points[:,j]),fine_mean=np.sum(w2*p2[:,j]),absolute_difference=abs(np.sum(weights*points[:,j])-np.sum(w2*p2[:,j]))))
                started=time.perf_counter(); objective=lambda z:-log_target(z,y[:17],mu,sd)
                opt=minimize(objective,mu.copy(),method='Nelder-Mead',options={'maxiter':500,'xatol':1e-8})
                covariance=hessian(objective,opt.x); chol=np.linalg.cholesky(covariance)
                lap=np.exp(rng.multivariate_normal(opt.x,covariance,10000))
                lap=lap[(lap[:,0]>.2)&(lap[:,0]<2.2)&(lap[:,1]>.12)&(lap[:,1]<1.3)][:6000]
                lap_seconds=time.perf_counter()-started
                started=time.perf_counter(); chains=[]; accepts=[]
                for chain in range(4):
                    initial=opt.x+2*chol@rng.normal(size=2)
                    initial=np.clip(initial,np.log([.201,.121]),np.log([2.199,1.299]))
                    z,acc=mh_chain(y[:17],mu,sd,initial,chol,90000+rep*10+chain)
                    chains.append(z); accepts.append(acc)
                chains=np.array(chains); mh_seconds=time.perf_counter()-started+lap_seconds
                ess,rhat=diagnostic(chains)
                for method,draws,cost in [('Grid',ref,grid_seconds),('Laplace',lap,lap_seconds),('MH',chains.reshape(-1,2),mh_seconds)]:
                    wd,mmd=(0.,0.) if method=='Grid' else distances(draws,ref,rng)
                    latent=np.array([sir(b,g,.97,.005,24)[17:,1]*100 for b,g in draws[:2000]])
                    loglik=-.5*((np.log(y[17:])[None,:]-np.log(latent))/.15)**2-np.log(y[17:])[None,:]-np.log(.15*np.sqrt(2*np.pi))
                    nll=float(np.mean(-logsumexp(loglik,axis=0)+np.log(len(latent))))
                    obsdraws=latent*np.exp(rng.normal(0,.15,latent.shape))
                    mean=obsdraws.mean(axis=0); low,median,high=np.quantile(obsdraws,[.025,.5,.975],axis=0)
                    rows.append(dict(scenario=scenario,rep=rep,prior=prior,method=method,nll=nll,
                        rmse=float(np.sqrt(np.mean((mean-y[17:])**2))),mae=float(np.mean(abs(median-y[17:]))),
                        coverage95=float(np.mean((low<=y[17:])&(y[17:]<=high))),width95=float(np.mean(high-low)),
                        beta_error=abs(draws[:,0].mean()-target[0]),gamma_error=abs(draws[:,1].mean()-target[1]),
                        wasserstein_std=wd,mmd2=mmd,seconds=cost,source_seconds=source_seconds if prior=='Transfer' else 0,
                        ess=ess if method=='MH' else np.nan,rhat=rhat if method=='MH' else np.nan,
                        acceptance=np.mean(accepts) if method=='MH' else np.nan,edge_mass=edge_mass,
                        prior_log_beta=mu[0],prior_log_gamma=mu[1],prior_sd_beta=sd[0],prior_sd_gamma=sd[1]))
                    for j,(b,g) in enumerate(draws):
                        samples.append(dict(scenario=scenario,rep=rep,prior=prior,method=method,draw=j,chain=j//1600 if method=='MH' else -1,beta=b,gamma=g))
                    predictions.extend(dict(scenario=scenario,rep=rep,prior=prior,method=method,t=17+j,y=y[17+j],mean=mean[j],median=median[j],lo95=low[j],hi95=high[j]) for j in range(8))
                print('posterior',scenario,rep,prior,'Rhat',round(rhat,3),flush=True)
    # Write the large draw table once and replace the destination atomically.
    for name,data in [('posterior_metrics',rows),('posterior_samples',samples),('posterior_predictions',predictions),('synthetic_observations',observations),('source_parameters',sources),('grid_convergence',convergence)]:
        temporary=tab/f'{name}.pending.csv'
        pd.DataFrame(data).to_csv(temporary,index=False)
        temporary.replace(tab/f'{name}.csv')

if __name__=='__main__': main()
