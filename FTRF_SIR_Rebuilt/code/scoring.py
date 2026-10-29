"""RMSE/MAE, predictive NLL, CRPS, WIS, coverage, width; no parameter posterior claims."""
import numpy as np
from scipy.stats import truncnorm
ALPHAS=np.array([.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9])
QUANTILES=np.unique(np.r_[.5,ALPHAS/2,1-ALPHAS/2])

def score_one(y,loc,scale):
    a,b=-loc/scale,(100-loc)/scale
    dist=truncnorm(a,b,loc=loc,scale=scale)
    mean=float(dist.mean()); median=float(dist.ppf(.5))
    lo=dist.ppf(ALPHAS/2); hi=dist.ppf(1-ALPHAS/2)
    intervals=hi-lo+2/ALPHAS*(lo-y)*(y<lo)+2/ALPHAS*(y-hi)*(y>hi)
    wis=(.5*abs(y-median)+np.sum(ALPHAS/2*intervals))/(len(ALPHAS)+.5)
    probs=(np.arange(399)+.5)/399
    errors=y-dist.ppf(probs)
    crps=2*np.mean(errors*(probs-(errors<0)))
    l95,u95=dist.ppf([.025,.975])
    return dict(mean=mean,median=median,se=(mean-y)**2,ae=abs(median-y),
                nll=-float(dist.logpdf(y)),crps=crps,wis=wis,covered95=float(l95<=y<=u95),
                width95=u95-l95,lo95=l95,hi95=u95),dict(zip(QUANTILES,dist.ppf(QUANTILES)))
