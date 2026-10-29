"""Conservative discrete SIR solver and historical-window mechanism features."""
import numpy as np
from scipy.optimize import least_squares
try:
    from numba import njit
except ImportError:
    def njit(*args,**kwargs):
        return lambda f:f

@njit(cache=True)
def sir(beta, gamma, s0, i0, weeks, substeps=8):
    """Exponential competing transfers preserve positivity and S+I+R=1."""
    s, i = float(s0), float(i0)
    out = np.empty((weeks+1,3))
    out[0,0],out[0,1],out[0,2]=s,i,1-s-i
    for t in range(weeks):
        for _ in range(substeps):
            infection = s * (-np.expm1(-beta*i/substeps))
            removal = i * (-np.expm1(-gamma/substeps))
            s, i = s-infection, i+infection-removal
        out[t+1,0],out[t+1,1],out[t+1,2]=s,i,1-s-i
    return out

def fit_forecast(y, max_h=4):
    """Rates are weekly; observed percentages are working state proxies only."""
    y = np.asarray(y, float)
    i0 = np.clip(y[0]/100, 1e-6, .7)
    gamma = .7
    def residual(p):
        trajectory = sir(p[0], gamma, p[1]*(1-i0), i0, len(y)-1)
        return (100*trajectory[:,1]-y)/max(np.std(y), .2)
    fit = least_squares(residual, [.95, .8], bounds=([.03,.03],[3.5,.999]),
                        max_nfev=35, ftol=1e-5, xtol=1e-5, gtol=1e-5)
    z = sir(fit.x[0],gamma,fit.x[1]*(1-i0),i0,len(y)-1+max_h)
    # Reanchor I to the latest observation and keep the fitted state feasible.
    latest_i = np.clip(y[-1]/100,1e-6,.99)
    latest_s = min(z[len(y)-1,0], 1-latest_i)
    future = sir(fit.x[0],gamma,latest_s,latest_i,max_h)[1:,1]*100
    return future, fit.x[0], latest_s, bool(fit.success)
