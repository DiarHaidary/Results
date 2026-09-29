"""Independent targeted checks. Python >=3.11; numpy, scipy and sympy.

These diagnostics do not prove the general inequalities and do not reproduce
the original multi-million-path Julia tables. Run from any working directory.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import json
import numpy as np
from scipy.integrate import quad
import sympy as sp

OUT = Path(__file__).resolve().parent / 'independent_audit_results.json'
SEED = 20260930
rng = np.random.default_rng(SEED)

def esym(lam):
    e=np.zeros(len(lam)+1); e[0]=1.
    for x in lam:
        e[1:] += x*e[:-1].copy()
    return e

def volume(lam,k):
    e=esym(lam)
    return (k+1)*e[k+1]/e[k] if k<len(lam) else 0.

def covariance(alpha,c,u,v):
    c=np.asarray(c); pi=c/c.sum()
    a=np.asarray(u)-pi@u; b=np.asarray(v)-pi@v
    def f(s):
        t=np.exp(s); ct=c/(1+t*c)
        phi=np.exp(-alpha*np.log1p(t*c).sum())
        return t*phi*(alpha**2*(ct@a)*(ct@b)+alpha*((ct*ct*a)@b))
    ans,err=quad(f,-55,55,epsabs=2e-12,epsrel=2e-10,limit=300)
    return float(ans),float(err)

def onestep(lam,k,alpha):
    lam=np.asarray(lam)
    deleted=[esym(np.delete(lam,i)) for i in range(len(lam))]
    c=lam*np.array([e[k-1] for e in deleted])
    u=lam/c
    rho=np.array([e[k]/e[k-1] for e in deleted])
    D,err=covariance(alpha,c,u,rho)
    lhs=k*((c/c.sum())@rho+D/(alpha*lam.sum()))
    return dict(spectrum=lam.tolist(),k=k,alpha=alpha,lhs=float(lhs),rhs=volume(lam,k),D=D,quad_error=err)

def spectral_chains(lam,alpha,paths=100000):
    vals=np.tile(lam,(paths,1)).astype(float)
    results=[]
    for k in range(1,len(lam)):
        total=vals.sum(axis=1)
        # Energy tilting is a mixture of one-coordinate Gamma size biases.
        ix=(rng.random(paths)[:,None] >= np.cumsum(vals/total[:,None],axis=1)).sum(axis=1)
        g=rng.gamma(alpha,1,size=vals.shape)
        g[np.arange(paths),ix]+=rng.exponential(size=paths)
        w=vals*np.sqrt(g); q=(vals*g).sum(axis=1)
        B=-w[:,:,None]*w[:,None,:]/q[:,None,None]
        B[:,np.arange(vals.shape[1]),np.arange(vals.shape[1])]+=vals
        vals=np.linalg.eigvalsh(B)[:,1:]
        vals=np.maximum(vals,0)
        trace=vals.sum(axis=1)
        results.append(dict(k=k,mean=float(trace.mean()),se=float(trace.std(ddof=1)/np.sqrt(paths)),volume=volume(lam,k)))
    return dict(spectrum=list(lam),alpha=alpha,paths=paths,steps=results)

def exact_coordinate(r,den):
    n=r+1; h=sp.Rational(1,den); eta=h**(n*n)
    L=sp.eye(n)
    for a in range(n):
        for b in range(a):
            L[a,b]=h**(b+1)/(a+b+2)
    A=L*sp.diag(*[eta**j for j in range(n)])*L.T
    minors={0:sp.Integer(1)}
    for mask in range(1,1<<n):
        ids=[j for j in range(n) if mask>>j&1]
        minors[mask]=A.extract(ids,ids).det()
    @lru_cache(None)
    def residual(mask):
        return sum(minors[mask|1<<j]/minors[mask] for j in range(n) if not mask>>j&1)
    @lru_cache(None)
    def mean(mask):
        if mask.bit_count()==r:
            return residual(mask)
        tr=residual(mask)
        return sum(minors[mask|1<<j]/minors[mask]/tr*mean(mask|1<<j) for j in range(n) if not mask>>j&1)
    ratio=mean(0)*A.inv().trace()
    return dict(r=r,h=f'1/{den}',expected_trace_times_inverse_trace=float(ratio),limit=2**r)

spectra=[[5,3,2,1],[10,1,.1,.01,.001],[4,4,1,.5,.25,.1],[100,10,1,1,1]]
results={'seed':SEED,'scope':'Targeted independent diagnostics; not formal verification.',
 'one_step':[onestep(lam,2,a) for lam in spectra for a in [.5,1.]],
 'nonmonotone_counterexample':dict(zip(['D','quad_error'],covariance(1,[.02139,.27215,.03634],[-1.6053,-.1726,-.0118],[2.5635,2.5023,1.1242]))),
 'chains':[spectral_chains(lam,a) for lam in spectra for a in [.5,1.]],
 'coordinate':[exact_coordinate(r,d) for r in [2,3] for d in [10,100]],
 'wishart_inverse_moment_identity':[dict(field=field,r=r,ell=ell,bartlett_product=float(np.prod([1/(ell+1-i) for i in range(1,ell+1)])),claimed=float(1/sp.factorial(ell))) for field in ['real','complex'] for r in range(1,7) for ell in range(1,r+1)]}
assert all(x['lhs']<=x['rhs']+1e-9 for x in results['one_step'])
assert results['nonmonotone_counterexample']['D']>0
assert all(abs(x['bartlett_product']-x['claimed'])<1e-14 for x in results['wishart_inverse_moment_identity'])
OUT.write_text(json.dumps(results,indent=2),encoding='utf-8')
print(f'Wrote {OUT}; {len(results["chains"])} chain experiments and {len(results["one_step"])} quadrature checks.')
