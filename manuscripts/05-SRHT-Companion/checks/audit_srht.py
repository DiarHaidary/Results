"""Independent scalar and whole-sketch checks for the SRHT companion.

Run: python checks/audit_srht.py
Requires mpmath and numpy. Outputs audit_results.json next to this script.
Numerical checks supplement the proofs and do not formally verify them.
"""
from pathlib import Path
import itertools
import json
import math
from fractions import Fraction as F

import mpmath as mp
import numpy as np

mp.mp.dps = 80
OUT = Path(__file__).resolve().parent

def allowances(q, d, alpha=1, beta=1):
    a=mp.mpf(alpha)*q*(2*q+1)
    b=mp.mpf(beta)*(d+8*q*q)+mp.mpf(alpha)*q*(2*q-1)
    return a+b,(mp.sqrt(a)+mp.sqrt(2*b))**2

def atom(n,m):
    rho=mp.mpf(m)/n
    return mp.exp(mp.loggamma(n)-mp.loggamma(m+1)-mp.loggamma(n-m)
                  +m*mp.log(rho)+(n-1-m)*mp.log1p(-rho))

def certificate(n,d,q,m,eps):
    if m==n: return mp.mpf(0)
    _,lam=allowances(q,d)
    rho=mp.mpf(m)/n
    root=(2*mp.sqrt((1-rho)*lam/m)+abs(1-2*rho)*lam/m)/(1-atom(n,m))
    return d*(root/eps)**(2*q)

def order(d,delta,factor):
    return int(mp.ceil(mp.log(factor*d/delta,4)))

def width(coeff,q,d,eps):
    D,L=allowances(q,d)
    return int(mp.ceil(coeff*(L if coeff in (20,72) else D)/eps**2))

def bisect_width(n,d,q,eps,delta):
    if certificate(n,d,q,n-1,eps)>delta:return n
    lo,hi=1,n-1
    while lo<hi:
        mid=(lo+hi)//2
        if certificate(n,d,q,mid,eps)<=delta:hi=mid
        else:lo=mid+1
    assert certificate(n,d,q,lo,eps)<=delta
    return lo

def scalar_checks():
    n=2**30;eps=mp.mpf('0.5');delta=mp.mpf('0.01')
    expected={10:(6,6,381888,369386,106080,102608,58396,4),
              100:(8,7,749952,717818,165120,157664,125750,6),
              1000:(10,9,1900800,1733312,473280,427412,346315,10),
              10000:(11,11,9894528,7894646,2748480,2192958,1423395,18)}
    main=[]
    for d, exp in expected.items():
        qo=order(d,delta,4);qs=order(d,delta,mp.mpf(6)/5)
        pairs=[(bisect_width(n,d,q,eps,delta),q) for q in range(1,41)]
        m,q=min(pairs)
        vals=(qo,qs,width(216,qo,d,eps),width(72,qo,d,eps),width(60,qs,d,eps),width(20,qs,d,eps),m,q)
        assert vals==exp,(d,vals,exp)
        assert certificate(n,d,q,m-1,eps)>delta
        main.append({'d':d,'q_original':qo,'q_star':qs,'counts':list(vals[2:6]),
                     'bisection_candidate':m,'selected_order':q,
                     'certificate_before':float(certificate(n,d,q,m-1,eps)),
                     'certificate_at_boundary':float(certificate(n,d,q,m,eps)),
                     'certificate_at_returned':float(certificate(n,d,q,m+10,eps)),
                     'certificate_at_canonical':float(certificate(n,d,qs,vals[5],eps))})
    # Independent exact comparison on an explicit finite grid of all binomial
    # atoms. This is not a proof beyond that grid.
    largest=mp.mpf(0);arg=None
    for nn in range(2,101):
        for m in range(1,nn):
            p=atom(nn,m);v=mp.mpf(m)*(1-mp.mpf(m)/nn)
            bound=min(mp.mpf('0.5'),mp.sqrt(mp.pi/(8*v)))
            ratio=p/bound
            assert ratio<=1+mp.mpf('1e-70')
            if ratio>largest:largest=ratio;arg=(nn,m)
    gate=[]
    for price in [1,252,math.factorial(10),math.comb(1000,10)]:
        dd=delta/price; qs=order(100,dd,mp.mpf(6)/5)
        m,q=min((bisect_width(n,100,q,eps,dd),q) for q in range(1,61))
        gate.append({'price':str(price),'q_star':qs,'canonical':width(20,qs,100,eps),
                     'coherence':width(60,qs,100,eps),'direct_plus_10':m+10,'direct_q':q})
    _,L=allowances(1,1)
    r2=mp.mpf(8)/5*mp.sqrt(L)
    r3=mp.mpf(64)/37*(mp.sqrt(L/3)+L/6)
    assert r3>r2
    return {'main':main,'gate':gate,'atom_grid_max_ratio':float(largest),'atom_argmax':arg,
            'nonmonotonic_root_example':{'n':4,'d':1,'q':1,'Lambda':float(L),'root_M2':float(r2),'root_M3':float(r3)},
            'kappa0':float(mp.exp(2-mp.euler)/2),
            'scope':'80-decimal boundary evaluations; bisection candidates are checked, not globally proved minimal'}

def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def transpose(a):return [list(row) for row in zip(*a)]
def trace(a):return sum(a[i][i] for i in range(len(a)))
def fro2(a):return sum(x*x for row in a for x in row)

def enumerate_n4():
    """Exact rational Gram energies and covariance maps, all n=4 laws."""
    n=4
    h=[[F(1,2)*v for v in row] for row in [[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]]
    # Explicit fixtures make these checks independent of Julia's seeded RNG.
    c=[[1,2,0],[-2,1,1],[0,-1,2],[2,0,-1]]
    k=matmul(c,transpose(c))
    b=[[2,1,0,1],[1,-1,2,0],[0,1,1,-2],[-1,0,2,1]]
    X=matmul(b,transpose(b))
    signs=list(itertools.product([-1,1],repeat=n));results=[]
    for two in [False,True]:
        for m in [1,2,3,4]:
            energy=F(0);cov=[[F(0) for _ in range(n)] for _ in range(n)]
            residual=0.;largest_violation=0.;count=0
            for x in signs:
                for y in signs if two else [(1,)*n]:
                    if two:
                        hy=[[h[i][j]*y[j] for j in range(n)] for i in range(n)]
                        q=matmul(hy,h)
                    else:q=h
                    q=[[q[i][j]*x[j] for j in range(n)] for i in range(n)]
                    for subset in itertools.combinations(range(n),m):
                        rr=[q[j] for j in subset]
                        g=matmul(transpose(rr),rr)
                        e=[[F(n,m)*g[i][j]-(i==j) for j in range(n)] for i in range(n)]
                        ex=matmul(matmul(e,X),e)
                        energy+=trace(matmul(matmul(matmul(k,e),k),e))
                        for i in range(n):
                            for j in range(n):cov[i][j]+=ex[i][j]
                        sf=np.asarray(rr,float)*math.sqrt(n/m);kf=np.asarray(k,float)
                        ky=kf@sf.T
                        res=kf-ky@np.linalg.pinv(sf@ky)@ky.T
                        rf=float(np.linalg.norm(res,'fro')**2)
                        residual+=rf
                        gram=float(trace(matmul(matmul(matmul(k,e),k),e)))
                        largest_violation=max(largest_violation,rf-gram)
                        count+=1
            energy/=count;cov=[[v/count for v in row] for row in cov]
            lam=F(n-m,m*(n-1));t=trace(k);ff=fro2(k);diag=sum(k[i][i]**2 for i in range(n))
            predicted=lam*((1-F(2,n))*t*t+(1-F(4,n))*ff+F(4,n)*diag) if two else lam*(t*t+ff-2*diag)
            assert energy==predicted
            predicted_cov=[[lam*((1-F(4,n))*X[i][j]+(1-F(2,n))*trace(X)*(i==j)+F(4,n)*X[i][i]*(i==j)) if two else lam*(X[i][j]+trace(X)*(i==j)-2*X[i][i]*(i==j)) for j in range(n)] for i in range(n)]
            assert cov==predicted_cov
            bound=2*F(1,m)*(1-F(m,n))*t*t
            assert residual/count<=float(bound)+1e-10
            assert largest_violation<1e-8
            results.append({'law':'two-round' if two else 'one-round','n':n,'M':m,'states':count,
                            'exact_energy':str(energy),'formula_exact':True,'covariance_exact':True,
                            'mean_nystrom_residual_float':residual/count,'bound':str(bound),
                            'maximum_deterministic_violation':largest_violation})
    return {'K':k,'X':X,'cases':results}

if __name__=='__main__':
    report={'scalar':scalar_checks(),'whole_sketch':enumerate_n4(),
            'verification':'finite scalar and exact rational law checks; not formal verification'}
    (OUT/'audit_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS: main scalar table, gate table, binomial atoms n<=100, exact n=4 energy/covariance laws, deterministic residual comparisons')
    print('Wrote',OUT/'audit_results.json')
