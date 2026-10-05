import math,json
from fractions import Fraction as F
from pathlib import Path
N=2**30; EPS=F(1,2); DELTA=F(1,100); SCALE=10**18
def sqr(x,upper):
    x=F(x); a=math.isqrt(x.numerator*SCALE*SCALE//x.denominator)
    if upper and a*a*x.denominator<x.numerator*SCALE*SCALE:a+=1
    return F(a,SCALE)
def lambda_float(d,q,k):
    return (math.sqrt(k*(2*q+1))+math.sqrt(2*(d+8*q*q+k*(2*q-1))))**2
def froot(d,q,M,which):
    rho=M/N; sr=math.sqrt((1-rho)/rho); th=abs(1-2*rho)/rho
    cb=max(.5,1-1/(2*math.sqrt(M*(1-rho))))
    L=[0.]+[min(1.,lambda_float(d,q,k)/N) for k in range(1,q+1)]
    if which=='jacobi_uniform':L=[0.]+[L[q]]*q
    if which=='existing':return (2*sr*math.sqrt(L[q])+th*L[q])/cb
    edge=[sr*math.sqrt(L[k+1]) for k in range(q)];diag=[th*x for x in L]
    x=[1.]+[0.]*q
    for j in range(q):
        x=[diag[k]*x[k]+(edge[k-1]*x[k-1] if k else 0)+(edge[k]*x[k+1] if k<q else 0) for k in range(q+1)]
    return sum(z*z for z in x)**(1/(2*q))/cb
def candidate(d,q,which):
    target=.5*(.01/d)**(1/(2*q)); lo=1;hi=min(N//4,math.ceil(60*(d+12*q*q)/.25))
    assert froot(d,q,hi,which)<=target
    while lo<hi:
        m=(lo+hi)//2
        if froot(d,q,m,which)<=target:hi=m
        else:lo=m+1
    return lo
def certificate(d,q,M,which):
    rho=F(M,N); sr=sqr((1-rho)/rho,True); th=abs(1-2*rho)/rho
    cb=max(F(1,2),1-F(1,2)/sqr(F(M)*(1-rho),False))
    L=[F(0)]+[min(F(1),(sqr(k*(2*q+1),True)+sqr(2*(d+8*q*q+k*(2*q-1)),True))**2/N) for k in range(1,q+1)]
    if which=='jacobi_uniform':L=[F(0)]+[L[q]]*q
    if which=='existing':moment=((2*sr*sqr(L[q],True)+th*L[q])/cb)**(2*q)
    else:
        edge=[sr*sqr(L[k+1],True) for k in range(q)];diag=[th*x for x in L];x=[F(1)]+[F(0)]*q
        for j in range(q):
            x=[diag[k]*x[k]+(edge[k-1]*x[k-1] if k else 0)+(edge[k]*x[k+1] if k<q else 0) for k in range(q+1)]
        moment=sum(z*z for z in x)/cb**(2*q)
    ratio=F(d)*moment/EPS**(2*q)/DELTA
    assert ratio<=1,(d,q,M,which,float(ratio))
    return dict(exact_rational_inequality_passed=True,failure_bound_over_delta=float(ratio),root=float(moment)**(1/(2*q)))
rows=[]
for d,q in [(10,6),(100,8),(1000,10),(10000,11)]:
    x={}
    for which in ['existing','jacobi_uniform','jacobi']:
        c=candidate(d,q,which); M=math.ceil(c*101/100)
        x[which]=dict(float_candidate_minimum=c,certified_M=M,**certificate(d,q,M,which))
    row=dict(n=N,d=d,q=q,epsilon=.5,delta=.01,universal_60_M=240*(d+12*q*q),evaluated=x,reduction_vs_evaluated_existing=1-x['jacobi']['certified_M']/x['existing']['certified_M'])
    rows.append(row)
    print(json.dumps(row))
Path(__file__).with_name('srht_terminal_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
