import math, json
from fractions import Fraction as F
from pathlib import Path

N=2**30
EPS=F(1,2)
DELTA=F(1,100)
SCALE=10**18

def sqrt_bound(x, upper):
    x=F(x)
    a=math.isqrt(x.numerator*SCALE*SCALE//x.denominator)
    if upper and a*a*x.denominator<x.numerator*SCALE*SCALE:
        a+=1
    return F(a,SCALE)

def lam(d,j,k):
    return (math.sqrt(k*(2*j+1))+math.sqrt(2*(d+8*j*j+k*(2*j-1))))**2

def root(d,q,M):
    rho=M/N
    u=math.sqrt((1-rho)/rho)
    th=abs(1-2*rho)/rho
    cb=max(.5,1-1/(2*math.sqrt(M*(1-rho))))
    x=[1.]+[0.]*q
    for j in range(1,q+1):
        ell=[0.]+[min(1.,lam(d,j,k)/N) for k in range(1,j+1)]
        edge=[u*math.sqrt(ell[k+1]) for k in range(j)]
        y=[0.]*(q+1)
        for k in range(j+1):
            y[k]=th*ell[k]*x[k]
            if k:y[k]+=edge[k-1]*x[k-1]
            if k<j:y[k]+=edge[k]*x[k+1]
        x=y
    return sum(t*t for t in x)**(1/(2*q))/cb

def candidate(d,q):
    target=.5*(.01/d)**(1/(2*q))
    lo=1
    hi=240*(d+12*q*q)
    assert root(d,q,hi)<=target
    while lo<hi:
        mid=(lo+hi)//2
        if root(d,q,mid)<=target:hi=mid
        else:lo=mid+1
    return lo

def certificate(d,q,M):
    rho=F(M,N)
    u=sqrt_bound((1-rho)/rho,True)
    th=abs(1-2*rho)/rho
    cb=max(F(1,2),1-F(1,2)/sqrt_bound(F(M)*(1-rho),False))
    x=[F(1)]+[F(0)]*q
    for j in range(1,q+1):
        ell=[F(0)]+[min(F(1),(sqrt_bound(k*(2*j+1),True)+sqrt_bound(2*(d+8*j*j+k*(2*j-1)),True))**2/N) for k in range(1,j+1)]
        edge=[u*sqrt_bound(ell[k+1],True) for k in range(j)]
        y=[F(0)]*(q+1)
        for k in range(j+1):
            y[k]=th*ell[k]*x[k]
            if k:y[k]+=edge[k-1]*x[k-1]
            if k<j:y[k]+=edge[k]*x[k+1]
        x=y
    moment=sum(t*t for t in x)/cb**(2*q)
    ratio=F(d)*moment/(EPS**(2*q)*DELTA)
    assert ratio<=1,(d,q,M,float(ratio))
    return dict(exact_rational_inequality_passed=True, failure_bound_over_delta=float(ratio))

rows=[]
for d,q,old in [(10,6,28080),(100,8,60271),(1000,10,169157),(10000,11,1011962)]:
    c=candidate(d,q)
    M=(c*101+99)//100
    row=dict(n=N,d=d,q=q,epsilon=.5,delta=.01,old_terminal_grade_M=old,float_candidate=c,certified_prefix_stage_M=M,reduction=1-M/old,**certificate(d,q,M))
    rows.append(row)
    print(json.dumps(row))
Path(__file__).with_name('srht_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
