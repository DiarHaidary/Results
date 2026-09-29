from pathlib import Path
import numpy as np
from scipy.special import gammaln,digamma
from scipy.optimize import minimize_scalar
from math import ceil,log,exp
import json
out={}
def lc(n,p,cx=True):
    a=1. if cx else .5
    b=n if cx else n/2
    return p*log(n)+gammaln(p+a)+gammaln(b)-gammaln(a)-gammaln(b+p)
def D(n,cx=True):
    return log(n)+digamma(2 if cx else 1.5)-digamma(n+1 if cx else n/2+1)
constants={}
for n in [2,3,4,5,10,64]:
    for cx in [False,True]:
        opt=minimize_scalar(lambda p: lc(n,p,cx),bounds=(0,1),method="bounded",options={"xatol":1e-14})
        constants[f"{n}-{'complex' if cx else 'real'}"]={"D":D(n,cx),"kappa_factor":exp(D(n,cx)),"gamma":-opt.fun,"lambda":opt.x,"fourth":exp(lc(n,2,cx))}
out["constants"]=constants
out["leading_ratios"]=[(196*32)**(1/3),(2704*32)**(1/3)]
out["kappa0"]=exp(2-np.euler_gamma)/2
out["binary_bases"]={"direct":10*2**.5/9,"fourth_real":(100/27)**(1/3),"fourth_cx":(800/243)**(1/3),"drift_real":(100*exp(1)/81)**(1/3),"drift_cx":(400/(81*exp(.5)))**(1/3),"gamma_cx":(200*exp(constants["2-complex"]["gamma"])/81)**(1/3),"gamma_real":(200*exp(constants["2-real"]["gamma"])/81)**(1/3),"analysis_benchmark":(200/81)**(1/3)}
def drift(n,q,eps,cx):
    b=q+1 if cx else 1
    L=(6*b*n**q*exp(q*D(n,cx))/eps**2)**(1/3)
    return b*ceil(L/b)+2*ceil(L)
def empirical(n,q,eps):
    def obj(p): return (log(6)+q*lc(n,p,True)-p*log(eps))/(p-1)
    opt=minimize_scalar(obj,bounds=(1+1e-7,2),method="bounded")
    v=min(opt.fun,obj(2))
    return (q+1)*ceil(exp(v))
def head_tail(n,q,eps,b,k):
    # Integer minima: exhaustive when manageable; otherwise enough neighbors
    # to verify the rounded display rather than certify a global integer minimum.
    from math import sqrt
    alpha=b*k
    rstar=(sqrt(192*n**q)/(eps*alpha))**(2/3)
    maxr=n**q-1
    lo=max(1,int(rstar)-100)
    hi=min(maxr,int(rstar)+101)
    candidates=set(range(lo,hi+1))|{1,maxr}
    return min(b*ceil(k*(r+2*log(24)))+2*ceil(sqrt(192*n**q/(r*eps**2))) for r in candidates)
def fourth(n,q,eps):
    return min(head_tail(n,q,eps,1,196*exp(q*lc(n,2,False))),head_tail(n,q,eps,q+1,196*exp(q*lc(n,2,True))))
def joint(n,q,eps):
    def lk(p): return (p*log(52)+q*lc(n,p,True))/(p-1)
    opt=minimize_scalar(lk,bounds=(1+1e-7,2),method="bounded",options={"xatol":1e-12})
    k=exp(min(opt.fun,lk(2)))
    return head_tail(n,q,eps,q+1,k)
rows=[]
for q in [8,12,16,20,24,30,40]:
    eps=.9**q
    rows.append({"q":q,"N":2**q,"drift_real":drift(2,q,eps,False),"drift_cx":drift(2,q,eps,True),"fourth":fourth(2,q,eps),"joint":joint(2,q,eps),"empirical":empirical(2,q,eps),"direct":2*ceil((6*2**q)**.5/eps)})
out["budget_rows"]=rows
out["crossovers"]={}
for epsname,fun in [("0.9^q",lambda q:.9**q),("0.1",lambda q:.1)]:
    winners=[q for q in range(2,61) if drift(2,q,fun(q),False)<=drift(2,q,fun(q),True)]
    out["crossovers"][epsname]=winners
out["moment_derivative_max_error"]=max(abs((lc(n,1+1e-5,cx)-lc(n,1-1e-5,cx))/(2e-5)-D(n,cx)) for n in [2,3,5,10] for cx in [False,True])
rng=np.random.default_rng(20260929)
mc=[]
# Cheap exact formulas in terms of local invariants, avoiding tensor expansion.
for q in [4,6,8,10,12,14]:
    samples=20000
    w=rng.normal(size=(samples,q,2))+1j*rng.normal(size=(samples,q,2))
    w=w/np.linalg.norm(w,axis=-1,keepdims=True)*np.sqrt(2)
    z0=np.prod(w[:,:,0],axis=1)
    X=abs(z0)**2
    bilinear=np.prod(np.sum(w*w,axis=-1),axis=1)
    N=2**q;delta=1/(N-1);scale=(1+1/(N-1))/2
    scalar=((1-delta*delta)*X+delta*delta*N)/((1-delta)*X+delta*N)
    g=np.empty((samples,2,2))
    g[:,0,0]=(N+bilinear.real)/2
    g[:,1,1]=(N-bilinear.real)/2
    g[:,0,1]=g[:,1,0]=bilinear.imag/2
    y=np.stack([z0.real,z0.imag],axis=1)
    small=delta*g+(1-delta)*y[:,:,None]*y[:,None,:]
    solved=np.linalg.solve(small,y[:,:,None])[:,:,0]
    overlap=np.einsum("si,si->s",y,solved)
    block=2*delta+(1-delta)*overlap
    mc.append({"q":q,"samples":samples,"scalar":float(scalar.mean()/scale),"scalar_se":float(scalar.std(ddof=1)/np.sqrt(samples)/scale),"block":float(block.mean()/scale),"block_se":float(block.std(ddof=1)/np.sqrt(samples)/scale),"proved":exp(-q*D(2,True)),"chernoff":exp(-q*constants["2-complex"]["gamma"])})
out["monte_carlo"]=mc
# Basic matrix identities: no attempt to replace the mathematical proof.
checks=[]
for N in [4,8,12]:
    for r in [1,N//2,N]:
        Q,_=np.linalg.qr(rng.normal(size=(N,N)))
        vals=np.r_[rng.uniform(.1,3,size=r),np.zeros(N-r)]
        A=(Q*vals)@Q.T
        Z=rng.normal(size=(N,2));W=rng.normal(size=(N,2))
        def nys(M,V): return M@V@np.linalg.pinv(V.T@M@V,rcond=1e-10)@V.T@M
        R=A-nys(A,Z)
        upd=R-nys(R,W)
        both=A-nys(A,np.c_[Z,W])
        checks.append(float(np.linalg.norm(upd-both)))
out["sequential_max_error"]=max(checks)
from itertools import product
N,b=3,2
rawB=rng.normal(size=(N,N));rawR=rng.normal(size=(N,N))
B=rawB@rawB.T;R=rawR@rawR.T
mixed=[]
for hashes in product(range(b),repeat=N):
    for signs in product([-1.,1.],repeat=N):
        H=np.zeros((b,N))
        H[list(hashes),np.arange(N)]=signs
        mixed.append(np.trace(H@B@H.T)*np.trace(H@R@H.T))
formula=np.trace(B)*np.trace(R)+2/b*(np.trace(B@R)-np.dot(np.diag(B),np.diag(R)))
out["countsketch_mixed_identity_error"]=float(abs(np.mean(mixed)-formula))
N=4
rawR=rng.normal(size=(N,N));R=rawR@rawR.T
errors=[]
for signs in product([-1.,1.],repeat=4):
    x=np.kron(signs[:2],signs[2:])
    errors.append(np.linalg.norm(x*(R@x)-np.diag(R))**2)
expected=np.linalg.norm(R,"fro")**2-np.dot(np.diag(R),np.diag(R))
out["product_diagonal_training_identity_error"]=float(abs(np.mean(errors)-expected))
assert out["moment_derivative_max_error"]<1.1e-9
assert out["sequential_max_error"]<1e-11
assert out["countsketch_mixed_identity_error"]<1e-10
assert out["product_diagonal_training_identity_error"]<1e-10
assert abs(constants["2-complex"]["gamma"]-(log(2)-1-log(log(2))))<1e-12
assert all(row["drift_real"]<=row["fourth"] and row["drift_cx"]<=row["joint"] for row in rows)
out["status"]="PASS: numerical identities and reported constants; not formal verification"
Path(__file__).with_name("verification_results.json").write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
