"""Independent checks for the SparseStack companion; no formal verification claim.

Run with Python >=3.10 and NumPy. Exact certificate arithmetic uses only stdlib.
"""
from fractions import Fraction as F
from math import isqrt, ceil, lgamma, log, prod, factorial
from itertools import combinations, permutations, product
from pathlib import Path
import json
import re
import numpy as np

ROOT = Path(__file__).resolve().parent
K = 10**12
def sqrt_up(x):
    x = F(x)
    n, d = x.numerator*K*K, x.denominator
    k = isqrt((n+d-1)//d)
    while k*k*d<n:
        k += 1
    return F(k, K)

intervals=[(F(j,40),F(j+1,40)) for j in range(20)]+[(F(1,2),F(1))]
cert=[]
for C,B in [(F(69),F(8)),(F(74),F(15,2)),(F(79),F(7))]:
    for left,right in intervals:
        n=max(int(C/B),int(C/(right*(B+F(1,5)))))+1
        g=(2*sqrt_up(F(n-1,n))+1)*sqrt_up(F(1,n))
        p=F(1)
        for j in range(8):
            u=F(2*j+1,16)
            theta=1-left+left*u
            psi=2*sqrt_up((1+right*u/C)*theta/C)+(2+sqrt_up(2))*sqrt_up(u*theta/(C*B))+theta/C+(1+sqrt_up(2)+g)*u/B
            p*=psi
        lo,hi=0,500000
        while lo+1<hi:
            mid=(lo+hi)//2
            if mid**8 >= p*10**48:
                hi=mid
            else:
                lo=mid
        assert p<F(1,256) and hi<500000 and 9*C<800
        cert.append(dict(C=str(C),B=str(B),left=str(left),right=str(right),buckets=n,T=hi))

# Compare every certificate number with the actual LaTeX table, not merely
# with the threshold. This detects transcription errors in the manuscript.
source=(ROOT.parent/'SparseStack_Constants_Nystrom_Adaptive.tex').read_text(encoding='utf-8')
start=source.index('Interval $[\\ell,r]$')
end=source.index('\\bottomrule',start)
rows=re.findall(r'^\$\[([^]]+)\]\$ & (\d+) & (\d+) & (\d+) & (\d+) & (\d+) & (\d+)\\\\',source[start:end],re.M)
assert len(rows)==21
for j,row in enumerate(rows):
    for pair in range(3):
        rec=cert[21*pair+j]
        assert [rec['buckets'],rec['T']]==list(map(int,row[1+2*pair:3+2*pair]))

# Exact integer ceilings for the displayed gate table.
table=[]
for q in [7,11,18,46,72,69,32]:
    row={'q':q}
    for C,B in [(F(81),F(17,2)),(F(69),F(8))]:
        s=ceil(B*(2*q+1)/F(1,2))
        b=ceil(C*(100+2*q+1)/(s*F(1,4)))
        row[str(C)]=dict(s=s,b=b,m=s*b)
    table.append(row)

# Exhaustive second moments for n=3,b=2,s=1. This tests the exact
# pair covariance, weighted AMM formula, diagonal/trace variances,
# mixed trace identity and convex-potential block inequality.
n,b=3,2
A=np.array([[2.,.4,.1],[.4,1.5,-.2],[.1,-.2,1.]])
eig,U=np.linalg.eigh(A)
X=(U*np.sqrt(eig))@U.T
second=[]; diagonal=[]; trace=[]; drift=[]; convex=[]; mixed=[]
B=A@A
for hashes in product(range(b),repeat=n):
    for signs in product([-1.,1.],repeat=n):
        H=np.zeros((b,n)); H[list(hashes),range(n)]=signs
        Z=H.T@H-np.eye(n)
        W=X@Z@X
        G=X@H.T
        P=G@np.linalg.pinv(G)
        R=X@(np.eye(n)-P)@X
        lhs=np.sum(R*R); rhs=np.trace(A@Z@A@Z)
        assert lhs <= rhs+1e-10
        second.append(rhs)
        diagonal.append(np.sum(np.diag(A@Z)**2))
        trace.append(np.trace(A@Z)**2)
        drift.append(np.trace(A-R))
        convex.append(np.trace(A@A)-np.trace(R@R))
        mixed.append(np.trace(H@B@H.T)*np.trace(H@A@H.T))
off=np.sum(A*A)-np.sum(np.diag(A)**2)
expected=(np.trace(A)**2+np.sum(A*A)-2*np.sum(np.diag(A)**2))/b
assert abs(np.mean(second)-expected)<1e-12
assert abs(np.mean(diagonal)-off/b)<1e-12
assert abs(np.mean(trace)-2*off/b)<1e-12
mixed_expected=np.trace(B)*np.trace(A)+2/b*(np.trace(B@A)-np.dot(np.diag(B),np.diag(A)))
assert abs(np.mean(mixed)-mixed_expected)<1e-12
kappa=1+2/b
assert np.mean(drift)>=np.trace(A@A)/(kappa*np.trace(A))
assert np.mean(convex)>=np.trace(A@A@A)/(kappa*np.trace(A))

# Full-rank RP set laws, including complex isometries.
rng=np.random.default_rng(20260930)
max_gate_ratio=0.; max_flat_error=0.
for N,k in [(4,2),(5,2),(5,3),(6,3)]:
    for _ in range(6):
        V=np.linalg.qr(rng.standard_normal((N,k))+1j*rng.standard_normal((N,k)))[0]
        T=rng.standard_normal((k,k))+1j*rng.standard_normal((k,k))
        for C in [np.eye(k),T.conj().T@T+.05*np.eye(k)]:
            M=V@C@V.conj().T
            laws={J:0. for J in combinations(range(N),k)}
            for order in permutations(range(N),k):
                R=M.copy(); p=1.
                for i in order:
                    d=float(R[i,i].real)
                    p*=max(d,0)/float(np.trace(R).real)
                    if d>1e-14:
                        R=R-np.outer(R[:,i],R[i,:])/d
                        R=(R+R.conj().T)/2
                    else:
                        p=0.; break
                laws[tuple(sorted(order))]+=p
            lamb=np.linalg.eigvalsh(C)[::-1]
            price=factorial(k)*float(np.linalg.det(C).real)/prod(sum(lamb[j:]) for j in range(k))
            assert abs(sum(laws.values())-1)<1e-8
            for J,p in laws.items():
                nu=abs(np.linalg.det(V[list(J)]))**2
                max_gate_ratio=max(max_gate_ratio,p/(price*nu))
                assert p <= price*nu+1e-10
                if np.allclose(C,np.eye(k)):
                    max_flat_error=max(max_flat_error,abs(p-nu))

# Stable tridiagonal propagation checks numerical headroom without
# interpreting a finite grid as an all-order theorem.
def majorant_root(d,q,eps,C=69,B=8):
    s=ceil(B*(2*q+1)/eps); b=ceil(C*(d+2*q+1)/(s*eps**2)); m=s*b; h=b-1
    nu=2*np.arange(q+1,dtype=float)
    alpha=(np.sqrt((m+nu)*(d+nu+1))+np.sqrt(h*nu*(d+nu-1)/2)+nu*np.sqrt(h)+np.sqrt(b*h*nu*(nu-1).clip(0)/2))/m
    alpha[0]=np.sqrt((d+1)/m)
    chi=(d+2*nu-1+2*np.sqrt(b*(nu-1).clip(0)*(d+nu-1))+(nu+1)*np.sqrt(b)+h*nu)/m
    chi[0]=0
    # Log domain prevents the tiny vacuum component from underflowing at q=3000.
    with np.errstate(divide='ignore'):
        la=np.log(alpha[:-1]); lc=np.log(chi)
    x=np.full(q+1,-np.inf); x[0]=0.
    for _ in range(2*q):
        y=lc+x
        y[1:]=np.logaddexp(y[1:],la+x[:-1])
        y[:-1]=np.logaddexp(y[:-1],la+x[1:])
        x=y
    return np.exp(x[0]/(2*q))/eps
majorant={}
for C in [81,69,40,35,30]:
    B=8.5 if C==81 else 8
    vals=[(majorant_root(d,q,eps,C,B),d,q,eps) for d,q,eps in product([1,10,100,10**4,10**6,10**8],[2,5,10,40,100,200],[1.,.5,.05])]
    if C in [69,35]:
        vals +=[(majorant_root(d,q,1.,C,B),d,q,1.) for d,q in product([10,10**3,10**6,10**9],[200,1000,3000])]
    worst=max(vals)
    majorant[str(C)]={'count':len(vals),'max_root':worst[0],'at':list(worst[1:])}

def renyi_cost(k,beta,delta=.01,alpha=None):
    lf=lgamma(k+1); bk=log(1+2*lf/k)
    alpha=min(beta,(k-1)/(k*bk)) if alpha is None else alpha
    bound=k*bk*alpha-(k-1)*log(alpha)
    p=1/beta
    return (p*log(1/delta)+p*bound+k*log(2)+p*(k*log(k)-lf))/(p-1)
k=10
lf=lgamma(k+1); bk=log(1+2*lf/k); beta=2/3
prices=[0.,log(252),lf,lgamma(1001)-lgamma(11)-lgamma(991),
        renyi_cost(k,beta,alpha=beta/(1+beta*bk))-log(100),
        renyi_cost(k,beta)-log(100),
        min(renyi_cost(k,bb) for bb in np.linspace(.01,.99,981))-log(100)]
q_from_prices=[max(1,ceil((log(100/.01)+p)/log(4))) for p in prices]
assert q_from_prices==[7,11,18,46,72,69,32],q_from_prices

result={'certificate_count':len(cert),'certificate_max_T':max(x['T'] for x in cert),'certificate_table_matches':True,'certificate':cert,'gate_q_from_prices':q_from_prices,'gate_table':table,'exact_enumeration':{'number_of_sketches':len(second),'weighted_second_moment':float(np.mean(second)),'formula':float(expected),'diagonal_MSE':float(np.mean(diagonal)),'trace_variance':float(np.mean(trace)),'mean_trace_drop':float(np.mean(drift)),'drift_bound':float(np.trace(A@A)/(kappa*np.trace(A))),'mean_square_potential_drop':float(np.mean(convex)),'square_potential_bound':float(np.trace(A@A@A)/(kappa*np.trace(A)))},'gate_instances':48,'max_gate_ratio':max_gate_ratio,'max_flat_gate_error':max_flat_error,'majorant_grid':majorant}
(ROOT/'independent_audit_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['certificate','gate_table']},indent=2))
