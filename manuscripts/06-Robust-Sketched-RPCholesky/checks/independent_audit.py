"""Independent numerical checks for the robust RPCholesky manuscript.

Run: python checks/independent_audit.py
These finite numerical checks are not a formal proof or a verification of
the imported initial-moment and sketch-certificate theorems.
"""
from pathlib import Path
from itertools import permutations
import json
import math
import numpy as np

RNG = np.random.default_rng(29092026)


def residual(A, J):
    if not J:
        return A.copy()
    return A - A[:, J] @ np.linalg.solve(A[np.ix_(J, J)], A[J, :])


def trace(A):
    return float(np.trace(A).real)


def qr(F):
    return np.linalg.qr(F, mode="reduced")[0]


def path_probability(A, w):
    R = A.copy()
    probability = 1.0
    for i in w:
        probability *= R[i, i].real / trace(R)
        R -= np.outer(R[:, i], R[i, :]) / R[i, i]
    return probability


def h4(r, beta):
    lf = math.lgamma(r + 1)
    b = math.log1p(2 * lf / r)
    alpha = min(beta, (r - 1) / (r * b))
    B = r * b * alpha - (r - 1) * math.log(alpha)
    golden = (1 + math.sqrt(5)) / 2
    cv = 1 + 2 * math.log(golden) + 1 / golden
    return B + r * math.log(r) - lf + beta * r * math.log(2) + (1-beta) * cv * r


def progress(u):
    return 1/u - math.log(u)


def count(r, eps, theta, H, a):
    c = math.log(2*math.e/a)
    g = math.ceil(math.log(max(1/theta, (H/(r*theta)-c)/math.log(5), 1)) / -math.log1p(-1/r))
    rho = (1-1/r)**g
    D = (1-rho)*r*c + rho*H/theta
    s = math.ceil(D/math.log(5))
    n = math.ceil((r/a)*(progress(eps)-progress(10/a-1+1/r)))
    return r+g+s+n


def gamma(eta):
    return (1-eta-eta*eta)/(1+eta)


def chi2(eta):
    return (eta/math.sqrt(2)+eta*eta/(2*(1-eta)))**2/(1-eta-eta*eta)**2


def best_count(r, eps, a=1, eta=None):
    values = []
    for beta in np.arange(0.001, 1.0, 0.001):
        H = h4(r, beta)
        theta = 1-beta
        if eta is None:
            H += r*math.log(1/a)
        else:
            H = (H+r*chi2(eta))/2
            theta /= 2
            a = gamma(eta)
        values.append(count(r, eps, theta, H, a))
    return min(values)


results = {"method": "independent NumPy implementation; deterministic seed 29092026", "checks": {}}
maximum_identity_error = 0.0
good_kernel_cases = 0
for complex_data in (False, True):
    for trial in range(40):
        D, N = 9, 6
        F = RNG.normal(size=(D, N))
        if complex_data:
            F = F + 1j*RNG.normal(size=(D, N))
        A = F.conj().T @ F
        Pi = RNG.normal(size=(11, D))/math.sqrt(11)
        B = (Pi @ F).conj().T @ (Pi @ F)
        w = tuple(RNG.choice(N, size=3, replace=False))
        Q = qr(F[:, w])
        G = Q.conj().T @ Pi.T @ Pi @ Q
        product = np.linalg.det(G).real
        for s in range(len(w)):
            product *= trace(residual(A, list(w[:s])))/trace(residual(B, list(w[:s])))
        ratio = path_probability(B, w)/path_probability(A, w)
        maximum_identity_error = max(maximum_identity_error, abs(ratio-product)/max(1,abs(product)))
        R = residual(A, list(w))
        E = F-Q@(Q.conj().T@F)
        P = Pi@Q@np.linalg.solve(G,Q.conj().T@Pi.T)
        RB = residual(B, list(w))
        err = np.linalg.norm(RB-(Pi@E).conj().T@(np.eye(11)-P)@(Pi@E))/max(1,np.linalg.norm(RB))
        maximum_identity_error = max(maximum_identity_error, float(err))
        XB = np.linalg.lstsq(Pi@F[:,w],Pi@F,rcond=None)[0]
        measured = np.linalg.norm(F-F[:,w]@XB)**2
        formula = trace(R)+np.linalg.norm(np.linalg.solve(G,Q.conj().T@Pi.T@Pi@E))**2
        maximum_identity_error = max(maximum_identity_error, abs(measured-formula)/max(1,abs(formula)))

        # Exact drift identities for the determinant and the trace.
        r = 2
        tau = sum(np.linalg.eigvalsh(A)[:-r])
        t = r/tau
        current = residual(A, list(w[:1]))
        T = trace(current)
        eigen = np.linalg.eigvalsh(current)
        pot = float(np.prod(1+t*eigen))
        expected_pot = expected_trace = 0.0
        for i in range(N):
            if current[i,i].real < 1e-10:
                continue
            nxt = current-np.outer(current[:,i],current[i,:])/current[i,i]
            weight = current[i,i].real/T
            expected_pot += weight*np.linalg.det(np.eye(N)+t*nxt).real
            expected_trace += weight*trace(nxt)
        predicted_pot = pot*sum(eigen/(1+t*eigen))/T
        predicted_trace = T-trace(current@current)/T
        assert np.isclose(expected_pot,predicted_pot,rtol=1e-10,atol=1e-10)
        assert np.isclose(expected_trace,predicted_trace,rtol=1e-10,atol=1e-10)

        # All three weighted tests, functional comparisons, and chi-square.
        eta = 0.25
        H = RNG.normal(size=(D,D)); H=(H+H.T)/2
        H *= (eta/2)/np.linalg.norm(H,2)
        Pi = np.linalg.cholesky(np.eye(D)+H).T
        J = list(w[:2]); Q=qr(F[:,J]); E=F-Q@(Q.conj().T@F)
        R=residual(A,J); RB=residual((Pi@F).conj().T@(Pi@F),J)
        active=np.flatnonzero(np.diag(R).real>1e-9)
        diagonal=np.diag(R).real[active]; T=trace(R)
        pi=diagonal/T
        decrease=np.array([np.vdot(R[:,i],R[:,i]).real/R[i,i].real for i in active])
        potential=np.array([np.linalg.det(np.eye(N)+t*(R-np.outer(R[:,i],R[i,:])/R[i,i])).real for i in active])
        weights=[pi,pi*decrease/(pi@decrease),pi*potential/(pi@potential)]
        norms=[]
        for i in active:
            u=E[:,i]/np.linalg.norm(E[:,i]); U=np.column_stack([Q,u])
            norms.append(np.linalg.norm(U.conj().T@Pi.T@Pi@U-np.eye(U.shape[1]),2))
        norms=np.array(norms)
        assert all(omega@(norms**4)<=eta**4/4+1e-12 for omega in weights)
        y=np.diag(RB).real[active]/diagonal
        ybar=pi@y
        q=pi*y/ybar
        assert q@potential <= (pi@potential)/gamma(eta)+1e-10
        assert q@decrease >= gamma(eta)*(pi@decrease)-1e-10
        assert sum(q*q/pi)-1 <= chi2(eta)+1e-10
        good_kernel_cases += 1

assert maximum_identity_error < 1e-10
results["checks"]["algebraic_identity_cases"] = 80
results["checks"]["maximum_relative_identity_error"] = maximum_identity_error
results["checks"]["good_kernel_cases"] = good_kernel_cases

expected_rows = [(10,1,[71,85,123,235]),(10,.1,[185,226,349,687]),
 (100,1,[722,848,1237,2360]),(100,.1,[1852,2261,3497,6881]),
 (1000,1,[7243,8497,12393,23610]),(1000,.1,[18545,22625,34998,68821])]
for r,eps,expected in expected_rows:
    calculated=[best_count(r,eps,a) for a in (1,.8,.5,.25)]
    assert calculated==expected,(r,eps,calculated)
results["checks"]["robust_count_table"] = "all 24 entries reproduced"

rows=[]
for r,eps,eta in [(10,1,.1),(10,1,.25),(100,1,.1),(100,1,.25),
                  (100,.1,.1),(100,.1,.25),(1000,.1,.1),(1000,.1,.25)]:
    K=best_count(r,eps/2,eta=eta)
    ell=max(2,math.ceil(2*math.log2(4*(K+1)*1e6/eps)+K*chi2(eta)/math.log(2)))
    sparsity=math.ceil(8.5*(2*ell+1)/eta)
    b=math.ceil(81*(K+1+2*ell+1)/(eta*eta*sparsity))
    rows.append({"r":r,"epsilon":eps,"eta":eta,"K":K,"ell":ell,
                 "SparseStack_rows":sparsity*b,"SRHT_uncapped_threshold":math.ceil(216*(K+1+12*ell*ell)/(eta*eta))})
results["checks"]["sketch_parameter_table"] = rows
assert all(chi2(x)<=2*x*x for x in np.linspace(.001,.25,1000))

# Exact finite enumeration of gate probabilities and union amplification.
gate_cases=0
union_cases=0
for case in range(5):
    N,k=6,3
    V=qr(RNG.normal(size=(N,k)))
    C0=RNG.normal(size=(k,k)); C=C0.T@C0+np.eye(k)
    M=V@C@V.T; values=np.linalg.eigvalsh(C)[::-1]
    Lambda=math.factorial(k)*np.linalg.det(C)/np.prod([sum(values[j:]) for j in range(k)])
    set_prob={}
    for w in permutations(range(N),k):
        key=tuple(sorted(w)); set_prob[key]=set_prob.get(key,0)+path_probability(M,w)
    assert abs(sum(set_prob.values())-1)<1e-10
    assert 1-1e-10<=Lambda<=math.factorial(k)+1e-10
    assert all(p<=Lambda*np.linalg.det(V[list(J),:])**2+1e-10 for J,p in set_prob.items())
    gate_cases+=1
    F=RNG.normal(size=(8,5)); A=F.T@F; k=2
    set_prob={}
    for w in permutations(range(5),k):
        key=tuple(sorted(w)); set_prob[key]=set_prob.get(key,0)+path_probability(A,w)
    fractional=sum(p*math.sqrt(trace(residual(A,list(J)))) for J,p in set_prob.items())
    union=sum(p*q*trace(residual(A,sorted(set(J)|set(L))))
              for J,p in set_prob.items() for L,q in set_prob.items())
    assert union<=fractional*fractional+1e-10
    union_cases+=1
results["checks"]["gate_exact_enumeration_cases"]=gate_cases
results["checks"]["union_exact_enumeration_cases"]=union_cases

def renyi_price(k,delta):
    p=np.linspace(1.01,30,3000)
    lf=math.lgamma(k+1); b=math.log1p(2*lf/k)
    alpha=np.minimum(1/p,(k-1)/(k*b))
    B=k*b*alpha-(k-1)*np.log(alpha)
    values=p/(p-1)*(math.log(1/delta)+B+(k/p)*math.log(2)+k*math.log(k)-lf)
    return float(min(values))

crossover={}
for delta in (.01,1e-6):
    first=next(k for k in range(2,200) if renyi_price(k,delta)<math.log(1/delta)+math.lgamma(k+1))
    crossover[str(delta)]=first
assert crossover=={"0.01":128,"1e-06":129}
results["checks"]["Renyi_crossover"]=crossover
gate_prices=[math.log(100/.01),math.log(100/.01)+math.log(100/19),
             math.log(100/.01)+math.lgamma(11),math.log(100)+renyi_price(10,.01),
             math.log(100/.01)+math.lgamma(1001)-math.lgamma(11)-math.lgamma(991)]
gate_rows=[]
for price in gate_prices:
    ell=math.ceil(price/math.log(4)); sparsity=math.ceil(8.5*(2*ell+1)/.5)
    b=math.ceil(81*(100+2*ell+1)/(.25*sparsity))
    gate_rows.append([ell,sparsity,b,sparsity*b])
assert gate_rows==[[7,255,147,37485],[8,289,132,38148],[18,629,71,44659],
                   [32,1105,49,54145],[46,1581,40,63240]]
results["checks"]["gate_parameter_table"]=gate_rows

floor_cases=0
for a in (.25,.5,.8,1):
    for r in (1,2,7):
        for k in (0,1,9):
            delta=.01; m=max(k+r+1,math.ceil(k/(a*delta)))
            head=(m-k)*(1-a)/(a*r)
            for s in range(k):
                tails=m-s
                assert 1/tails<=1/(a*(r*head+tails))+1e-14
            eig=sorted([head]*r+[1.0]*m,reverse=True)
            tau=sum(eig[r:])
            assert tau>0 and ((m-k)/a)/tau>=1/a-delta-1e-10
            floor_cases+=1
results["checks"]["floor_edge_cases"]=floor_cases
results["status"]="passed"
out=Path(__file__).with_name("independent_audit_results.json")
out.write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
print(json.dumps(results,indent=2))
