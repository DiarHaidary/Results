from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

e_lower=sum((F(1,factorial(k)) for k in range(21)),F(0))
e_upper=e_lower+F(22,21*factorial(21))
target=F(1,100*64)

def moment(m,r,q,eigs,nu,e_bound):
    b=e_bound*(27*r*q+192*q*q+24*q)/m
    total=F(0)
    for eigen in eigs:
        w=[eigen/m]+[F((k+1)*(m-k),m*m)*nu for k in range(1,q)]
        state=[F(1)]+[F(0)]*q
        for _ in range(2*q):
            nxt=[F(0)]*(q+1)
            for k,x in enumerate(state):
                if k: nxt[k-1]+=x
                if k<q:nxt[k+1]+=x*w[k]
                if k:nxt[k]+=x*b
            state=nxt
        total+=state[0]
    return total

data=[]
for r,q,vals,label in [(4,3,[F(15)]*4,'full 2x2 product'),
                      (10,3,[F(35)]*10,'rank10 common factor'),
                      (2,3,[F(133,12),F(49,6)],'anisotropic rank2 source example')]:
    old=F(4*r+4)
    inputs=[([old]*r,old),(vals,old),([max(vals)]*r,max(vals)),(vals,max(vals))]
    counts={'full 2x2 product':[14514,13668,13606,13606],
            'rank10 common factor':[24869,23649,23447,23447],
            'anisotropic rank2 source example':[10581,10115,10403,10107]}[label]
    for case,((eigs,nu),m) in enumerate(zip(inputs,counts)):
        current_upper=moment(m,r,q,eigs,nu,e_upper)
        previous_lower=moment(m-1,r,q,eigs,nu,e_lower)
        assert current_upper<=target
        assert previous_lower>target
        data.append(dict(model=label,case=case,m=m,
            positive_margin=float(target-current_upper),
            previous_excess=float(previous_lower-target)))
Path(__file__).with_name('kr_spectrum_results.json').write_text(json.dumps(dict(
    status='passed',method='Exact rational positive-path recurrence with Taylor bounds on e',
    cases=data),indent=2),encoding='utf8')
print('Passed all12 exact rational candidate/previous-integer checks.')
