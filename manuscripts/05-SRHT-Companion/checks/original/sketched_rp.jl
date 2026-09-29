# Illustrative SRHT parameters for the path-local sketched-RP theorem (flat transforms)
setprecision(BigFloat, 256)
aq(q) = q*(2q+1); bq(q,d) = (d + 8q^2) + q*(2q-1)
Lam(q,d) = (sqrt(BigFloat(aq(q))) + sqrt(2*BigFloat(bq(q,d))))^2
lg4(x) = log(BigFloat(x))/log(BigFloat(4))
for (r, eps, ratio, eta, N) in ((5, 0.5, 1e3, 0.25, 10^4), (5, 0.5, 1e3, 0.1, 10^4))
    g = (1-eta)/(1+eta)
    ell = 1 + ceil(Int, (r/g)*(1/eps + log(1/eps) + log(ratio)))
    del0 = eps / (N*ell*BigFloat(g)^(-(ell-1))*ratio)   # additive term <= eps*tau
    q = ceil(Int, lg4(6ell/(5del0)))
    M = ceil(Int, 20Lam(q, ell)/eta^2)
    println("r=$r eps=$eps trK/tau=$ratio eta=$eta N=$N gamma=$g ell=$ell log10 delta0=", Float64(log10(del0)), " q=$q M=$M  M/ell=", M/ell)
end
