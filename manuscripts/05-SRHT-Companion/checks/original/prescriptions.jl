# Prescription tables for the two-round SRHT companion (diagnostic recomputation).
setprecision(BigFloat, 256)
function logfact(k::Integer)
    k < 1000 && return sum(log(BigFloat(i)) for i in 1:k; init=BigFloat(0))
    x = BigFloat(k) + 1
    s = (x - big"0.5")*log(x) - x + log(2*BigFloat(pi))/2
    s += 1/(12x) - 1/(360x^3) + 1/(1260x^5) - 1/(1680x^7) + 1/(1188x^9)
    return s
end
logbinom(a, b) = logfact(a) - logfact(b) - logfact(a - b)
function patom(n, M)          # P{Bin(n-1,M/n)=M}
    rho = BigFloat(M)/n
    exp(logbinom(n-1, M) + M*log(rho) + (n-1-M)*log1p(-rho))
end
aq(q; al=1) = al*q*(2q+1)
bq(q, d; al=1, be=1) = be*(d + 8q^2) + al*q*(2q-1)
Dq(q, d; al=1, be=1) = be*d + (8be + 4al)*q^2
Lam(q, d) = (sqrt(BigFloat(aq(q))) + sqrt(2*BigFloat(bq(q, d))))^2
function root(n, M, q, d)
    rho = BigFloat(M)/n; L = Lam(q, d)
    (2*sqrt((1-rho)*L/M) + abs(1-2rho)*L/M) / (1 - patom(n, M))
end
cert(n, M, q, d, eps) = d * (root(n,M,q,d)/eps)^(2q)
n = 2^30; eps = big"0.5"; delta = big"0.01"
lg4(x) = log(BigFloat(x))/log(BigFloat(4))
qorig(d, del) = ceil(Int, lg4(4d/del))
qstar(d, del) = ceil(Int, lg4(6d/(5del)))
qfree(d, del) = ceil(Int, lg4(d/del))
r(x) = ceil(Int, x/eps^2)
function direct(d, del; qmax=40)
    best = (typemax(Int), 0)
    for q in 1:qmax
        cert(n, n-1, q, d, eps) > del && continue
        lo, hi = 1, n-1
        while lo < hi
            mid = (lo + hi) >>> 1
            if cert(n, mid, q, d, eps) <= del; hi = mid; else lo = mid + 1; end
        end
        if lo < best[1]; best = (lo, q); end
    end
    return best
end
println("kappa0 = ", exp(2 - Base.MathConstants.eulergamma)/2)
println("\nMAIN TABLE  d | q q* | Dq* Lam_q* | 216D_q | 72Lam_q | 60D_q* | 20Lam_q* | direct(min) q | +10 cert | indep")
for d in (10, 100, 1000, 10000)
    q = qorig(d, delta); qs = qstar(d, delta)
    Mmin, qd = direct(d, delta)
    Mf = Mmin + 10
    println(d, " | ", q, " ", qs, " | ", Dq(qs,d), " ", Float64(Lam(qs,d)), " | ", r(216Dq(q,d)), " | ", r(72Lam(q,d)), " | ",
      r(60Dq(qs,d)), " | ", r(20Lam(qs,d)), " | ", Mmin, " ", qd, " | ", Mf, " ", Float64(cert(n,Mf,qd,d,eps)),
      " (at min: ", Float64(cert(n,Mmin,qd,d,eps)), ", at min-1: ", Float64(cert(n,Mmin-1,qd,d,eps)), ") | ", min(n,4qd))
end
println("\nexact-root failure bounds at canonical 20Lam_q*:")
for d in (10, 100, 1000, 10000)
    qs = qstar(d, delta); M = r(20Lam(qs,d)); println(d, "  ", Float64(cert(n,M,qs,d,eps)), "  p_atom=", Float64(patom(n,M)))
end
println("\nq-slack: 72Lam at q=ceil(log4(d/delta)) vs q orig")
for d in (1, 10, 100, 1000, 10000)
    q = qorig(d, delta); qf = qfree(d, delta)
    a = r(72Lam(q,d)); b = r(72Lam(qf,d)); c = r(216Dq(qf,d))
    println(d, " q=", q, " qf=", qf, " 72Lam(q)=", a, " 72Lam(qf)=", b, " saving=", round(100*(1-b/a), digits=1), "%  216D(qf)=", c, " 216D(q)=", r(216Dq(q,d)), " 20Lam(q*)=", r(20Lam(qstar(d,delta),d)), " q*=", qstar(d,delta))
end
println("\nPLOT DATA (d, 216D_q, 72Lam_q, 60D_q*, 20Lam_q*, direct+10, qdirect)")
for d in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
    q = qorig(d, delta); qs = qstar(d, delta); Mmin, qd = direct(d, delta)
    println("(", d, ", ", r(216Dq(q,d)), ", ", r(72Lam(q,d)), ", ", r(60Dq(qs,d)), ", ", r(20Lam(qs,d)), ", ", Mmin+10, ", ", qd, ")")
end
# sup over d of 20Lam_{q*}/d at delta=0.01, and of 60D/d
mx = 0.0; arg = 0; mx2 = 0.0
for d in 1:300000
    v = Float64(20Lam(qstar(d,delta), d)/d); w = Float64(60Dq(qstar(d,delta), d)/d)
    if v > mx; global mx = v; global arg = d; end
    if w > mx2; global mx2 = w; end
end
println("\nsup_{d<=3e5} 20Lam_q*/d = $mx at d=$arg ;  sup 60D_q*/d = $mx2 (60*193 = ", 60*193, ")")
