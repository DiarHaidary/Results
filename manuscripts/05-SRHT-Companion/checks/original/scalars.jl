# Scalar facts used in the proofs (diagnostic re-checks)
setprecision(BigFloat, 256)
println("80<81: 2/sqrt20+1/20 < 1/2 : ", 2/sqrt(big(20)) + big(1)/20, " < 0.5")
println("288<289: sqrt2/3+1/36 = ", sqrt(big(2))/3 + big(1)/36)
println("646 retune: 12/sqrt646+18/646 = ", 12/sqrt(big(646)) + big(18)/646)
println("sqrt(pi)/24 = ", sqrt(big(pi))/24, "  vs 1/12 = ", big(1)/12)
println("kappa0 = e^{2-gamma}/2 = ", exp(2 - big(Base.MathConstants.eulergamma))/2)
function bgamma(x::BigFloat); z = BigFloat(); ccall((:mpfr_gamma, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, 0); z; end
g(p) = 2^p*bgamma(p + big"0.5")/sqrt(big(pi))
for e in (big"1e-2", big"1e-4", big"1e-8"); println("g(1+e)^(1/e) at e=$e: ", g(1+e)^(1/e)); end
# atom lemma: P{Bin(n,rho)=M} <= min(1/2, sqrt(pi/(8v))), v = M(1-M/n)
worst = 0.0; wh = (0,0)
for n in 2:400, M in 1:n-1
    rho = M/n; v = M*(1-rho)
    p = Float64(binomial(big(n), M) * big(rho)^M * (1-big(rho))^(n-M))
    b = min(0.5, sqrt(pi/(8v)))
    if p/b > worst; global worst = p/b; global wh = (n, M); end
end
println("atom lemma: max p/bound over 2<=n<=400 = ", worst, " at ", wh)
# p_{n,M} = P{Bin(n-1,rho)=M} = P{Bin(n,rho)=M}
n, M = 37, 11; rho = big(M)/n
println("atom identity diff: ", binomial(big(n-1),M)*rho^M*(1-rho)^(n-1-M) - binomial(big(n),M)*rho^M*(1-rho)^(n-M))
# coupling constant a = E(K-M)_+ = v P(Bin(n-1,rho)=M)
for (n, M) in ((10,3), (40,17), (200,5))
    rho = big(M)/n
    a = sum(max(k-M,0)*binomial(big(n),k)*rho^k*(1-rho)^(n-k) for k in 0:n)
    v = M*(1-rho); pY = binomial(big(n-1),M)*rho^M*(1-rho)^(n-1-M)
    println("coupling (n,M)=($n,$M): a - v p = ", Float64(a - v*pY))
end
# D <= Lambda <= 3D on a grid of coherences
mn = Inf; mx = 0.0
for al in (1.0, 1.5, 2.0, 4.0, 10.0), be in (1.0, 2.0, 5.0), d in (1, 10, 1000), q in 1:30
    a = al*q*(2q+1); b = be*(d+8q^2) + al*q*(2q-1); D = be*d + (8be+4al)*q^2
    @assert abs(a + b - D) < 1e-9*D
    L = (sqrt(a) + sqrt(2b))^2
    global mn = min(mn, L/D); global mx = max(mx, L/D)
end
println("Lambda/D range on grid: [", mn, ", ", mx, "]")
# rank-two comparison delta/2 > 8 delta^2/9 iff delta < 9/16
println("delta/2 - 8delta^2/9 at 9/16: ", big(9)/32 - 8*(big(9)/16)^2/9)
# q_*^2 <= 16 d at delta = 0.01
ok = all(d -> (q = ceil(Int, log(big(120d))/log(big(4))); q^2 <= 16d), 1:10^6)
println("q*^2 <= 16d for d <= 1e6 at delta=0.01: ", ok)
# q_orig = q_f + 1 always, q_* in {q_f, q_f+1}
chk = true
for d in 1:5000, del in (0.3, 0.1, 0.01, 1e-3, 1e-6)
    lg(x) = log(big(x))/log(big(4))
    qf = ceil(Int, lg(d/big(del))); qo = ceil(Int, lg(4d/big(del))); qs = ceil(Int, lg(6d/(5big(del))))
    global chk &= (qo == qf + 1) && (qs in (qf, qf+1))
end
println("q_orig = q_f+1 and q_* in {q_f,q_f+1}: ", chk)
