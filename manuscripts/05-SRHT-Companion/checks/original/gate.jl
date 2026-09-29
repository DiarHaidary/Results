# RP-gate and rank-two tables (diagnostic recomputation); reuses prescriptions.jl definitions
include_string(Main, replace(read("prescriptions.jl", String), r"println\(\"kappa0.*"s => ""))
d = 100; k = 10; N = 1000
prices = [("independent target", big(1)), ("two flat bands 5+5", big(factorial(10) ÷ (factorial(5)^2))),
          ("arbitrary spectrum 10!", big(factorial(10))), ("union bound binom(1000,10)", binomial(big(1000), 10))]
println("GATE d=$d eps=1/2 delta=0.01 k=$k N=$N")
for (name, G) in prices
    del0 = delta / G
    qs = qstar(d, del0); M = r(20Lam(qs, d)); M60 = r(60Dq(qs, d))
    Mmin, qd = direct(d, del0; qmax=80)
    println(name, " | Gamma=", G, " | q*=", qs, " | 20Lam=", M, " | 60D=", M60, " | direct=", Mmin+10, " q=", qd, " cert=", Float64(cert(n, Mmin+10, qd, d, eps)/del0))
end
println("\nRANK TWO: delta0 = delta/2 vs 8 delta^2/9")
for d in (10, 100, 1000, 10000)
    q1 = qstar(d, delta/2); q2 = qstar(d, 8delta^2/9)
    println(d, " | q*(d/2)=", q1, " M=", r(20Lam(q1,d)), " | q*(8d^2/9)=", q2, " M=", r(20Lam(q2,d)))
end
