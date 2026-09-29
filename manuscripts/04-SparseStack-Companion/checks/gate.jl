# RP-gate parameter tables and Renyi comparison (exact integer arithmetic where it matters).
using LinearAlgebra, Random
setprecision(256)
d = 100; ε = big(1)//2; δ = big(1)//100
qof(logx) = max(1, ceil(Int, logx / log(big(4)) - big(1)//10^30))   # ceil(log4(x)), x = exp(logx)
function qexact(x)                  # smallest q>=1 with 4^q >= x, x rational/integer
    q = 1; while big(4)^q < x; q += 1; end; q
end
function params(q, C, B)
    v = 2q + 1; D = d + v
    s = ceil(Int, B * v / ε)
    b = ceil(Int, C * D / (s * ε^2))
    (q, s, b, s * b)
end
# Renyi-3/2 budgets
β = big(2)//3
bk(k) = log(1 + 2 * log(factorial(big(k))) / k)
function Bnote(k); b = bk(k); k * β * b / (1 + β * b) + (k - 1) * log(b + 1 / β); end
function Bra01(k); b = bk(k); α = min(β, (k - 1) / (k * b)); k * b * α - (k - 1) * log(α); end
frakd(k, Bf) = 3 * Bf(k) + 2k * log(big(2)) + 3 * (k * log(big(k)) - log(factorial(big(k))))
k = 10
println("log 10! = ", Float64(log(factorial(big(10)))))
println("frak d (note form) = ", Float64(frakd(k, Bnote)), "   (RA-01 form) = ", Float64(frakd(k, Bra01)))
for (C, B, name) in [(81, big(17)//2, "v3.0 (81,17/2)"), (69, big(8), "new (69,8)")]
    println("== ", name)
    for (lab, Λ) in [("independent", big(1)), ("5+5 levels", div(factorial(big(10)), factorial(big(5))^2)),
                     ("10!", factorial(big(10))), ("union binom(1000,10)", binomial(big(1000), 10))]
        println(rpad(lab, 24), params(qexact(d * Λ / δ), C, B))
    end
    for (lab, Bf) in [("Renyi note-form", Bnote), ("Renyi RA-01 form", Bra01)]
        logx = log(big(d)) + frakd(k, Bf) - 3 * log(δ)       # d / delta0 with delta0 = delta^3 e^{-frak d}
        println(rpad(lab, 24), params(qof(logx), C, B))
    end
end
# crossover: Renyi beats factorial iff frak d_k + 2 log(1/delta) < log k!
for (lab, Bf) in [("note", Bnote), ("RA-01", Bra01)]
    kk = 2
    while frakd(kk, Bf) + 2 * log(1 / δ) >= log(factorial(big(kk))); kk = kk < 1000 ? kk + 1 : Int(ceil(kk * 1.01)); end
    println("Renyi($lab) first beats k! (delta=0.01) near k = ", kk)
    kk = 2
    while frakd(kk, Bf) >= log(factorial(big(kk))); kk = kk < 1000 ? kk + 1 : Int(ceil(kk * 1.01)); end
    println("   ignoring the delta^3 cost, near k = ", kk)
end
# Lambda_C for spike
lam(λ) = (kk = length(λ); λs = sort(λ, rev = true); factorial(kk) * prod(λs) / prod(sum(λs[j+1:end]) for j in 0:kk-1))
η = 3.0; println("spike: ", lam([1 + η; ones(9)]), " vs ", 10 * (1 + η) / (10 + η))
# brute-force check of mu_C(J) <= Lambda_C nu_V(J)
function rpsetlaw(A, k)
    N = size(A, 1); law = Dict{Vector{Int},Float64}()
    function rec(R, path, p)
        if length(path) == k
            key = sort(path); law[key] = get(law, key, 0.0) + p; return
        end
        t = real(tr(R))
        for i in 1:N
            Rii = real(R[i, i]); Rii <= 1e-14 * t && continue
            Rn = R - R[:, i] * R[i:i, :] / Rii
            rec(Rn, [path; i], p * Rii / t)
        end
    end
    rec(A, Int[], 1.0); law
end
Random.seed!(7); worst = 0.0
for trial in 1:40
    N = rand(4:6); kk = rand(2:3)
    V = Matrix(qr(randn(ComplexF64, N, kk)).Q)[:, 1:kk]
    G = randn(ComplexF64, kk, kk); Cm = G * G' + 0.05I; Cm = (Cm + Cm') / 2
    law = rpsetlaw(V * Cm * V', kk); ΛC = lam(real(eigvals(Hermitian(Cm))))
    for (J, μ) in law
        ν = abs2(det(V[J, :])); global worst = max(worst, μ / (ΛC * ν))
    end
end
println("max over random cases of mu/(Lambda_C nu) = ", worst)
