# Constants for the manuscript (final formulas). Functional conversion calK_a(r,eps;theta,H);
# per-coordinate two-sided: H -> H + r log(1/a);  sketched: theta=(1-beta)/2, H=(H4+r chi^2)/2, a=gamma.
using Printf
const φg = (1 + sqrt(5.0)) / 2
const cvol = 1 + 2log(φg) + 1 / φg
logfact(k) = k == 0 ? 0.0 : sum(log.(1:k))
br(r) = log(1 + 2logfact(r) / r)
Bdef(r, β) = r == 1 ? 0.0 : (b = br(r); α = min(β, (r - 1) / (r * b)); r * b * α - (r - 1) * log(α))
loga(r) = r * log(r) - logfact(r)
H4(r, β) = Bdef(r, β) + loga(r) + β * r * log(2) + (1 - β) * cvol * r
F(u) = 1 / u - log(u)
function calK(r, ε, θ, H; a=1.0)
    ca = log(2 * exp(1.0) / a)
    arg = max(1 / θ, (H / (r * θ) - ca) / log(5), 1.0)
    g = ceil(Int, log(arg) / (-log(1 - 1 / r)))
    ρ = (1 - 1 / r)^g
    D = (1 - ρ) * r * ca + ρ * H / θ
    s = ceil(Int, D / log(5))
    ua = 10 / a - 1 + 1 / r
    n = ceil(Int, (r / a) * (F(ε) - F(ua)))
    return r + g + s + n, (r, g, s, n)
end
γη(η) = (1 - η - η^2) / (1 + η)
χ2(η) = (η / sqrt(2) + η^2 / (2(1 - η)))^2 / (1 - η - η^2)^2
function best(f)
    b = (typemax(Int), (0, 0, 0, 0), 0.0)
    for β in range(0.005, 0.995, length=991)
        K, ph = f(β); K < b[1] && (b = (K, ph, β))
    end
    b
end
Kexact(r, ε) = best(β -> calK(r, ε, 1 - β, H4(r, β)))
Ktwo(r, ε, a) = best(β -> calK(r, ε, 1 - β, H4(r, β) + r * log(1 / a); a=a))
Ksk(r, ε, η) = best(β -> calK(r, ε, (1 - β) / 2, (H4(r, β) + r * χ2(η)) / 2; a=γη(η)))
println("== validation (a=1): 71 722 7243 72539 726105")
for r in (10, 100, 1000, 10^4, 10^5)
    K, ph, β = Kexact(r, 1.0); @printf("r=%d K=%d %s beta=%.3f\n", r, K, string(ph), β)
end
println("== two-sided K_a")
for r in (10, 100, 1000), ε in (1.0, 0.1)
    @printf("r=%5d eps=%.1f :", r, ε)
    for a in (1.0, 0.8, 0.5, 0.25)
        K, ph, _ = Ktwo(r, ε, a); @printf("  a=%.2f K=%7d", a, K)
    end
    println()
end
println("== gamma, chi2")
for η in (0.05, 0.1, 0.25); @printf("eta=%.2f gamma=%.4f chi2=%.5f log(1/gamma)=%.4f\n", η, γη(η), χ2(η), log(1 / γη(η))); end
# q so that trA*Pfail <= (ε/2) τ with Pfail <= 2(K+1) e^{K chi2/2} 2^{-q/2}; target = ε ηA/2
qneed(K, η, target) = max(2, ceil(Int, 2 * (log2(2 * (K + 1) / target) + K * χ2(η) / (2log(2)))))
sstack(d, q, η) = (s = ceil(Int, (17 / 2) * (2q + 1) / η); b = ceil(Int, 81 * (d + 2q + 1) / η^2 / s); (s, b, s * b))
srht(d, q, η) = ceil(Int, 216 * (d + 12q^2) / η^2)
println("== sketched: K = Ksk(r, eps/2, eta), target eps*etaA/2, etaA=1e-6, d=K+1 (real) ")
for r in (10, 100, 1000), ε in (1.0, 0.1), η in (0.1, 0.25)
    K, ph, β = Ksk(r, ε / 2, η); K0, _, _ = Kexact(r, ε)
    q = qneed(K, η, ε * 1e-6 / 2)
    s, b, m = sstack(K + 1, q, η); M = srht(K + 1, q, η)
    @printf("r=%5d eps=%.1f eta=%.2f  Kexact(eps)=%7d Ksk=%7d q=%5d  SS s=%6d m=%10d (m/K=%.0f)  SRHT M=%12d (M/K=%.0f)\n", r, ε, η, K0, K, q, s, m, m / K, M, M / K)
end
println("== SRHT best eta (r=100, eps=0.1, etaA=1e-6)")
for r in (10, 100, 1000)
    bestM = (typemax(Int), 0.0, 0, 0)
    for η in range(0.01, 0.25, length=97)
        K, _, _ = Ksk(r, 0.05, η); q = qneed(K, η, 0.1 * 1e-6 / 2); M = srht(K + 1, q, η)
        M < bestM[1] && (bestM = (M, η, K, q))
    end
    @printf("r=%d best SRHT M=%d eta=%.3f K=%d q=%d M/K^1.5=%.1f\n", r, bestM..., bestM[1] / bestM[3]^1.5)
end
println("== chi-square-only scale K/eps^2 vs ours")
for r in (10, 100), ε in (1.0, 0.1)
    K, _, _ = Kexact(r, ε); @printf("r=%d eps=%.1f K=%d K/eps^2=%.0f\n", r, ε, K, K / ε^2)
end
