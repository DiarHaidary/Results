include("common.jl")
println("Anchor (5,3,2,1), k=2:")
for (α, nm) in ((0.5, "real"), (1.0, "complex"))
    L, R, D, T1, T2 = onestep([5.0, 3, 2, 1], 2, α)
    @printf("  %-8s LHS=%.9f  V_2=%.9f  D=%.6f  T1=%.3e T2=%.3e  margin (RHS-LHS)/k/ES*... -D/E[S]=%.6f\n", nm, L, R, D, T1, T2, -D / (α * 11))
end
println("k=3 one-step for (5,3,2,1):")
for α in (0.5, 1.0)
    L, R, _ = onestep([5.0, 3, 2, 1], 3, α); @printf("  α=%.1f  LHS=%.6f  V_3=%.6f\n", α, L, R)
end
# random scan: signs of T1,T2 and monotonicity
rng = MersenneTwister(1)
nviol = 0; ncase = 0; worst = -Inf
for trial in 1:1500
    d = rand(rng, 3:9); k = rand(rng, 2:d-1); α = rand(rng, (0.05, 0.2, 0.5, 1.0, 2.0, 5.0))
    λ = exp.(rand(rng, d) .* rand(rng, (1.0, 5.0, 15.0)))
    rand(rng) < 0.3 && (λ[2] = λ[1])
    L, R, D, T1, T2 = onestep(λ, k, α)
    global ncase += 1
    sc = abs(T1) + abs(T2) + 1e-300
    if T1 > 1e-12 * sc || T2 > 1e-12 * sc || L > R * (1 + 1e-10)
        global nviol += 1
    end
    global worst = max(worst, (L - R) / R)
end
@printf("scan: %d cases, %d sign/S_k violations, max relative (LHS-RHS)/RHS = %.2e\n", ncase, nviol, worst)

# Lemma G with heterogeneous shapes: Monte-Carlo-free check via Laplace formula
function lemmaG_gap(c, a, u, v)   # returns D' = E[(y.ũ)(y.ṽ)/R] with π ∝ a c
    π_ = a .* c ./ sum(a .* c); ut = u .- sum(π_ .* u); vt = v .- sum(π_ .* v)
    Φ(t) = exp(-sum(a .* log1p.(t .* c)))
    lapint(t -> (ct = c ./ (1 .+ t .* c); Φ(t) * (sum(a .* ct .* ut) * sum(a .* ct .* vt) + sum(a .* ct .^ 2 .* ut .* vt))))
end
# brute-force MC check of the heterogeneous Laplace formula
function gam(rng, a)  # Marsaglia-Tsang for a>=1, boost for a<1
    if a < 1
        return gam(rng, a + 1) * rand(rng)^(1 / a)
    end
    d = a - 1 / 3; c = 1 / sqrt(9d)
    while true
        x = randn(rng); v = (1 + c * x)^3
        v <= 0 && continue
        uu = rand(rng)
        log(uu) < 0.5x^2 + d - d * v + d * log(v) && return d * v
    end
end
c = [0.3, 0.7, 1.5, 4.0]; a = [0.4, 1.3, 0.7, 2.0]; u = [0.1, 0.5, 0.6, 2.0]; v = [3.0, 1.0, 0.2, -1.0]
π_ = a .* c ./ sum(a .* c); ut = u .- sum(π_ .* u); vt = v .- sum(π_ .* v)
acc = 0.0; N = 2_000_000
for _ in 1:N
    y = [c[i] * gam(rng, a[i]) for i in 1:4]
    global acc += dot(y, ut) * dot(y, vt) / sum(y)
end
@printf("heterogeneous Lemma G: Laplace D'=%.5f  MC D'=%.5f\n", lemmaG_gap(c, a, u, v), acc / N)
nv = 0; nvr = 0
for trial in 1:2000
    d = rand(rng, 2:8)
    c = sort(exp.(3 .* randn(rng, d))); a = exp.(randn(rng, d)); u = sort(randn(rng, d)); v = sort(randn(rng, d), rev = true)
    g = lemmaG_gap(c, a, u, v); g > 1e-10 * (norm(u) * norm(v) * sum(a .* c)) && (global nv += 1)
    g2 = lemmaG_gap(reverse(c), a, u, v); g2 > 1e-10 * (norm(u) * norm(v) * sum(a .* c)) && (global nvr += 1)
end
println("Lemma G heterogeneous shapes: violations = $nv / 2000;  with c reversed: $nvr / 2000")
