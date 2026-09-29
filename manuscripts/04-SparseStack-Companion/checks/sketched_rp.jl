# Numerical diagnostic (not a proof): exact RP vs sketched RP (SparseStack, Gaussian) on a
# low-rank-plus-small-tail Gram matrix A = F'F. Reports E tr R_A(J_k)/tau and regression excess.
using LinearAlgebra, Random, Statistics
Random.seed!(20260929)
p, n, r = 300, 400, 5
lam = vcat(ones(r), 2.5e-6 .* (1 .+ rand(p - r)))
Uf = Matrix(qr(randn(p, p)).Q); Vf = Matrix(qr(randn(n, n)).Q)[:, 1:p]
F = Uf * Diagonal(sqrt.(lam)) * Vf'
A = Symmetric(F' * F); τ = sum(sort(eigvals(A), rev = true)[r+1:end]); ηA = τ / tr(A)
function rp(Bm, k)
    R = Matrix(Bm); J = Int[]
    for _ in 1:k
        d = max.(diag(R), 0.0); t = sum(d); t <= 0 && break
        i = findfirst(cumsum(d) .>= rand() * t); i === nothing && (i = argmax(d))
        push!(J, i); R .-= R[:, i] * R[i:i, :] / R[i, i]
    end
    J
end
resid(J) = (Q = Matrix(qr(F[:, J]).Q)[:, 1:length(J)]; norm(F - Q * (Q' * F))^2)
cs(n, b) = (H = zeros(b, n); for i in 1:n; H[rand(1:b), i] = rand((-1.0, 1.0)); end; H)
ss(m, s) = vcat([cs(p, m ÷ s) for _ in 1:s]...) / sqrt(s)
println("eta_A = ", ηA)
for k in (10, 20)
    ex = mean(resid(rp(A, k)) for _ in 1:150) / τ
    for (lab, mk) in (("m=2k", 2), ("m=4k", 4))
        m = mk * k
        vals = Float64[]; regs = Float64[]; valsG = Float64[]
        for _ in 1:150
            P = ss(m, 4); PF = P * F; J = rp(Symmetric(PF' * PF), k)
            push!(vals, resid(J))
            X = pinv(PF[:, J]) * PF; push!(regs, norm(F - F[:, J] * X)^2 / resid(J))
            G = randn(m, p) / sqrt(m); GF = G * F; push!(valsG, resid(rp(Symmetric(GF' * GF), k)))
        end
        println("k=$k $lab: exact RP ", round(ex, digits = 3), "  SparseStack RP ", round(mean(vals) / τ, digits = 3),
                "  Gaussian RP ", round(mean(valsG) / τ, digits = 3), "  same-sketch regression excess ", round(mean(regs), digits = 3),
                "  1+k/(m-k-1) = ", round(1 + k / (m - k - 1), digits = 3))
    end
end
