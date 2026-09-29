# Diagnostic only: RP on the sketched Gram B = F'Pi'Pi F versus exact RP on A = F'F.
using LinearAlgebra, Random, Printf, Statistics, SparseArrays
Random.seed!(11)
function make_F(D, n, r; headscale=1e3, tailrank=200)
    Uh = Matrix(qr(randn(D, r)).Q)[:, 1:r]; Vh = Matrix(qr(randn(n, r)).Q)[:, 1:r]
    H = Uh * Diagonal(headscale .* (1:r) .^ -1.0) * Vh'
    T = randn(D, tailrank) * randn(tailrank, n) ./ sqrt(tailrank * n)   # spread-out tail
    # add some coherent "spiky" tail columns
    F = H + T
    return F
end
function rp_select(G::Matrix{Float64}, k::Int)
    n = size(G,1); R = copy(G); J = Int[]
    for t in 1:k
        d = max.(diag(R), 0.0); T = sum(d); T <= 1e-12*tr(G) && break
        u = rand() * T; i = findfirst(cumsum(d) .>= u); i === nothing && (i = argmax(d))
        push!(J, i); c = R[:, i] / sqrt(R[i, i]); R .-= c * c'
    end
    return J
end
resid_true(F, J) = (isempty(J) ? sum(abs2, F) : (Q = Matrix(qr(F[:, J]).Q)[:, 1:length(J)]; sum(abs2, F - Q * (Q' * F))))
function sparsestack(m, D, s)
    b = m ÷ s; rows = Int[]; cols = Int[]; vals = Float64[]
    for g in 1:s, i in 1:D
        push!(rows, (g - 1) * b + rand(1:b)); push!(cols, i); push!(vals, rand((-1.0, 1.0)) / sqrt(s))
    end
    return Matrix(sparse(rows, cols, vals, s * b, D))
end
D, n, r = 300, 400, 5
F = make_F(D, n, r); A = F' * F; ev = sort(eigvals(Symmetric(A)), rev=true)
tau = sum(ev[r+1:end]); @printf("tr A / tau_r = %.1f (eta_A = %.2e)\n", tr(A) / tau, tau / tr(A))
trials = 150
for k in (10, 20, 40)
    ex = mean(resid_true(F, rp_select(A, k)) for _ in 1:trials) / tau
    @printf("k=%3d exact RP: E tr R_A / tau = %.3f\n", k, ex)
    for mult in (1, 2, 4), kind in (:gauss, :sstack)
        m = mult * k; errs = Float64[]; regr = Float64[]
        for _ in 1:trials
            Pi = kind == :gauss ? randn(m, D) ./ sqrt(m) : sparsestack(max(m, 8), D, 4)
            PF = Pi * F; B = PF' * PF; J = rp_select(B, k)
            e = resid_true(F, J); push!(errs, e)
            XB = (PF[:, J]) \ PF                      # same-sketch regression
            push!(regr, sum(abs2, F - F[:, J] * XB) / e)
        end
        @printf("   m=%3d %-6s sketched RP: E tr R_A / tau = %.3f   same-sketch regression / optimal: median %.3f  max %.3f\n",
                m, String(kind), mean(errs) / tau, median(regr), maximum(regr))
    end
end
