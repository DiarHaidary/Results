# Diagnostic (numerical, not a proof): sketched RP on B = (Pi F)^*(Pi F) versus exact RP on A = F^*F.
# Records along the sketched path the quantities used in the proof:
#   xbar = trR_B/trR_A, the chi-square V(J) = sum_i pi_i (x_i-1)^2 with x_i = (R_B)_ii/(R_A)_ii,
#   the weighted bracket for the trace-decrease weights, and the same-sketch reconstruction excess.
using LinearAlgebra, Random, Statistics, SparseArrays, Printf
Random.seed!(2026)
function make_F(D, n, r)
    Uh = Matrix(qr(randn(D, r)).Q)[:, 1:r]; Vh = Matrix(qr(randn(n, r)).Q)[:, 1:r]
    H = Uh * Diagonal(1e3 .* (1:r) .^ -1.0) * Vh'
    T = randn(D, 200) * randn(200, n) ./ sqrt(200 * n)
    return H + T
end
function sparsestack(m, D, s)
    b = max(1, m ÷ s); rows = Int[]; cols = Int[]; vals = Float64[]
    for g in 1:s, i in 1:D
        push!(rows, (g - 1) * b + rand(1:b)); push!(cols, i); push!(vals, rand((-1.0, 1.0)) / sqrt(s))
    end
    return Matrix(sparse(rows, cols, vals, s * b, D))
end
pick(d) = (T = sum(d); u = rand() * T; i = findfirst(cumsum(d) .>= u); i === nothing ? argmax(d) : i)
schur!(R, i) = (c = R[:, i] / sqrt(R[i, i]); R .-= c * c'; R)
function run_exact(A, K)
    R = copy(A); J = Int[]
    for t in 1:K
        d = max.(diag(R), 0.0); sum(d) <= 1e-13 * tr(A) && break
        i = pick(d); push!(J, i); schur!(R, i)
    end
    return tr(R)
end
function run_sketched(A, B, K)
    RA = copy(A); RB = copy(B); J = Int[]
    Vmax = 0.0; xdev = 0.0; wdev = 0.0; chi = 0.0; rdev = 0.0; chisum = 0.0
    for t in 1:K
        dA = max.(diag(RA), 0.0); dB = max.(diag(RB), 0.0)
        TA = sum(dA); TB = sum(dB); TB <= 1e-13 * tr(B) && break
        mask = dA .> 1e-12 * TA
        x = zeros(length(dA)); x[mask] = dB[mask] ./ dA[mask]
        π = dA ./ TA
        V = sum(π[mask] .* (x[mask] .- 1) .^ 2); Vmax = max(Vmax, V)
        xdev = max(xdev, abs(TB / TA - 1))
        c = zeros(length(dA)); c[mask] = vec(sum(abs2, RA[:, mask], dims=1)) ./ dA[mask]   # decrease weights
        w = dA .* c
        wdev = max(wdev, abs(sum(w .* x) / sum(w) - 1)); xb = TB/TA; c2 = sum(π[mask] .* (x[mask] .- xb) .^ 2)/xb^2; chi = max(chi, c2); chisum += log1p(c2); rdev = max(rdev, abs(sum(w .* x)/sum(w)/xb - 1))
        i = pick(dB); push!(J, i); schur!(RA, i); schur!(RB, i)
    end
    return tr(RA), J, Vmax, xdev, wdev, chi, rdev, chisum
end
D, n, r = 2000, 400, 5
F = make_F(D, n, r); A = Symmetric(F' * F) |> Matrix
ev = sort(eigvals(Symmetric(A)), rev=true); tau = sum(ev[r+1:end])
@printf("n=%d D=%d r=%d  tr A/tau_r = %.1f\n", n, D, r, tr(A) / tau)
trials = 100
for K in (20, 40)
    ex = mean(run_exact(A, K) for _ in 1:trials) / tau
    @printf("K=%d exact RP: E trR/tau = %.4f\n", K, ex)
    for mult in (2, 4, 8), kind in (:gauss, :sstack)
        m = mult * K
        C2 = Float64[]; RD = Float64[]; CS = Float64[]; errs = Float64[]; Vs = Float64[]; xd = Float64[]; wd = Float64[]; exc = Float64[]
        for _ in 1:trials
            Pi = kind == :gauss ? randn(m, D) ./ sqrt(m) : sparsestack(m, D, 8)
            PF = Pi * F; B = PF' * PF
            e, J, V, x1, w1, c2, rd, cs = run_sketched(A, B, K); push!(C2, c2); push!(RD, rd); push!(CS, cs)
            push!(errs, e); push!(Vs, V); push!(xd, x1); push!(wd, w1)
            XB = PF[:, J] \ PF
            push!(exc, sum(abs2, F - F[:, J] * XB) / e)
        end
        @printf("  m=%4d %-6s E trR_A/tau=%.4f | path max chi2 V: median %.3f | max|trB/trA-1| median %.3f | decrease-bracket dev median %.3f | reuse excess median %.3f (K/(m-K-1)+1=%.3f)\n      per-step chi2 Var(x)/xbar^2 path-max median %.4f | sum log(1+chi2) median %.3f | decrease-weight ratio dev path-max median %.4f\n",
            m, String(kind), mean(errs) / tau, median(Vs), median(xd), median(wd), median(exc), 1 + K / (m - K - 1), median(C2), median(CS), median(RD))
    end
end
