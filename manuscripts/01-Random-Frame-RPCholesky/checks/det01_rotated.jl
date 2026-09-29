# Haar-rotated RPCholesky at budget k = r on the spectrum (1, η, ..., η^r), n = r+1.
# Exact over the pivot sequence (dynamic program over principal minors), Monte Carlo over the frame. BigFloat.
using LinearAlgebra, Random, Statistics, Printf, Base.Threads
PREC = parse(Int, get(ENV, "PREC", "512")); setprecision(BigFloat, PREC)
function haarbig(n, cplx, rng)
    G = cplx ? Complex{BigFloat}.(randn(rng, ComplexF64, n, n)) : BigFloat.(randn(rng, n, n))
    F = qr(G); Q = Matrix(F.Q); d = diag(F.R)
    return Q * Diagonal(d ./ abs.(d))
end
function ErR(A, r)
    n = size(A, 1)
    mins = Vector{BigFloat}(undef, 2^n)
    for mask in 0:2^n-1
        S = [i for i in 1:n if (mask >> (i - 1)) & 1 == 1]
        mins[mask+1] = isempty(S) ? one(BigFloat) : real(det(A[S, S]))
    end
    P = Dict{Int,BigFloat}(0 => one(BigFloat))
    for s in 0:r-1
        Pn = Dict{Int,BigFloat}()
        for (m, p) in P
            outs = [i for i in 1:n if (m >> (i - 1)) & 1 == 0]
            ws = [mins[(m|(1 << (i - 1)))+1] for i in outs]
            T = sum(ws)
            for (j, i) in enumerate(outs)
                m2 = m | (1 << (i - 1)); Pn[m2] = get(Pn, m2, zero(BigFloat)) + p * ws[j] / T
            end
        end
        P = Pn
    end
    return sum(p * mins[end] / mins[m+1] for (m, p) in P)
end
η = parse(BigFloat, get(ENV, "ETA", "1e-4"))
for cplx in (false, true), r in 2:6
    n = r + 1; λ = [η^i for i in 0:r]; τ = λ[end]
    e = zeros(BigFloat, n + 1); e[1] = 1
    for x in λ, k in n:-1:1; e[k+1] += x * e[k]; end
    Vr = (r + 1) * e[r+2] / e[r+1]
    F = parse(Int, get(ENV, "FRAMES", "2000"))
    out = zeros(F)
    @threads for t in 1:F
        rng = MersenneTwister(97 * t + 13 * r + (cplx ? 5 : 0))
        Q = haarbig(n, cplx, rng)
        A = Q * Diagonal(λ) * Q'; A = (A + A') / 2
        out[t] = Float64(real(ErR(A, r)) / τ)
    end
    @printf("%s r=%d  mean E_RP tr R_r/τ = %.4f ± %.4f  median %.3f  max %.2f   V_r/τ = %.6f   2^r = %d\n",
            cplx ? "complex" : "real   ", r, mean(out), std(out) / sqrt(F), median(out), maximum(out), Float64(Vr / τ), 2^r)
end
