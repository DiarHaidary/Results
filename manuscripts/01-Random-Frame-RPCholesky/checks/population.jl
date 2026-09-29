# Gaussian-row population chain: prior V_0 = C, tilt (1+z*Vz)/(1+tr V), V' = V - Vzz*V/(1+z*Vz).
include("common.jl")
using Base.Threads
function popstep(V, cplx, rng)
    r = size(V, 1); F = eigen(Hermitian(V)); v = max.(F.values, 0.0); U = F.vectors
    T = sum(v)
    z = cplx ? randn(rng, ComplexF64, r) : complex.(randn(rng, r))
    if rand(rng) < T / (1 + T)            # tilt component z*Vz: boost coordinate i ∝ v_i (in eigenbasis)
        x = rand(rng) * T; i = 1; acc = v[1]
        while acc < x && i < r; i += 1; acc += v[i]; end
        m2 = cplx ? (randexp(rng) + randexp(rng)) : (randn(rng)^2 + randn(rng)^2 + randn(rng)^2)
        ph = cplx ? cis(2π * rand(rng)) : (rand(rng) < 0.5 ? -1.0 : 1.0)
        z[i] = sqrt(m2) * ph
    end
    y = U * z
    Vy = V * y
    Vn = V - Vy * Vy' / (1 + real(dot(y, Vy)))
    return Hermitian((Vn + Vn') / 2)
end
function popmeans(C, cplx, K, N)
    nt = nthreads(); S = zeros(K, nt); S2 = zeros(K, nt)
    @threads for t in 1:nt
        rng = MersenneTwister(7 * t + (cplx ? 1 : 0))
        for _ in 1:cld(N, nt)
            V = Hermitian(complex.(C))
            for k in 1:K
                V = popstep(Matrix(V), cplx, rng); x = real(tr(V)); S[k, t] += x; S2[k, t] += x^2
            end
        end
    end
    M = cld(N, nt) * nt; m = vec(sum(S, dims = 2)) ./ M
    m, sqrt.((vec(sum(S2, dims = 2)) ./ M .- m .^ 2) ./ M)
end
Ek(c, k) = (e = esym(c); sum(Float64(e[j+1]) / factorial(big(k - j)) for j in 0:min(length(c), k)))
bound(c, k) = (k + 1) * Ek(c, k + 1) / Ek(c, k) - 1
K = 6
for (nm, C) in (("I_2", Matrix(1.0I, 2, 2)), ("diag(1e4,1)", Matrix(Diagonal([1e4, 1.0]))), ("1e4 I_2", Matrix(1e4I, 2, 2)),
               ("diag(10,1,0.1)", Matrix(Diagonal([10, 1, 0.1]))), ("1e4 I_3", Matrix(1e4I, 3, 3)))
    r = size(C, 1); c = diag(C)
    for cplx in (false, true)
        m, s = popmeans(C, cplx, K, 1_000_000)
        println(nm, cplx ? " complex" : " real   ")
        for k in 1:K
            @printf("   k=%d  E tr V_k = %.4f ± %.4f   bound %.4f   %s\n", k, m[k], s[k], bound(Float64.(c), k), k >= r ? @sprintf("r/(k+1-r)=%.4f", r / (k + 1 - r)) : "")
        end
    end
end
# Wishart(r+2) identity: E det(I + x W^{-1}) = L_r(-x) = Σ C(r,l) x^l / l!, W ~ Wish_r(r+2, I) real
rng = MersenneTwister(3)
for r in 2:5
    x = Float64(r); N = 400_000; acc = 0.0; acct = 0.0
    for _ in 1:N
        G = randn(rng, r + 2, r); W = G' * G
        acc += det(I + x * inv(W)); acct += tr(inv(W))
    end
    L = sum(binomial(r, l) * x^l / factorial(l) for l in 0:r)
    @printf("r=%d  MC E det(I+rW^-1) = %.3f  L_r(-r) = %.3f   E tr W^-1 = %.3f (exact r = %d)\n", r, acc / N, L, acct / N, r)
end
