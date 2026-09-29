# Haar-averaged RPCholesky via the spectral Markov chain (Fact 1 + Gamma representation), vs volume sampling.
include("common.jl")

using Base.Threads
gamma_half(rng) = randn(rng)^2 / 2
function step(λ, α, rng)
    d = length(λ); T = sum(λ)
    r = rand(rng) * T; J = 1; acc = λ[1]
    while acc < r && J < d; J += 1; acc += λ[J]; end
    g = α == 0.5 ? [gamma_half(rng) for _ in 1:d] : [randexp(rng) for _ in 1:d]
    g[J] += randexp(rng)                    # size bias: Gamma(α) -> Gamma(α+1)
    u = sqrt.(g)
    w = λ .* u; q = dot(u, w)
    B = Symmetric(Matrix(Diagonal(λ)) - w * w' / q)
    μ = eigvals(B)
    deleteat!(μ, argmin(abs.(μ)))           # remove the kernel direction u
    return max.(μ, 0.0)
end
function chain_means(λ0, α, N)
    n = length(λ0); nt = nthreads()
    S = zeros(n - 1, nt); S2 = zeros(n - 1, nt)
    @threads for t in 1:nt
        rng = MersenneTwister(SEED + 1000 * t + round(Int, 10α))
        for _ in 1:cld(N, nt)
            λ = copy(λ0)
            for k in 1:n-1
                λ = step(λ, α, rng); x = sum(λ)
                S[k, t] += x; S2[k, t] += x^2
            end
        end
    end
    M = cld(N, nt) * nt
    m = vec(sum(S, dims = 2)) ./ M; v = vec(sum(S2, dims = 2)) ./ M .- m .^ 2
    return m, sqrt.(v ./ M)
end
