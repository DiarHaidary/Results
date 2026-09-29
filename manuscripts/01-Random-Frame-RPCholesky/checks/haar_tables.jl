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
SEED = parse(Int, get(ENV, "SEED", "0"))
N = parse(Int, get(ENV, "NS", "4000000"))
for λ in ([5.0, 3, 2, 1], [10.0, 1, 0.1, 0.01, 0.001], [4.0, 4, 1, 0.5, 0.25, 0.1], [100.0, 10, 1, 1, 1])
    mr, sr = chain_means(λ, 0.5, N); mc, sc = chain_means(λ, 1.0, N)
    println("spectrum ", λ)
    for k in 1:length(λ)-1
        V = Vk(λ, k)
        exact2 = k == 2 ? (onestep(λ, 2, 0.5)[1], onestep(λ, 2, 1.0)[1]) : (NaN, NaN)
        @printf("  k=%d  real %.5f ± %.5f   complex %.5f ± %.5f   volume %.5f   ratios %.4f %.4f   exact(k=2) %.6f %.6f\n",
                k, mr[k], sr[k], mc[k], sc[k], V, mr[k] / V, mc[k] / V, exact2...)
    end
end
