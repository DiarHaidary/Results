# α-family of spectral chains: α→0 successive sampling (diagonal RP), α=1/2 real Haar, α=1 complex Haar, α→∞ volume.
include("common.jl")
using Base.Threads
function gam(rng, a)
    if a < 1
        return gam(rng, a + 1) * rand(rng)^(1 / a)
    end
    d = a - 1 / 3; c = 1 / sqrt(9d)
    while true
        x = randn(rng); v = (1 + c * x)^3
        v <= 0 && continue
        log(rand(rng)) < 0.5x^2 + d - d * v + d * log(v) && return d * v
    end
end
function stepα(λ, α, rng)
    d = length(λ); T = sum(λ); r = rand(rng) * T; J = 1; acc = λ[1]
    while acc < r && J < d; J += 1; acc += λ[J]; end
    g = [gam(rng, α + (i == J)) for i in 1:d]
    u = sqrt.(g); w = λ .* u; q = dot(u, w)
    μ = eigvals(Symmetric(Matrix(Diagonal(λ)) - w * w' / q)); deleteat!(μ, argmin(abs.(μ)))
    return max.(μ, 0.0)
end
function diagRP(λ0, rng)   # successive sampling
    λ = copy(λ0); out = Float64[]
    for k in 1:length(λ0)-1
        T = sum(λ); r = rand(rng) * T; J = 1; acc = λ[1]
        while acc < r && J < length(λ); J += 1; acc += λ[J]; end
        deleteat!(λ, J); push!(out, sum(λ))
    end
    out
end
function means(λ0, α, N)
    n = length(λ0); nt = nthreads(); S = zeros(n - 1, nt)
    @threads for t in 1:nt
        rng = MersenneTwister(31 * t + round(Int, 100α))
        for _ in 1:cld(N, nt)
            if α == 0
                S[:, t] .+= diagRP(λ0, rng)
            else
                λ = copy(λ0)
                for k in 1:n-1
                    λ = stepα(λ, α, rng); S[k, t] += sum(λ)
                end
            end
        end
    end
    vec(sum(S, dims = 2)) ./ (cld(N, nt) * nt)
end
for λ in ([5.0, 3, 2, 1], [100.0, 10, 1, 1, 1])
    println("spectrum ", λ, "  volume ", [round(Vk(λ, k), digits = 5) for k in 1:length(λ)-1])
    for α in (0.0, 0.1, 0.5, 1.0, 2.0, 8.0, 50.0)
        m = means(λ, α, 2_000_000)
        println("  α=", α, "  ", round.(m, digits = 4), (α > 0 ? "  exact one-step k=2: $(round(onestep(λ,2,α)[1],digits=5))" : ""))
    end
end
