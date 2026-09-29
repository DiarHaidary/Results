include("common.jl")
function lemmaG_gap(c, a, u, v)
    π_ = a .* c ./ sum(a .* c); ut = u .- sum(π_ .* u); vt = v .- sum(π_ .* v)
    Φ(t) = exp(-sum(a .* log1p.(t .* c)))
    lapint(t -> (ct = c ./ (1 .+ t .* c); Φ(t) * (sum(a .* ct .* ut) * sum(a .* ct .* vt) + sum(a .* ct .^ 2 .* ut .* vt))))
end
rng = MersenneTwister(5); nv = 0; worst = -Inf; wc = nothing
for trial in 1:3000
    d = rand(rng, 3:6)
    c = exp.(4 .* randn(rng, d)); a = fill(rand(rng, (0.5, 1.0, 0.1)), d); u = sort(randn(rng, d)); v = sort(randn(rng, d), rev = true)
    g = lemmaG_gap(c, a, u, v) / (norm(u .- mean(u)) * norm(v .- mean(v)) * sum(a .* c))
    if g > worst; global worst = g; global wc = (c, a, u, v); end
    g > 1e-9 && (global nv += 1)
end
println("non-monotone c: violations $nv / 3000, worst normalized ", worst); println(wc)
