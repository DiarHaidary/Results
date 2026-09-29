include("chain_lib.jl")
SEED = 7777
λ = [100.0, 10, 1, 1, 1]
m, s = chain_means(λ, 0.5, 16_000_000)
println("real 16M: ", m, " ± ", s, "   exact k=2 ", onestep(λ,2,0.5)[1], "  V1 ", Vk(λ,1))
