include("chain_lib.jl")
SEED = 424242
N = 8_000_000
for λ in ([5.0, 3, 2, 1], [10.0, 1, 0.1, 0.01, 0.001], [4.0, 4, 1, 0.5, 0.25, 0.1], [100.0, 10, 1, 1, 1])
    mr, sr = chain_means(λ, 0.5, N); mc, sc = chain_means(λ, 1.0, N)
    println("spectrum ", λ)
    for k in 1:length(λ)-1
        V = Vk(λ, k)
        if k == 2
            er = onestep(λ, 2, 0.5)[1]; ec = onestep(λ, 2, 1.0)[1]
            @printf("  k=%d & %.6g & %.6g (exact) & %.6g (exact) & %.4f & %.4f   [MC %.6g±%.2g, %.6g±%.2g]\n", k, V, er, ec, er / V, ec / V, mr[k], sr[k], mc[k], sc[k])
        else
            @printf("  k=%d & %.6g & %.6g ± %.2g & %.6g ± %.2g & %.4f & %.4f\n", k, V, mr[k], sr[k], mc[k], sc[k], mr[k] / V, mc[k] / V)
        end
    end
end
