include("common.jl")
rng = MersenneTwister(11)
n = 6; λ = [5.0, 3, 2, 1, 0.5, 0.1]; Q = qr(randn(rng, n, n)).Q * I; A = Q * Diagonal(λ) * Q'; A = (A + A') / 2
# (1) det-tilted Gaussian sketch: E[det(Ω'AΩ) tr R_Ω] / E[det(Ω'AΩ)] = V_m
for m in 1:4
    num = 0.0; den = 0.0
    for _ in 1:400_000
        Ω = randn(rng, n, m); G = Ω' * A * Ω; dG = det(G)
        R = A - A * Ω * (G \ (Ω' * A)); num += dG * tr(R); den += dG
    end
    @printf("m=%d  det-tilted E tr R = %.4f   V_m = %.4f   (E det / (m! e_m) = %.4f)\n", m, num / den, Vk(λ, m), den / 400_000 / (factorial(m) * Float64(esym(λ)[m+1])))
end
# (2) dilation: redundant tight frame Z (n×K), weights q; Ã = D^{1/2} Z' A Z D^{1/2}; compare residual traces for a pivot set
K = 10; Z0 = randn(rng, n, K); q = rand(rng, K) .+ 0.5
S0 = Z0 * Diagonal(q) * Z0'; Z = inv(sqrt(Symmetric(S0))) * Z0      # now Σ q_i z_i z_i' = I
@printf("tightness error %.2e\n", norm(Z * Diagonal(q) * Z' - I))
At = Diagonal(sqrt.(q)) * Z' * A * Z * Diagonal(sqrt.(q))
S = [2, 7, 9]; ZS = Z[:, S]
R = A - A * ZS * ((ZS' * A * ZS) \ (ZS' * A))
Rt = At - At[:, S] * (At[S, S] \ At[S, :])
@printf("tr R = %.10f  tr(Ã/Ã_SS) = %.10f ; max |diag(Ã/S) - q_i z_i'Rz_i| = %.2e\n", tr(R), tr(Rt), maximum(abs.(diag(Rt) .- q .* [Z[:, i]' * R * Z[:, i] for i in 1:K])))
println("nonzero spectra: ", round.(sort(filter(x -> x > 1e-10, eigvals(Symmetric(R))), rev = true), digits = 8))
println("                 ", round.(sort(filter(x -> x > 1e-10, eigvals(Symmetric(Rt))), rev = true), digits = 8))
