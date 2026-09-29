# (R) raw independent-column GN: exact coefficient B_{p,k} vs bound mu2*Cnew*k/p
using LinearAlgebra, Random, Printf
rng = MersenneTwister(290931)
Cnew = 112/(3pi) - 4/3
for (name, lam) in (("rank-one S0", [1.0]), ("S0 = I2/2", [0.5, 0.5]), ("S0 = I4/4", fill(0.25, 4)))
    mu2 = 1 + 2*sum(lam.^2)
    for (k, p) in ((4, 20), (4, 40), (8, 36), (8, 72))
        N = 20000; B = 0.0; tinv = 0.0
        for _ in 1:N
            Xi = randn(rng, p, k); w = [sum(lam .* randn(rng, length(lam)).^2) for _ in 1:p]
            H = Xi' * Diagonal(w) * Xi; Hi = inv(H)
            B += tr(Hi * Xi' * Diagonal(w.^2) * Xi * Hi); tinv += tr(Hi)
        end
        @printf("%-12s k=%d p=%3d: B = %.4f  mu2*E tr H^-1 = %.4f  bound mu2*Cnew*k/p = %.4f  ratio bound/B = %.1f\n",
                name, k, p, B/N, mu2*tinv/N, mu2*Cnew*k/p, mu2*Cnew*k/p/(B/N))
    end
end
