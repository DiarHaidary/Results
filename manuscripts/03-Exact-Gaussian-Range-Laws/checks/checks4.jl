# Head-only ratio Phi(y) > 1 on (0,1); surrogate minimal right inverse; raw GN constants
using LinearAlgebra, Random, Printf
Phi(y) = (2y + (2/pi)*(1-y^2)*log(1/y))/(1+y^2)
ys = range(1e-8, 1-1e-8, length=2_000_001)
@printf("min over (0,1) of Phi(y)-1 = %.3e (at y=%.6f); Phi(0.5)=%.4f Phi(0.1)=%.4f Phi(1e-3)=%.4f\n",
        minimum(Phi.(ys) .- 1), ys[argmin(Phi.(ys))], Phi(0.5), Phi(0.1), Phi(1e-3))
g(y) = (2/pi)*(1+y)*log(1/y) - (1-y)
@printf("min g(y)/(1-y) = %.4f\n", minimum(g.(ys[1:end-1]) ./ (1 .- ys[1:end-1])))
# surrogate: ||D^2 G1' (G1 D^2 G1')^{-1}||_F >= ||pinv(G1)||_F
rng = MersenneTwister(1)
worst = Inf
for _ in 1:10000
    r, k = rand(rng, 1:4), rand(rng, 5:9); G1 = randn(rng, r, k); D = Diagonal(exp.(randn(rng, k)))
    X = D^2*G1'*inv(G1*D^2*G1'); global worst = min(worst, norm(X) - norm(pinv(G1)))
end
@printf("min ||D(G1 D)^+||_F - ||G1^+||_F over 10000 trials = %.3e (>= 0)\n", worst)
Cnew = 112/(3pi) - 4/3
@printf("C_new = %.6f, 3C_new = %.4f, 1 + 3C_new/32 = %.4f\n", Cnew, 3Cnew, 1 + 3Cnew/32)
@printf("complex check: E||G1^+||^2 for complex r x k: r/(k-r)\n")
let r=2, k=6, N=200000
    s = 0.0
    for _ in 1:N
        G = (randn(rng, r, k) + im*randn(rng, r, k))/sqrt(2); s += real(tr(inv(G*G')))
    end
    @printf("MC %.4f vs r/(k-r) = %.4f\n", s/N, r/(k-r))
end
