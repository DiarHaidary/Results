# Head-only counterexample: exact one-column Nystrom error formulas, and projective-uniformity checks
using LinearAlgebra, Random, Statistics, Printf
rng = MersenneTwister(290930)
PhiS(y) = (2y + (2/pi)*(1-y^2)*log(1/y))/(1+y^2)   # structured E err / sqrt(mu1 mu2), y = sqrt(mu2/mu1) <= 1
for (mu1, mu2) in ((1.0, 0.25), (1.0, 0.01), (1.0, 1e-4))
    N = 2_000_000; eg = 0.0; es = 0.0
    for _ in 1:N
        g1, g2, h1, h2, z1, z2 = randn(rng, 6)
        # Nystrom k=1 error for A=diag(mu1,mu2): mu1 mu2 (c1^2+c2^2)/(mu1 c1^2 + mu2 c2^2)
        c1, c2 = g1*h1, g2*h2
        es += mu1*mu2*(c1^2 + c2^2)/(mu1*c1^2 + mu2*c2^2)
        eg += mu1*mu2*(z1^2 + z2^2)/(mu1*z1^2 + mu2*z2^2)
    end
    y = sqrt(mu2/mu1)
    @printf("mu=(%g,%g): Gaussian MC %.4f exact %.4f | structured MC %.4f exact %.4f | ratio %.3f\n",
            mu1, mu2, eg/N, sqrt(mu1*mu2), es/N, sqrt(mu1*mu2)*PhiS(y), PhiS(y))
end
# spherical common-Schmidt probes: direction law of V'omega uniform (second and fourth moments of the direction)
let n1=6, n2=3, R=3
    E = Matrix(qr(randn(rng, n1, n1)).Q); F = Matrix(qr(randn(rng, n2, n2)).Q); lam = [0.6, 0.4]
    U = zeros(n1*n2, R)
    for a in 1:R, b in 1:2; U[:, a] .+= sqrt(lam[b]) .* kron(E[:, (a-1)*2+b], F[:, b]); end
    Ns = 200000; X = zeros(R, Ns)
    for t in 1:Ns
        g = randn(rng, n1); h = randn(rng, n2); g ./= norm(g); h ./= norm(h)
        x = U'*kron(g, h); X[:, t] = x/norm(x)
    end
    m4 = mean(X[1, :].^4)
    @printf("spherical factors: ||E uu' - I/R|| = %.4f ; E u1^4 = %.4f (uniform on S^%d: %.4f)\n",
            opnorm(X*X'/Ns - I/R), m4, R-1, 3/(R*(R+2)))
end
