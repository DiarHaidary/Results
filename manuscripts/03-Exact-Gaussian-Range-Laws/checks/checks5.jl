# (O) near-isotropic left covariance: exact GN error formula; (P) rectangular exact recovery from 2 n1 queries
# (Q) Phi(y) closed form vs numerical integral against product-of-Cauchy density
using LinearAlgebra, Random, Statistics, Printf
rng = MersenneTwister(290931)
println("== (O) GN with left covariance Sigma on the column support ==")
let R=6, k=2, p=8, Nc=5, Nt=200000
    M = randn(rng, R, Nc)                  # A = U M with U = I_R (support coordinates)
    Qc = Matrix(qr(randn(rng, R, R)).Q)    # adapted basis: first k columns = Q
    η = 0.3; Δ = randn(rng, R, R); Δ = (Δ + Δ')/2; Δ .*= η/opnorm(Δ); Σ = I + Δ
    Σq = Qc'*Σ*Qc; Mq = Qc'*M
    S11 = Σq[1:k,1:k]; S12 = Σq[1:k,k+1:end]; S22 = Σq[k+1:end,k+1:end]
    D = S11\S12; Sc = S22 - S12'*(S11\S12); Ec = Mq[k+1:end,:]
    formula = norm(Ec)^2 + norm(D*Ec)^2 + tr(inv(S11))/(p-k-1)*tr(Ec'*Sc*Ec)
    bound = (1 + η^2/(1-η)^2 + (1+η)/(1-η)*k/(p-k-1))*norm(Ec)^2
    C = cholesky(Symmetric(Σq)).L; acc = 0.0
    for _ in 1:Nt
        Ψt = randn(rng, p, R)*C'            # rows ~ N(0, Σq)  (Psi^T U in adapted coords)
        X = Ψt[:,1:k]; Ahat_top = X \ (Ψt*Mq)  # Q (Psi^T Q)^+ Psi^T A, coordinates on Q
        err = norm(Mq[1:k,:] - Ahat_top)^2 + norm(Ec)^2
        acc += err
    end
    @printf("MC %.4f  formula %.4f  bound %.4f\n", acc/Nt, formula, bound)
end
println("\n== (P) rectangular exact recovery from n1 right + n1 left product queries ==")
function cs_support(n1, n2, R, lam)
    s = length(lam); E = Matrix(qr(randn(rng, n1, n1)).Q); F = Matrix(qr(randn(rng, n2, n2)).Q)
    U = zeros(n1*n2, R)
    for a in 1:R, b in 1:s; U[:, a] .+= sqrt(lam[b]) .* kron(E[:, (a-1)*s+b], F[:, b]); end
    U
end
let n1=15, n2=3, R=5
    UL = cs_support(n1, n2, R, [0.6, 0.4]); UR = cs_support(n1, n2, R, [0.5, 0.3, 0.2])
    A = UL*Diagonal([5., 3., 2., 1., .5])*UR'
    LR = kron(Matrix(1.0I, n1, n1), randn(rng, n2)); LL = kron(Matrix(1.0I, n1, n1), randn(rng, n2))
    Y = A*LR; Ahat = Y*pinv(LL'*Y)*(LL'*A)
    @printf("rel err %.2e\n", norm(A - Ahat)/norm(A))
end
println("\n== (Q) Phi(y) closed form vs quadrature ==")
Phi(y) = (2y + (2/pi)*(1-y^2)*log(1/y))/(1+y^2)
fprod(t) = (2/pi^2)*log(abs(t))/(t^2-1)     # density of product of two standard Cauchy
for y in (0.5, 0.1, 0.01)
    # integrate y(1+t^2)/(1+y^2 t^2) f(t) over t>0, times 2, substitution t = e^s
    h = 1e-4; s = -40:h:40; acc = 0.0
    for si in s
        t = exp(si); ft = abs(t-1) < 1e-8 ? 1/pi^2 : fprod(t)
        acc += y*(1+t^2)/(1+y^2*t^2)*ft*t*h
    end
    @printf("y=%.2f: quadrature %.6f closed form %.6f\n", y, 2acc, Phi(y))
end
