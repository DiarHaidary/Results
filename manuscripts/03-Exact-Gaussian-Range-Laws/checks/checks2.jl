# Additional diagnostics: rank ceiling sharpness (hook construction), Grassmannian image dimension,
# single-core TT family, rectangular one-query powers, compressed matrices, frozen designs.
using LinearAlgebra, Random, Statistics, Printf
rng = MersenneTwister(290926)

println("== (H) hook construction: dependent-factor product probe with V'omega ~ N(0, I_{m+1}) ==")
let ns = [3, 2, 4]
    q = length(ns); m = sum(ns .- 1); N = prod(ns)
    # hook coordinates: multi-indices with at most one index different from the last one
    idx(ms) = LinearIndices(Tuple(reverse(ns)))[reverse(ms)...]   # kron ordering: first factor slowest
    hook = Int[]
    push!(hook, idx(ns))
    for j in 1:q, i in 1:ns[j]-1
        ms = copy(ns); ms[j] = i; push!(hook, idx(ms))
    end
    @assert length(unique(hook)) == m + 1
    Vt(om) = om[hook]
    Nsamp = 100000; X = zeros(m + 1, Nsamp)
    maxdev = 0.0
    for t in 1:Nsamp
        z0 = randn(rng); zs = [randn(rng, n - 1) for n in ns]
        xs = Vector{Vector{Float64}}(undef, q)
        xs[1] = [zs[1]; z0]
        for j in 2:q; xs[j] = [zs[j] ./ z0; 1.0]; end
        om = foldl(kron, xs)
        target = zeros(m + 1); target[1] = z0; c = 2
        for j in 1:q, i in 1:ns[j]-1; target[c] = zs[j][i]; c += 1; end
        maxdev = max(maxdev, norm(Vt(om) - target) / norm(target))
        X[:, t] = Vt(om)
    end
    @printf("dims %s: m+1 = %d; max |V'omega - z|/|z| = %.2e ; ||cov - I|| = %.3f\n", string(ns), m + 1, maxdev, opnorm(X*X'/Nsamp - I))
end

println("\n== (I) Grassmannian image dimension of product-probe ranges vs k(R-k) ==")
function grrank(ns, k, R)
    N = prod(ns); V = Matrix(qr(randn(rng, N, R)).Q)[:, 1:R]   # generic R-dim subspace
    xs = [[randn(rng, n) for n in ns] for _ in 1:k]
    Y = V' * hcat([foldl(kron, x) for x in xs]...)
    P = I - Y*pinv(Y)
    cols = Vector{Vector{Float64}}()
    for c in 1:k, j in eachindex(ns), i in 1:ns[j]
        e = zeros(ns[j]); e[i] = 1; ys = copy(xs[c]); ys[j] = e
        dY = zeros(R, k); dY[:, c] = V' * foldl(kron, ys)
        push!(cols, vec(P*dY))
    end
    rank(hcat(cols...); rtol = 1e-9)
end
for (ns, k, R) in (([2,2,2], 1, 4), ([2,2,2], 1, 5), ([2,2,2], 2, 5), ([2,2,2], 2, 6), ([2,2,2], 3, 6), ([2,2,2], 3, 7), ([3,3], 2, 6), ([3,3], 2, 7))
    m = sum(ns .- 1)
    @printf("dims %-8s k=%d R=%d: image dim %2d, k*m = %2d, dim Gr = k(R-k) = %2d, R <= k+m ? %s\n",
            string(ns), k, R, grrank(ns, k, R), k*m, k*(R-k), R <= k + m)
end

println("\n== (J) single-core TT family ==")
let q=3, d=[3, 4, 3], chi=2, t=2, Nsamp=60000
    # omega(i1,i2,i3) = chi^{-(q-1)/2} sum H1(1,i1,a) H2(a,i2,b) H3(b,i3,1)
    function ttprobe()
        H1 = randn(rng, d[1], chi); H2 = randn(rng, chi, d[2], chi); H3 = randn(rng, chi, d[3])
        om = zeros(d[1], d[2], d[3])
        for a in 1:chi, b in 1:chi, i1 in 1:d[1], i2 in 1:d[2], i3 in 1:d[3]
            om[i1, i2, i3] += H1[i1, a]*H2[a, i2, b]*H3[b, i3]
        end
        om ./ chi
    end
    # kron ordering: vec with first index slowest -> use permutedims
    vecc(T) = vec(permutedims(T, (3, 2, 1)))
    # support: u_a = X ⊗ phi_a at site 2, X = e_1 ⊗ e_1 on sites (1,3), phi_a = e_a in R^4, a = 1..3 (s=1)
    R = 3
    U = zeros(prod(d), R)
    for a in 1:R
        T = zeros(d...); T[1, a, 1] = 1.0; U[:, a] = vecc(T)
    end
    X = zeros(R, Nsamp); S2 = zeros(prod(d), prod(d))
    for s in 1:Nsamp
        om = vecc(ttprobe()); X[:, s] = U'*om
        if s <= 20000; S2 .+= om*om'; end
    end
    @printf("||E omega omega' - I|| (20000 samples) = %.3f ; mean ||omega||^2/N = %.3f\n", opnorm(S2/20000 - I), tr(S2/20000)/prod(d))
    Dn = X ./ sqrt.(sum(X.^2, dims=1))
    @printf("direction second moment ||E uu' - I/R|| = %.4f (0 for uniform)\n", opnorm(Dn*Dn'/Nsamp - I/R))
    @printf("fourth moment of first projected coordinate / (Var)^2 = %.2f (Gaussian 3; mixture larger)\n", mean(X[1, :].^4)/mean(X[1, :].^2)^2)
end

println("\n== (K) rectangular one-query power iteration and compressed matrices ==")
function cs_support(n1, n2, R, lam)
    s = length(lam); E = Matrix(qr(randn(rng, n1, n1)).Q); F = Matrix(qr(randn(rng, n2, n2)).Q)
    U = zeros(n1*n2, R)
    for a in 1:R, b in 1:s; U[:, a] .+= sqrt(lam[b]) .* kron(E[:, (a-1)*s+b], F[:, b]); end
    U, sum(lam[b]*F[:, b]*F[:, b]' for b in 1:s)
end
let n1=15, n2=3, R=5
    UL, _ = cs_support(n1, n2, R, [0.6, 0.4]); UR, _ = cs_support(n1, n2, R, [0.5, 0.3, 0.2])
    sig = [5.0, 3.0, 2.0, 1.0, 0.5]; A = UL*Diagonal(sig)*UR'
    hL = randn(rng, n2); hR = randn(rng, n2)
    LL = kron(Matrix(1.0I, n1, n1), hL); LR = kron(Matrix(1.0I, n1, n1), hR)
    y0 = A*kron(randn(rng, n1), randn(rng, n2))
    wL = norm(LL'*y0)^2/norm(y0)^2
    x0 = A'*kron(randn(rng, n1), randn(rng, n2)); wR = norm(LR'*x0)^2/norm(x0)^2
    y = copy(y0); z = copy(y0)
    for t in 1:4
        x = A'*(LL*(LL'*y)/wL)        # one query to A^T
        y = A*(LR*(LR'*x)/wR)         # one query to A
        z = A*(A'*z)
    end
    @printf("4 cycles of A A^T, two product queries per cycle: rel err %.2e\n", norm(y - z)/norm(z))
    Bc = LL'*A*LR/sqrt(wL*wR)
    @printf("singular values of compressed %dx%d matrix: %s\n", n1, n1, string(round.(svdvals(Bc)[1:6], digits=10)))
    # PSD logdet transfer
    U, S0 = cs_support(n1, n2, R, [0.7, 0.3]); Ap = U*Diagonal([4.0, 2.0, 1.0, 0.5, 0.1])*U'
    h = randn(rng, n2); L = kron(Matrix(1.0I, n1, n1), h); w = h'*S0*h; B = L'*Ap*L/w; rho = 0.3; Nn = n1*n2
    @printf("logdet(A+rho I_N) - [(N-n1)log rho + logdet(B+rho I_n1)] = %.2e\n",
            logdet(Ap + rho*I) - ((Nn - n1)*log(rho) + logdet(B + rho*I)))
    # eigenvector lift: A-eigenvector = A L x / (lambda sqrt(w))
    F = eigen(Symmetric(B)); x = F.vectors[:, end]; lam = F.values[end]
    u = Ap*L*x/(lam*sqrt(w))
    @printf("lifted eigenvector: ||A u - lambda u|| = %.2e, ||u|| = %.6f\n", norm(Ap*u - lam*u), norm(u))
end

println("\n== (L) fully frozen GN (both environments shared) vs Gaussian GN ==")
let n1=16, n2=3, R=8, k=4, p=10, Nt=6000
    UL, _ = cs_support(n1, n2, R, [0.8, 0.2]); UR, _ = cs_support(n1, n2, R, [0.5, 0.5])
    sig = [4.0, 2.0, 1.0, 0.7, 0.5, 0.3, 0.2, 0.1]; A = UL*Diagonal(sig)*UR'; Nn = n1*n2
    gn(Y, Psi) = Y*pinv(Psi'*Y)*(Psi'*A)
    ef = Float64[]; eg = Float64[]
    for _ in 1:Nt
        hR = randn(rng, n2); hL = randn(rng, n2)
        Om = kron(randn(rng, n1, k), hR); Psi = kron(randn(rng, n1, p), hL)
        push!(ef, norm(A - gn(A*Om, Psi))^2)
        push!(eg, norm(A - gn(A*randn(rng, Nn, k), randn(rng, Nn, p)))^2)
    end
    @printf("fully frozen: mean %.4f median %.4f | Gaussian: mean %.4f median %.4f\n", mean(ef), median(ef), mean(eg), median(eg))
end

println("\n== (M) exact recovery from m+1 Gaussian product queries (any PSD input of rank m+1) ==")
let ns = [2, 2, 2, 2]
    m = sum(ns .- 1); N = prod(ns); R = m + 1
    U = Matrix(qr(randn(rng, N, R)).Q)[:, 1:R]; A = U*Diagonal(rand(rng, R) .+ 0.1)*U'
    Z = hcat([foldl(kron, [randn(rng, n) for n in ns]) for _ in 1:R]...)
    Y = A*Z; C = Y*pinv(Z'*Y)*Y'
    @printf("N = %d, rank %d, %d product queries: rel err %.2e\n", N, R, R, norm(A - C)/norm(A))
end

println("\n== (N) 6r+8 budget error factor ==")
for r in (1, 2, 5, 10, 100)
    k = 2r + 2; p = 2k + 2
    @printf("r=%3d: k+p = %d (6r+8 = %d), factor = %.4f\n", r, k + p, 6r + 8, (1 + k/(p - k - 1))*(1 + r/(k - r - 1)))
end
