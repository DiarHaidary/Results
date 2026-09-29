# Numerical diagnostics for "Exact Gaussian range laws for structured sketches".
# Diagnostics only; the proofs are in the manuscript.
using LinearAlgebra, Random, Statistics, Printf
rng = MersenneTwister(20260929)
γE = Base.MathConstants.eulergamma

println("== constants ==")
kG = exp(2 - γE) / 2
@printf("kappa_G = %.10f ; 2 sqrt(6 kappa_G) = %.7f\n", kG, 2sqrt(6kG))
@printf("112/pi - 4 = %.4f ; C_new = 112/(3pi) - 4/3 = %.4f ; 3 C_new = %.4f\n", 112/pi - 4, 112/(3pi) - 4/3, 3*(112/(3pi) - 4/3))
@printf("pi^2/4 = %.4f ; pi^2/2 = %.4f\n", pi^2/4, pi^2/2)

# common-Schmidt support u_a = sum_b sqrt(lam_b) e_{a,b} ⊗ f_b  (kron(e,f): f index fast)
function cs_support(n1, n2, R, lam)
    s = length(lam); @assert R*s <= n1
    E = Matrix(qr(randn(rng, n1, n1)).Q); F = Matrix(qr(randn(rng, n2, n2)).Q)
    U = zeros(n1*n2, R)
    for a in 1:R, b in 1:s
        U[:, a] .+= sqrt(lam[b]) .* kron(E[:, (a-1)*s+b], F[:, b])
    end
    S0 = sum(lam[b]*F[:, b]*F[:, b]' for b in 1:s)
    U, S0
end
nys(A, Z) = (Y = A*Z; Y*pinv(Z'*Y)*Y')
function trunc_r(C, r)
    ev, Q = eigen(Symmetric(C)); idx = sortperm(ev, rev=true)[1:r]
    Q[:, idx]*Diagonal(max.(ev[idx], 0))*Q[:, idx]'
end

println("\n== (A) common-Schmidt KR vs Gaussian Nystrom ==")
let n1=12, n2=3, R=6, r=2, k=4, N=20000
    lam = [0.7, 0.3]; U, S0 = cs_support(n1, n2, R, lam)
    μ = [10.0, 5.0, 1.0, 0.5, 0.3, 0.2]; A = U*Diagonal(μ)*U'; τ = sum(μ[r+1:end])
    ekr = Float64[]; eg = Float64[]; ekr_r = Float64[]; eg_r = Float64[]
    for _ in 1:N
        Z = hcat([kron(randn(rng, n1), randn(rng, n2)) for _ in 1:k]...)
        G = randn(rng, n1*n2, k)
        Ck = nys(A, Z); Cg = nys(A, G)
        push!(ekr, tr(A - Ck)); push!(eg, tr(A - Cg))
        push!(ekr_r, tr(A - trunc_r(Ck, r))); push!(eg_r, tr(A - trunc_r(Cg, r)))
    end
    @printf("full: KR mean %.4f  Gauss mean %.4f ; rank-r: KR %.4f Gauss %.4f ; bound (1+r/(k-r-1))tau = %.4f\n",
        mean(ekr), mean(eg), mean(ekr_r), mean(eg_r), (1 + r/(k-r-1))*τ)
    for q in (0.1, 0.5, 0.9)
        @printf("  q=%.1f: KR %.4f Gauss %.4f\n", q, quantile(ekr_r, q), quantile(eg_r, q))
    end
end

println("\n== (B) Hurwitz pair ==")
let N=200000
    M1 = [1.0 0; 0 1]/sqrt(2); M2 = [0.0 -1; 1 0]/sqrt(2)
    sym(M) = (M + M')/2
    @printf("sym(M1'M1)-I/2: %.2e, sym(M2'M2)-I/2: %.2e, sym(M1'M2): %.2e\n",
        norm(sym(M1'M1) - I/2), norm(sym(M2'M2) - I/2), norm(sym(M1'M2)))
    bins = zeros(Int, 8)
    for _ in 1:N
        g = randn(rng, 2); h = randn(rng, 2); x = [g'*M1*h, g'*M2*h]
        θ = mod(atan(x[2], x[1]), pi); bins[min(8, 1 + floor(Int, θ/(pi/8)))] += 1
    end
    println("angle histogram (should be 0.125 each): ", round.(bins ./ N, digits=4))
end

println("\n== (C) head-only counterexample: Var log|slope| ==")
let N=400000
    ls = [log(abs((randn(rng)*randn(rng))/(randn(rng)*randn(rng)))) for _ in 1:N]
    lg = [log(abs(randn(rng)/randn(rng))) for _ in 1:N]
    @printf("structured %.4f (pi^2/2=%.4f), Gaussian %.4f (pi^2/4=%.4f)\n", var(ls), pi^2/2, var(lg), pi^2/4)
end

println("\n== (D) dimension counts: Segre cone and Grassmannian image ==")
# Jacobian rank of (x_1..x_q) -> V' (x_1⊗...⊗x_q)  and of k columns -> chart of Gr(k,R)
function kronall(xs); y = xs[1]; for x in xs[2:end]; y = kron(y, x); end; y; end
function jacrank(f, p; h=1e-6)
    f0 = f(p); J = zeros(length(f0), length(p))
    for i in eachindex(p); pp = copy(p); pp[i] += h; pm = copy(p); pm[i] -= h; J[:, i] = (f(pp) - f(pm))/(2h); end
    s = svdvals(J); count(>(1e-6*s[1]), s)
end
for ns in ([2,2], [2,2,2], [2,2,2,2,2], [3,4], [3,3,3])
    Ntot = prod(ns); m = sum(ns .- 1)
    R = Ntot  # full space, V = I
    split(p) = (idx = cumsum([0; ns]); [p[idx[j]+1:idx[j+1]] for j in eachindex(ns)])
    r1 = jacrank(p -> kronall(split(p)), randn(rng, sum(ns)))
    @printf("n=%s: Segre cone dim (Jacobian rank) %d vs m+1=%d\n", string(ns), r1, m+1)
end
let ns=[2,2,2], k=4
    Ntot = 8; m = 3; R = 8
    function chart(p)
        cols = [kronall([p[(i-1)*6+1:(i-1)*6+2], p[(i-1)*6+3:(i-1)*6+4], p[(i-1)*6+5:(i-1)*6+6]]) for i in 1:k]
        X = hcat(cols...); X2 = X[k+1:end, :]*inv(X[1:k, :]); vec(X2)
    end
    rk = jacrank(chart, randn(rng, 6k))
    @printf("n=[2,2,2], k=%d: Grassmannian image dim %d ; km=%d ; k(R-k)=%d\n", k, rk, k*m, k*(R-k))
end

println("\n== (E) one-query identity, compressed matrix, n1 recovery ==")
let n1=12, n2=4, R=4
    lam = [0.6, 0.3, 0.1]; U, S0 = cs_support(n1, n2, R, lam)
    A = U*Diagonal([10., 3., 1., 0.2])*U'
    h = randn(rng, n2); L = kron(Matrix(1.0I, n1, n1), h); w = h'*S0*h
    @printf("||U'LL'U - wI|| = %.2e\n", norm(U'*L*L'*U - w*I))
    v = A*randn(rng, n1*n2); Mv = reshape(v, n2, n1)'
    q = kron(Mv*h/w, h)
    @printf("one-query rel err %.2e ; w from response rel err %.2e\n", norm(A*q - A*v)/norm(A*v), abs(norm(L'*v)^2/norm(v)^2 - w)/w)
    u = randn(rng, n1*n2); qu = kron(reshape(u, n2, n1)'*h/w, h)
    @printf("unsupported vector: rel err %.3f (identity not expected)\n", norm(A*qu - A*u)/norm(A*u))
    B = L'*A*L/w
    println("eig(B) top 5: ", round.(sort(eigvals(Symmetric(B)), rev=true)[1:5], digits=10))
    Z = L; C = nys(A, Z)
    @printf("Nystrom from n1=%d queries e_i⊗h: rel err %.2e ; |tr B - tr A| = %.2e\n", n1, norm(A - C)/norm(A), abs(tr(B) - tr(A)))
    # powers with frozen environment
    y = A*kron(randn(rng, n1), randn(rng, n2)); z = copy(y)
    for t in 1:5
        y = A*kron(reshape(y, n2, n1)'*h/w, h); z = A*z
    end
    @printf("5 powers via one query each: rel err %.2e\n", norm(y - z)/norm(z))
end
let n1=6, n2=2, n3=2   # third order, entangled f_b, frozen product environment
    E = Matrix(qr(randn(rng, n1, n1)).Q); F = Matrix(qr(randn(rng, 4, 4)).Q)
    lam = [0.7, 0.3]; U3 = zeros(24, 3)
    for a in 1:3, b in 1:2; U3[:, a] .+= sqrt(lam[b]).*kron(E[:, (a-1)*2+b], F[:, b]); end
    A3 = U3*Diagonal([5., 2., 1.])*U3'
    worst = 0.0
    for _ in 1:1000
        h2, h3 = randn(rng, n2), randn(rng, n3); hh = kron(h2, h3); L3 = kron(Matrix(1.0I, n1, n1), hh)
        v3 = A3*randn(rng, 24); w3 = norm(L3'*v3)^2/norm(v3)^2
        q3 = kron(kron(L3'*v3/w3, h2), h3)
        worst = max(worst, norm(A3*q3 - A3*v3)/norm(A3*v3))
    end
    @printf("third order, 1000 frozen product environments: worst rel err %.2e\n", worst)
end

println("\n== (F) generalized Nystrom: Gaussian / calibrated / frozen / raw ==")
let n1=16, n2=3, R=10, k=4, p=10, N=6000
    lam = [0.6, 0.4]; U, S0 = cs_support(n1, n2, 8, lam)   # rank-8 PSD-type symmetric test with same left/right support
    σ = [1.0, 0.8, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05]; A = U*Diagonal(σ)*U'
    function gn(A, Ω, Ψ); Y = A*Ω; Y*pinv(Ψ'*Y)*(Ψ'*A); end
    res = Dict(:gauss=>Float64[], :cal=>Float64[], :frozen=>Float64[], :raw=>Float64[])
    pred = Float64[]
    for _ in 1:N
        Ω = hcat([kron(randn(rng, n1), randn(rng, n2)) for _ in 1:k]...)
        Y = A*Ω; Q = Matrix(qr(Y).Q)
        push!(pred, (1 + k/(p-k-1))*norm(A - Q*(Q'*A))^2)
        Ψg = randn(rng, n1*n2, p); push!(res[:gauss], norm(A - gn(A, Ω, Ψg))^2)
        hs = [randn(rng, n2) for _ in 1:p]; Ψ = hcat([kron(randn(rng, n1), hs[i]) for i in 1:p]...)
        push!(res[:raw], norm(A - gn(A, Ω, Ψ))^2)
        v = Y[:, 1]; Mv = reshape(v, n2, n1)'
        ws = [norm(Mv*hs[i])^2/norm(v)^2 for i in 1:p]
        push!(res[:cal], norm(A - gn(A, Ω, Ψ*Diagonal(1 ./ sqrt.(ws))))^2)
        hL = randn(rng, n2); Ψf = hcat([kron(randn(rng, n1), hL) for _ in 1:p]...)
        push!(res[:frozen], norm(A - gn(A, Ω, Ψf))^2)
    end
    @printf("prediction mean (1+k/(p-k-1))E||(I-QQ')A||^2 = %.4f (se %.4f)\n", mean(pred), std(pred)/sqrt(N))
    for key in (:gauss, :cal, :frozen, :raw)
        @printf("  %-7s mean %.4f (se %.4f) median %.4f\n", string(key), mean(res[key]), std(res[key])/sqrt(N), median(res[key]))
    end
end

println("\n== (G) calibrated trace estimator ==")
let n1=200, n2=3, R=60
    lam = [0.5, 0.3, 0.2]; U, S0 = cs_support(n1, n2, R, lam)
    ev = 1.0 ./ (1:R); A = U*Diagonal(ev)*U'; T = tr(A)
    function Tcal(m, s)
        Z = hcat([kron(randn(rng, n1), randn(rng, n2)) for _ in 1:m]...); Y = A*Z; C = Y*pinv(Z'*Y)*Y'
        v = Y[:, 1]; Mv = reshape(v, n2, n1)'; acc = 0.0
        for j in 1:s
            g, hj = randn(rng, n1), randn(rng, n2); om = kron(g, hj); wj = norm(Mv*hj)^2/norm(v)^2
            acc += (om'*(A*om) - om'*(C*om))/wj
        end
        tr(C) + acc/s
    end
    for (m, s) in ((5, 5), (10, 10), (20, 20))
        ests = [Tcal(m, s) for _ in 1:2000]
        @printf("m=s=%d: mean/T %.4f (se %.4f), Var/T^2 %.3e, bound 2kG/(ms) %.3e\n", m, mean(ests)/T, std(ests)/sqrt(2000)/T, var(ests)/T^2, 2kG/(m*s))
    end
    e = [abs(Tcal(R, 3) - T)/T for _ in 1:20]
    @printf("m = R = %d: max relative error over 20 runs %.2e\n", R, maximum(e))
end

println("\n== (H) left-scaling counterexample ==")
for t in (0.5, 2.0)
    Ψt = [1.0 1; 1 -1]; E0 = Diagonal([1.0, t]); A = Matrix(1.0I, 2, 2); Y = [1.0, 0]
    G1 = Y*pinv(reshape(Ψt*Y, 2, 1))*(Ψt*A); G2 = Y*pinv(reshape(E0*Ψt*Y, 2, 1))*(E0*Ψt*A)
    @printf("t=%.1f: unweighted err^2 %.4f ; weighted err^2 %.4f ; formula %.4f\n", t, norm(A - G1)^2, norm(A - G2)^2, 1 + ((1 - t^2)/(1 + t^2))^2)
end

println("\n== (I) single-core tensor-train family ==")
let q=3, d=4, χ=2, N=100000
    # omega(i1,i2,i3) = χ^{-1} sum_{α,β} H1(i1,α) H2(α,i2,β) H3(β,i3); site t = 2
    # support: u_a = X ⊗ φ_a at site 2 (s=1), X = e_1⊗e_1 on sites 1,3 ; φ_a = e_a, a=1,2
    bins = zeros(Int, 8); m2 = 0.0
    for _ in 1:N
        H1 = randn(rng, d, χ); H2 = randn(rng, χ, d, χ); H3 = randn(rng, χ, d)
        ω = zeros(d, d, d)
        for i1 in 1:d, i2 in 1:d, i3 in 1:d
            ω[i1, i2, i3] = sum(H1[i1, α]*H2[α, i2, β]*H3[β, i3] for α in 1:χ, β in 1:χ)/χ
        end
        x = [ω[1, 1, 1], ω[1, 2, 1]]
        θ = mod(atan(x[2], x[1]), pi); bins[min(8, 1 + floor(Int, θ/(pi/8)))] += 1
        m2 += sum(abs2, ω)
    end
    @printf("E||omega||^2 / d^q = %.4f ; angle histogram %s\n", m2/N/d^q, string(round.(bins ./ N, digits=4)))
end
