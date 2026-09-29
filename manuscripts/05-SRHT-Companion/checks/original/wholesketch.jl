# Exact enumeration of the real two-round and one-round flat Walsh laws (diagnostic).
using LinearAlgebra, Random
function hadamard(n); H = reshape([1],1,1); while size(H,1) < n; H = [H H; H -H]; end; H; end
subsets(n, M) = [findall(==(1), digits(b, base=2, pad=n)) for b in 0:(2^n-1) if count_ones(b) == M]
signs(n) = [[(b >> i) & 1 == 1 ? -1 : 1 for i in 0:n-1] for b in 0:(2^n-1)]
nysres(K, S) = (Y = K*S'; K - Y*pinv(S*K*S')*Y')
function run(n, M, K; oneround=false, T=Float64, R0=nothing)
    H = T.(hadamard(n)) ./ (T == Float64 ? sqrt(n) : T(2))
    In = Matrix{T}(I, n, n); cnt = 0; Etr = zero(T); Enys = 0.0; Ecov = zeros(T,n,n); Ediag = zero(T)
    for x in signs(n), y in (oneround ? [ones(Int,n)] : signs(n)), Ts in subsets(n, M)
        Q = oneround ? H*Diagonal(T.(x)) : H*Diagonal(T.(y))*H*Diagonal(T.(x))
        P = Q[Ts, :]; G = (T(n)/M) .* (P'*P)
        Etr += tr(K*(G-In)*K*(G-In)); Enys += norm(nysres(Float64.(K), sqrt(n/M) .* Float64.(P)))^2
        if R0 !== nothing
            Ecov += (G-In)*R0*(G-In); h = diag(R0*G); Ediag += sum((h - diag(R0)).^2)
        end
        cnt += 1
    end
    Etr /= cnt; Enys /= cnt; Ecov ./= cnt; Ediag /= cnt
    lam = T(n-M)/(M*(n-1)); tK = tr(K); fK = sum(abs2, K); dK = sum(abs2, diag(K))
    one = lam*(tK^2 + fK - 2dK); two = lam*((1-T(4)/n)*fK + (1-T(2)/n)*tK^2 + (T(4)/n)*dK)
    simple = 2*(1-M/n)*Float64(tK)^2/M
    out = (Etr=Float64(Etr), pred=Float64(oneround ? one : two), one=Float64(one), two=Float64(two), nys=Enys, simple=simple)
    if R0 !== nothing
        pc = lam*((1-T(4)/n)*R0 + (1-T(2)/n)*tr(R0)*In + (T(4)/n)*Diagonal(diag(R0)))
        pd = lam*((1-T(2)/n)*sum(abs2,R0) + sum(abs2,diag(R0)))
        out = merge(out, (covdiff=Float64(maximum(abs.(Ecov-pc))), diagdiff=Float64(Ediag-pd)))
    end
    out
end
Random.seed!(7)
n = 4; Bk = rand(-3:3, n, 3); K = Rational{BigInt}.(Bk*Bk'); Br = rand(-2:2,n,n); R0 = Rational{BigInt}.(Br*Br')
for M in 1:3
    o = run(n, M, K; T=Rational{BigInt}, R0=R0); println("n=4 M=$M two-round: ", o)
    o = run(n, M, K; T=Rational{BigInt}, oneround=true); println("n=4 M=$M one-round: ", o)
end
n = 8
Bk = rand(-3:3, n, 3); K = Float64.(Bk*Bk'); Br = rand(-2:2,n,n); R0 = Float64.(Br*Br')
for M in (1, 3, 6)
    println("n=8 M=$M rank3 two: ", run(n, M, K; R0=R0)); println("n=8 M=$M rank3 one: ", run(n, M, K; oneround=true))
end
# brackets are not ordered: spread vs spike
for (name, Kx) in (("ones/n", fill(1.0/n, n, n)), ("e1e1'", Matrix(Diagonal([1.0; zeros(n-1)]))))
    o2 = run(n, 3, Kx); o1 = run(n, 3, Kx; oneround=true)
    println(name, ": two-round E=", o2.Etr, " nys=", o2.nys, " | one-round E=", o1.Etr, " nys=", o1.nys, " | simple=", o2.simple)
end
