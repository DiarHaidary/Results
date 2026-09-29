# One-round covariance map / diagonal-training error check (exact rational, n=4; Float n=8)
using LinearAlgebra, Random
function hadamard(n); H = reshape([1],1,1); while size(H,1) < n; H = [H H; H -H]; end; H; end
subsets(n, M) = [findall(==(1), digits(b, base=2, pad=n)) for b in 0:(2^n-1) if count_ones(b) == M]
signs(n) = [[(b >> i) & 1 == 1 ? -1 : 1 for i in 0:n-1] for b in 0:(2^n-1)]
function chk(n, m, R0, T)
    H = T.(hadamard(n)) ./ (T == Float64 ? sqrt(n) : T(2)); In = Matrix{T}(I,n,n)
    Ec = zeros(T,n,n); Ed = zero(T); c = 0
    for x in signs(n), Ts in subsets(n,m)
        Q = H*Diagonal(T.(x)); P = Q[Ts,:]; G = (T(n)/m) .* (P'*P)
        Ec += (G-In)*R0*(G-In); Ed += sum(abs2, diag(R0*G) - diag(R0)); c += 1
    end
    Ec ./= c; Ed /= c; lam = T(n-m)/(m*(n-1))
    pc = lam*(tr(R0)*In + R0 - 2Diagonal(diag(R0))); pd = lam*(sum(abs2,R0) - sum(abs2,diag(R0)))
    println("n=$n m=$m covdiff=", Float64(maximum(abs.(Ec-pc))), " diagdiff=", Float64(Ed-pd))
end
Random.seed!(3)
Br = rand(-2:2,4,4); R0 = Rational{BigInt}.(Br*Br'); for m in 1:3; chk(4,m,R0,Rational{BigInt}); end
Br = rand(-2:2,8,8); R0 = Float64.(Br*Br'); for m in (1,3,5); chk(8,m,R0,Float64); end
