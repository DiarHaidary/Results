# check the imported prefix asymptotics  η^s / T_s(S) -> 1{s+1 ∉ S}  (T_s = residual trace after pivoting S, |S| = s < n-1)
using LinearAlgebra, Printf
function build(n, h)
    η = h^(n^2); L = Matrix{Rational{BigInt}}(I, n, n)
    for a in 1:n, b in 1:a-1; L[a, b] = h^b / (a + b); end
    L * Diagonal([η^i for i in 0:n-1]) * L', η
end
for n in (3, 4, 5), h in (big(1)//10, big(1)//1000)
    A, η = build(n, h); worst = 0.0
    for mask in 0:2^n-1
        S = [i for i in 1:n if (mask >> (i-1)) & 1 == 1]; s = length(S); s >= n - 1 && continue
        Sc = setdiff(1:n, S)
        T = isempty(S) ? tr(A) : tr(A[Sc, Sc] - A[Sc, S] * (A[S, S] \ A[S, Sc]))
        val = Float64(η^s / T); target = (s + 1) in S ? 0.0 : 1.0
        worst = max(worst, abs(val - target))
    end
    @printf("n=%d h=%s  max |η^s/T_s(S) - 1{s+1∉S}| = %.3e\n", n, string(h), worst)
end
