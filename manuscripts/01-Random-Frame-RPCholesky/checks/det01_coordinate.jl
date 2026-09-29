# Exact rational computation of coordinate RPCholesky on the scale-separated matrix A = L D L^T,
# L_ab = h^b/(a+b) (a>b), D = diag(1, η, ..., η^r), η = h^(n^2), n = r+1.  Reports E tr R_r * tr(A^{-1})
# (= E tr R_r / λ_min up to a factor 1+O(η)) and the volume ratio V_r(A) * tr(A^{-1}) / (r+1).
using LinearAlgebra, Printf
function det01(r, h::Rational{BigInt})
    n = r + 1; η = h^(n^2)
    L = Matrix{Rational{BigInt}}(I, n, n)
    for a in 1:n, b in 1:a-1
        L[a, b] = h^b / (a + b)
    end
    D = Diagonal([η^(i) for i in 0:r])
    A = L * D * L'
    minors = Dict{UInt,Rational{BigInt}}()
    idx(mask) = [i for i in 1:n if (mask >> (i - 1)) & 1 == 1]
    for mask in UInt(0):UInt(2^n - 1)
        S = idx(mask)
        minors[mask] = isempty(S) ? one(Rational{BigInt}) : det(A[S, S])
    end
    P = Dict{UInt,Rational{BigInt}}(UInt(0) => 1)
    for s in 0:r-1
        Pn = Dict{UInt,Rational{BigInt}}()
        for (m, p) in P
            dS = minors[m]
            outs = [i for i in 1:n if (m >> (i - 1)) & 1 == 0]
            T = sum(minors[m | (UInt(1) << (i - 1))] for i in outs) / dS
            for i in outs
                m2 = m | (UInt(1) << (i - 1))
                Pn[m2] = get(Pn, m2, 0) + p * (minors[m2] / dS) / T
            end
        end
        P = Pn
    end
    detA = minors[UInt(2^n - 1)]
    er = sum(minors[m] for m in keys(minors) if count_ones(m) == r)
    ErR = sum(p * detA / minors[m] for (m, p) in P)
    ratio = ErR * er / detA            # E tr R_r * tr(A^{-1})
    vol = (r + 1) * detA / er * er / detA  # V_r * tr A^{-1} = r+1 exactly
    return ratio
end
for r in 1:6
    for h in (big(1) // 10, big(1) // 100)
        (r >= 5 && h == big(1)//100) && continue
        t = @elapsed x = det01(r, h)
        @printf("r=%d  h=%s  E tr R_r * tr(A^-1) = %.6f   (2^r = %d, volume value (r+1) = %d)  [%.1fs]\n", r, string(h), Float64(x), 2^r, r + 1, t)
    end
end
