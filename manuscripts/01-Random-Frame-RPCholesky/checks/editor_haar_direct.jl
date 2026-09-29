# Editor re-check: direct Haar-rotated RP (not the spectral chain), exact expectation over pivots, spectrum (5,3,2,1), k=2.
using LinearAlgebra, Random
Random.seed!(7)
function haar(n, T)
    G = randn(T, n, n); F = qr(G); Q = Matrix(F.Q); d = diag(F.R); Q * Diagonal(d ./ abs.(d))
end
schur(B,i) = B - B[:,i]*transpose(B[i,:])/B[i,i]
function E2(B)
    t = real(tr(B)); acc = 0.0
    for i in 1:size(B,1)
        pi_ = real(B[i,i])/t; B1 = schur(B,i); t1 = real(tr(B1))
        for j in 1:size(B,1)
            j==i && continue
            acc += pi_*real(B1[j,j])/t1*real(tr(schur(B1,j)))
        end
    end
    acc
end
Λ = Diagonal([5.0,3,2,1])
for T in (Float64, ComplexF64)
    N = 400_000; s = 0.0; s2 = 0.0
    for _ in 1:N
        Q = haar(4,T); x = E2(Q*Λ*Q'); s += x; s2 += x^2
    end
    m = s/N; se = sqrt((s2/N - m^2)/N)
    println(T, "  E tr R_2 = ", round(m,digits=5), " ± ", round(se,digits=5))
end
