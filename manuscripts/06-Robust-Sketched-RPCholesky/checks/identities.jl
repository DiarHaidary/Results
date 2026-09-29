# Sanity checks of the exact identities used in Section 5 (numerical).
using LinearAlgebra, Random, Printf
Random.seed!(3)
D, N, m, t = 60, 40, 30, 5
F = randn(D, N) * Diagonal(1 ./ (1:N)) ; A = F'F
Pi = randn(m, D) / sqrt(m); B = F' * Pi' * Pi * F
schur(M, J) = M - M[:, J] * (M[J, J] \ M[J, :])
w = [3, 17, 8, 25, 1]; J = w[1:t]
Q = Matrix(qr(F[:, J]).Q)[:, 1:t]; E = (I - Q * Q') * F; G = Q' * Pi' * Pi * Q
P = Pi * Q * (G \ (Q' * Pi'))
RB = schur(B, J); RA = schur(A, J)
@printf("Lemma sk-residual: %.2e\n", norm(RB - (Pi * E)' * (I - P) * (Pi * E)) / norm(RB))
# path likelihood identity
function pathprob(M, w)
    p = 1.0
    for s in 0:length(w)-1
        R = s == 0 ? M : schur(M, w[1:s]); p *= R[w[s+1], w[s+1]] / tr(R)
    end
    p
end
ratio = pathprob(B, w) / pathprob(A, w)
rhs = det(G) * prod(tr(schur(A, w[1:s])) / tr(schur(B, w[1:s])) for s in 1:t-1) * tr(A) / tr(B)
@printf("likelihood identity: ratio=%.6f rhs=%.6f\n", ratio, rhs)
# reuse formula
XB = (Pi * F[:, J]) \ (Pi * F)
lhs = norm(F - F[:, J] * XB)^2; r2 = tr(RA) + norm(G \ (Q' * Pi' * Pi * E))^2
@printf("reuse identity: %.6f vs %.6f\n", lhs, r2)
