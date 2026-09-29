# Numerical diagnostics (Monte Carlo, not proofs) for Sections on Nystrom bounds.
using LinearAlgebra, Random, Statistics
Random.seed!(20260929)
cs(n, b) = (H = zeros(b, n); for i in 1:n; H[rand(1:b), i] = rand((-1.0, 1.0)); end; H)
sparsestack(n, s, b) = vcat([cs(n, b) for _ in 1:s]...) / sqrt(s)
nys_res(A, Om) = (Y = A * Om'; A - Y * pinv(Om * Y) * Y')      # Om: m x n
sym(X) = (X + X') / 2
# 1. one-block drift versus the block bound tr(R^2)/(kappa_b tr R)
function drift_ratio(R, b; N = 20000)
    dr = mean(tr(R) - tr(sym(nys_res(R, cs(size(R, 1), b)))) for _ in 1:N)
    bd = tr(R * R) / ((1 + 2 / b) * tr(R))
    (round(dr, digits = 3), round(bd, digits = 3), round(dr / bd, digits = 2))
end
W12 = (G = randn(12, 12); G * G' / 12)
u = normalize(randn(12)); R1 = u * u'
println("drift A=I12,b=3:      ", drift_ratio(Matrix(1.0I, 12, 12), 3))
println("drift Wishart12,b=3:  ", drift_ratio(W12, 3))
println("drift A=I40,b=8:      ", drift_ratio(Matrix(1.0I, 40, 40), 8))
println("drift rank-1 n=12,b=3:", drift_ratio(R1, 3))
# 2. exact identity E tr(AZAZ) and pointwise ||R||_F^2 <= tr(AZAZ)
function frob_check(A, s, b; N = 20000)
    n = size(A, 1); m = s * b; acc = 0.0; accR = 0.0; viol = 0
    for _ in 1:N
        P = sparsestack(n, s, b); Z = P' * P - I
        t = tr(A * Z * A * Z); r = norm(nys_res(A, P))^2
        acc += t; accR += r; viol += (r > t * (1 + 1e-9) + 1e-9)
    end
    exact = (tr(A)^2 + norm(A)^2 - 2 * sum(diag(A) .^ 2)) / m
    (MC = round(acc / N, digits = 3), exact = round(exact, digits = 3), bound2T2m = round(2 * tr(A)^2 / m, digits = 3),
     ER2 = round(accR / N, digits = 4), blockbound = round((1 + 2 / b) * tr(A)^2 / s, digits = 3), violations = viol)
end
G = randn(30, 30); A30 = G * Diagonal(0.7 .^ (0:29)) * G' / 30
println("frob A30 (decaying), s=4,b=6: ", frob_check(A30, 4, 6))
println("frob Wishart12,     s=3,b=4: ", frob_check(W12, 3, 4))
println("frob I_20,          s=2,b=5: ", frob_check(Matrix(1.0I, 20, 20), 2, 5))
# 3. diagonal training MSE
function diag_check(Rm, s, b; N = 40000)
    n = size(Rm, 1); m = s * b
    mse = mean(sum(abs2, diag(Rm * (P = sparsestack(n, s, b); P' * P)) - diag(Rm)) for _ in 1:N)
    (MC = round(mse, digits = 4), exact = round((norm(Rm)^2 - sum(diag(Rm) .^ 2)) / m, digits = 4))
end
println("diag training W12, s=3,b=4: ", diag_check(W12, 3, 4))
