# Energy-biased Gaussian probing in the full space (no rotation, dense non-diagonal A) vs the Haar one-step value.
using LinearAlgebra, Random, Statistics, Printf
function run(λ, k, cplx, N; n_extra=2, seed=1)
    rng = MersenneTwister(seed); n = length(λ) + n_extra
    λf = vcat(λ, zeros(n_extra))
    G = randn(rng, n, n); Q = Matrix(qr(G).Q); A = Q * Diagonal(λf) * Q'; A = (A + A') / 2
    acc = 0.0; acc2 = 0.0
    for _ in 1:N
        R = copy(A)
        for s in 1:k
            # sample w ~ N(0,I) tilted by w'Rw: mixture over eigen-directions of R
            F = eigen(Symmetric(R)); v = max.(F.values, 0); T = sum(v)
            x = rand(rng) * T; i = 1; a = v[1]
            while a < x && i < n; i += 1; a += v[i]; end
            z = cplx ? randn(rng, ComplexF64, n) : complex.(randn(rng, n))
            m2 = cplx ? randexp(rng) + randexp(rng) : randn(rng)^2 + randn(rng)^2 + randn(rng)^2
            z[i] = sqrt(m2) * (cplx ? cis(2π * rand(rng)) : sign(randn(rng)))
            w = F.vectors * z
            Rw = R * w; R = Hermitian(R - Rw * Rw' / real(dot(w, Rw))) |> Matrix
            R = real.(R) .* (cplx ? 1 : 1)
            if cplx; end
        end
        t = real(tr(R)); acc += t; acc2 += t^2
    end
    m = acc / N; m, sqrt((acc2 / N - m^2) / N)
end
m, s = run([5.0, 3, 2, 1], 2, false, 400_000)
@printf("real Gaussian energy-biased probing, k=2: %.5f ± %.5f  (Haar one-step exact 4.41527, volume 4.46341)\n", m, s)
