using LinearAlgebra, Random, Statistics, Printf

# elementary symmetric polynomials e_0..e_d (BigFloat for safety)
function esym(λ)
    d = length(λ); e = zeros(BigFloat, d + 1); e[1] = 1
    for x in λ, k in d:-1:1
        e[k+1] += big(x) * e[k]
    end
    return e   # e[j+1] = e_j
end
Vk(λ, k) = (E = esym(λ); E[k+1] == 0 ? 0.0 : Float64((k + 1) * E[k+2] / E[k+1]))

# c_i, u_i, rho_i for the one-step inequality at budget k (positive spectrum, length >= k)
function coeffs(λ, k)
    d = length(λ); c = zeros(d); u = zeros(d); ρ = zeros(d)
    for i in 1:d
        Ei = esym(λ[setdiff(1:d, i)])
        c[i] = Float64(λ[i] * Ei[k]); u[i] = Float64(1 / Ei[k]); ρ[i] = Float64(Ei[k+1] / Ei[k])
    end
    return c, u, ρ
end

# ∫_0^∞ f(t) dt via t = exp(s), trapezoid (integrands decay exponentially in s at both ends)
function lapint(f; lo = -90.0, hi = 90.0, h = 0.01)
    acc = 0.0
    for s in lo:h:hi
        t = exp(s); acc += f(t) * t
    end
    return acc * h
end

# One-step quantities: returns (LHS = E_tilt V_{k-1}(M_u), RHS = V_k, D, Term1, Term2)
function onestep(λ, k, α)
    c, u, ρ = coeffs(λ, k)
    π_ = c ./ sum(c); ub = sum(π_ .* u); ρb = sum(π_ .* ρ)
    ut = u .- ub; ρt = ρ .- ρb
    Φ(t) = exp(-α * sum(log1p.(t .* c)))
    T1 = lapint(t -> (ct = c ./ (1 .+ t .* c); Φ(t) * α^2 * sum(ct .* ut) * sum(ct .* ρt)))
    T2 = lapint(t -> (ct = c ./ (1 .+ t .* c); Φ(t) * α * sum(ct .^ 2 .* ut .* ρt)))
    D = T1 + T2
    ES = α * sum(λ)                     # E[S] = α Σ c_i u_i = α Σ λ_i
    ESZ = D + ub * ρb * α * sum(c)      # E[SZ] = D + ū ρ̄ E[R]
    return k * ESZ / ES, k * ρb, D, T1, T2
end
