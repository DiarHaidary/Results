# Exact DP (numerical illustration) of E tr R_K / tau_r on A = diag(H*I_r, I_m) for four pivot rules:
# exact RP, two-sided adversary (a pi <= q <= pi/a), upper-only adversary (q <= pi/a), lower-only adversary (q >= a pi).
# Every adversary gives the heads the least mass its constraint allows.
using Printf
function curve(r, m, H, Kmax, a, mode)
    P = Dict((r, m) => 1.0); out = Float64[]
    for k in 1:Kmax
        Q = Dict{Tuple{Int,Int},Float64}()
        for ((h, b), p) in P
            T = h * H + b
            πh = h * H / T; πu = b / T
            qh = mode == :exact ? πh :
                 mode == :two   ? max(a * πh, 1 - πu / a) :
                 mode == :upper ? max(0.0, 1 - πu / a) :
                 mode == :lower ? a * πh : error()
            h == 0 && (qh = 0.0); b == 0 && (qh = 1.0)
            h > 0 && qh > 0 && (Q[(h - 1, b)] = get(Q, (h - 1, b), 0.0) + p * qh)
            b > 0 && qh < 1 && (Q[(h, b - 1)] = get(Q, (h, b - 1), 0.0) + p * (1 - qh))
        end
        P = Q
        push!(out, sum(p * (h * H + b) for ((h, b), p) in P) / m)
    end
    out
end
for (r, m, a, Kmax) in ((1, 4000, 0.5, 1000), (4, 4000, 0.5, 1000))
    H = (m - Kmax) * (1 - a) / (a * r)
    cs = Dict(md => curve(r, m, H, Kmax, a, md) for md in (:exact, :two, :upper, :lower))
    println("# r=$r m=$m a=$a H=$H   columns: K exact two upper lower  (upper-only closed form (1-K/m)/a at K=Kmax: $((1-Kmax/m)/a))")
    for K in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1000)
        @printf("%d %.4f %.4f %.4f %.4f\n", K, cs[:exact][K], cs[:two][K], cs[:upper][K], cs[:lower][K])
    end
end
