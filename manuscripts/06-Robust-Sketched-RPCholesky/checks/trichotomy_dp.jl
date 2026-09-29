# Diagnostic: adversarial two-sided gamma-approximate RP vs one-sided (DRF-06) on diag(H*1_r, 1_m).
# State (a,b) = heads left, units left. tau_r = m. Exact DP of E tr R_k.
function etr(r, m, H, K, γ; mode=:two)
    P = Dict((r,m)=>1.0)
    for k in 1:K
        Q = Dict{Tuple{Int,Int},Float64}()
        for ((a,b),p) in P
            T = a*H + b
            if T == 0; Q[(a,b)] = get(Q,(a,b),0.0)+p; continue; end
            πh = a*H/T; πu = b/T
            if mode == :exact
                qh = πh
            elseif mode == :two      # adversary starves heads within two-sided band
                qu = min(1 - γ*πh, πu/γ); qh = 1 - qu
            elseif mode == :one      # one-sided q_i <= π_i/γ : may put zero on heads if allowed
                qu = min(1.0, πu/γ); qh = 1 - qu
            end
            if a > 0 && qh > 0; Q[(a-1,b)] = get(Q,(a-1,b),0.0) + p*qh; end
            if b > 0 && (1-qh) > 0; Q[(a,b-1)] = get(Q,(a,b-1),0.0) + p*(1-qh); end
        end
        P = Q
    end
    sum(p*(a*H+b) for ((a,b),p) in P)/m   # ratio to tau
end
for (r,m,H) in [(1,400,1e6),(4,400,1e6),(8,2000,1e8)]
    println("r=$r m=$m H=$H")
    for K in [r, 2r, 4r, 8r, 16r, 32r]
        e = etr(r,m,H,K,1.0; mode=:exact)
        t = etr(r,m,H,K,0.5; mode=:two)
        o = etr(r,m,H,K,0.5; mode=:one)
        println("  K=$K  exact=$(round(e,digits=4))  two-sided γ=.5: $(round(t,digits=4))  one-sided γ=.5: $(round(o,digits=4))")
    end
end
println("--- DRF-06 regime: H = (m-k)(1-γ)/γ ---")
for (r,m,γ,Kmax) in [(1,4000,0.5,1000),(4,4000,0.5,1000),(1,4000,0.8,1000)]
    H = (m-Kmax)*(1-γ)/γ
    println("r=$r m=$m γ=$γ H=$H")
    for K in [2r, 4r, 8r, 16r, 64r, 250, Kmax]
        e = etr(r,m,H,K,1.0; mode=:exact)
        t = etr(r,m,H,K,γ; mode=:two)
        o = etr(r,m,H,K,γ; mode=:one)
        println("  K=$K  exact=$(round(e,digits=4))  two-sided=$(round(t,digits=4))  one-sided=$(round(o,digits=4))")
    end
end
