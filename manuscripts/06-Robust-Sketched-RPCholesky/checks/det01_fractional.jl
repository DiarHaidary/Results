# Diagnostic (not a proof). DET-01 scale-separated family (j=1, n=r+1):
# A = L D L^*, L_ab = h^b/(a+b) (a>b), D = diag(1,eta,...,eta^(r-1),eta^r), eta = h^(n^2).
# Exact enumeration of all RP histories of r pivots; report E T_r/tau, E (T_r/tau)^theta,
# and the union of two independent r-pivot chains: E tr(A/(S1 u S2))/tau.
using LinearAlgebra
setprecision(BigFloat, 4096)
function family(r, h)
    n = r + 1
    L = Matrix{BigFloat}(I, n, n)
    for a in 1:n, b in 1:a-1
        L[a,b] = h^b / (a + b)
    end
    eta = h^(n^2)
    d = [eta^(i-1) for i in 1:r]; push!(d, eta^r)
    A = L * Diagonal(d) * L'
    return A
end
function lambda_min(A)
    B = inv(A); x = ones(BigFloat, size(A,1)); lam = big(0)
    for it in 1:60
        y = B*x; lam = norm(y)/norm(x); x = y/norm(y)
    end
    return 1/lam
end
function enumerate_rp(A, r)
    n = size(A,1)
    setprob = Dict{Vector{Int},BigFloat}(); setT = Dict{Vector{Int},BigFloat}()
    function rec(R, S, p, depth)
        if depth == r
            key = sort(S)
            setprob[key] = get(setprob, key, big(0)) + p
            setT[key] = tr(R)
            return
        end
        T = tr(R)
        for i in 1:n
            i in S && continue
            Rii = R[i,i]; Rii <= 0 && continue
            Rn = R - R[:,i]*R[i,:]'/Rii
            rec(Rn, vcat(S,i), p*Rii/T, depth+1)
        end
    end
    rec(copy(A), Int[], big(1), 0)
    return setprob, setT
end
for r in 2:4, h in (big"1e-2", big"1e-3")
    A = family(r, h); tau = lambda_min(A)
    P, Tv = enumerate_rp(A, r)
    m1 = sum(P[k]*Tv[k] for k in keys(P))/tau
    mh = sum(P[k]*(Tv[k]/tau)^(big(1)/2) for k in keys(P))
    mq = sum(P[k]*(Tv[k]/tau)^(big(1)/4) for k in keys(P))
    uni = sum(P[k]^2*Tv[k] for k in keys(P))/tau   # two independent chains: union = all coords unless same set
    println("r=$r h=$(Float64(h)): E T/tau=$(Float64(m1)) (2^r=$(2^r)),  E(T/tau)^(1/2)=$(Float64(mh)),  E(T/tau)^(1/4)=$(Float64(mq)),  union-of-2 E T/tau=$(Float64(uni)),  bound (E T^(1/2))^2/tau=$(Float64(mh^2))")
end
