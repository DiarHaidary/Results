# Editor re-check: rounded counterexample of Remark (Lemma G hypotheses), V_k values, Theorem F bounds.
setprecision(BigFloat, 128)
function Dval(c,u,v,α; h=0.01)
    π_ = c ./ sum(c); ut = u .- sum(π_ .* u); vt = v .- sum(π_ .* v)
    acc = big(0)
    for s in -60:h:80
        t = exp(big(s)); ct = c ./ (1 .+ t .* c); Φ = prod((1 .+ t .* c) .^ (-α))
        acc += Φ * (α^2 * sum(ct .* ut) * sum(ct .* vt) + α * sum(ct .^ 2 .* ut .* vt)) * t
    end
    Float64(acc*h)
end
c = big.([0.02139, 0.27215, 0.03634]); u = big.([-1.6053, -0.1726, -0.0118]); v = big.([2.5635, 2.5023, 1.1242])
println("rounded counterexample D = ", Dval(c,u,v,big(1)), "   (h=0.005: ", Dval(c,u,v,big(1);h=0.005), ")")
# sorted c (monotone) with same u,v should give D<=0
println("monotone c=(0.02139,0.03634,0.27215): D = ", Dval(big.([0.02139,0.03634,0.27215]),u,v,big(1)))
function esym(l)
    e = zeros(Float64, length(l)+1); e[1]=1
    for x in l
        for j in length(l):-1:1; e[j+1] += x*e[j]; end
    end
    e
end
V(l,k) = (k+1)*esym(l)[k+2]/esym(l)[k+1]
for sp in ([5,3,2,1.0],[10,1,0.1,0.01,0.001],[4,4,1,0.5,0.25,0.1],[100,10,1,1,1.0])
    println(sp, "  V_k = ", [round(V(sp,k),sigdigits=6) for k in 1:length(sp)-1])
end
Em(e,m) = sum(e[j+1]/factorial(m-j) for j in 0:min(length(e)-1,m))
for (C,k) in (([1.0,1],2),([1.0,1],4),([1e4,1],2),([1e4,1e4],2),([1e4,1e4],4),([10,1,0.1],3),([1e4,1e4,1e4],3),([1e4,1e4,1e4],5))
    e = esym(C); println("C=",C," k=",k,"  bound=", round((k+1)*Em(e,k+1)/Em(e,k)-1, digits=4), "  r/(k+1-r)=", round(length(C)/(k+1-length(C)),digits=4))
end
