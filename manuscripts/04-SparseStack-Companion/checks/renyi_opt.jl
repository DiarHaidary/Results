# Renyi-order optimisation of the likelihood price versus the factorial price.
lf(k) = sum(log, 1:k; init=0.0)
bk(k) = log(1 + 2*lf(k)/k)
Bk(k,β) = (b=bk(k); α=min(β,(k-1)/(k*b)); k*b*α - (k-1)*log(α))
# log(1/delta0) needed so that the selected frame fails w.p. <= delta
function renyi_cost(k, β, δ)
    p = 1/β
    logM = p*Bk(k,β) + k*log(2) + p*(k*log(k) - lf(k))       # log E_nu L^p
    (p/(p-1))*log(1/δ) + logM/(p-1)
end
fact_cost(k, δ) = log(1/δ) + lf(k)
βs = range(0.01, 0.99, length=981)
for δ in (0.01,)
  for k in [2,5,10,20,50,100,200,500,1000,10^4,10^5,2*10^5]
    best = minimum(renyi_cost(k,β,δ) for β in βs); bβ = βs[argmin([renyi_cost(k,β,δ) for β in βs])]
    println("k=$k  factorial=$(round(fact_cost(k,δ),digits=1))  renyi(p=3/2)=$(round(renyi_cost(k,2/3,δ),digits=1))  renyi(opt)=$(round(best,digits=1)) at p=$(round(1/bβ,digits=2))")
  end
end
# crossover with optimised order
k = 2
while minimum(renyi_cost(k,β,0.01) for β in βs) >= fact_cost(k,0.01); global k += 1; end
println("optimised-order Renyi first beats factorial at k = ", k)
# parameters at k=10 with optimised order
d=100; ε=0.5
c = minimum(renyi_cost(10,β,0.01) for β in βs)
q = ceil(Int, (log(d) + c)/log(4))
for (C,B) in [(81,8.5),(69,8.0)]
  v=2q+1; D=d+v; s=ceil(Int,B*v/ε-1e-12); b=ceil(Int,C*D/(s*ε^2)-1e-12); println((C,B)," opt-order q=$q s=$s b=$b m=$(s*b)")
end
