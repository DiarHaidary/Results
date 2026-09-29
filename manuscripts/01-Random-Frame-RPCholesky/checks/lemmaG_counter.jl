# Recheck the non-monotone-c violation of Lemma G with finer quadrature (alpha = 1).
setprecision(BigFloat, 128)
c = big.([0.02139032248858256, 0.272153109666727, 0.036338982346237496])
u = big.([-1.6053297579356425, -0.17257653131395223, -0.011796877238440431])
v = big.([2.5635013365259156, 2.502265905665405, 1.1242039214144768])
α = big(1)
π_ = c ./ sum(c); ut = u .- sum(π_ .* u); vt = v .- sum(π_ .* v)
function integ(h)
    acc = big(0)
    for s in -60:h:80
        t = exp(big(s)); ct = c ./ (1 .+ t .* c); Φ = prod((1 .+ t .* c) .^ (-α))
        acc += Φ * (α^2 * sum(ct .* ut) * sum(ct .* vt) + α * sum(ct .^ 2 .* ut .* vt)) * t
    end
    acc * h
end
for h in (0.02, 0.01, 0.005)
    println("h=$h  D = ", Float64(integ(big(h))))
end
