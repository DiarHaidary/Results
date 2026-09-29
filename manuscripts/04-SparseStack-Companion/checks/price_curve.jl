# Data for the price-comparison figure: log(1/delta_0) in nats for the factorial gate
# and for the Renyi-3/2 route (published envelope constant B_{k,beta}), delta = 0.01.
lf(k) = k < 2 ? 0.0 : (k <= 10^6 ? sum(log, 1:k) : k*log(k) - k + 0.5*log(2pi*k))
bk(k) = log(1 + 2*lf(k)/k)
Bk(k, β) = (b = bk(k); α = min(β, (k-1)/(k*b)); k*b*α - (k-1)*log(α))
frakd(k; β=2/3) = 3*Bk(k,β) + 2k*log(2) + 3*(k*log(k) - lf(k))
δ = 0.01
for k in [2,5,10,20,50,100,300,1000,3000,10^4,3*10^4,10^5,1.3*10^5,2*10^5,10^6]
    k = round(Int,k)
    f = log(1/δ) + lf(k); r = 3log(1/δ) + frakd(k)
    println(k, "  fact=", round(f,digits=2), "  renyi=", round(r,digits=2), "  ratio=", round(r/f,digits=3))
end
