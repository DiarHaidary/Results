lf(r)=Float64(sum(log(big(j)) for j in 1:r; init=big(0.0)))
function B(k,β); b=log(1+2*lf(k)/k); α=min(β,(k-1)/(k*b)); k*b*α-(k-1)*log(α); end
ren(k,δ)=minimum(p->p/(p-1)*(log(1/δ)+B(k,1/p)+k/p*log(2)+(k*log(k)-lf(k))), 1.01:0.01:40)
for δ in (0.01,1e-6)
  k=2; while ren(k,δ) >= log(1/δ)+lf(k); k+=1; end
  println("delta=$δ crossover k=$k ren=",ren(k,δ)," kfact=",log(1/δ)+lf(k))
end
for k in (10,1000); println("k=$k ren=",ren(k,0.01)," kf=",log(100)+lf(k)); end
