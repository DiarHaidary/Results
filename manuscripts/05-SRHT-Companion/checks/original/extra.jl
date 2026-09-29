# Extra diagnostics for the companion (not proofs)
setprecision(BigFloat, 256)
aq(q) = q*(2q+1); bq(q,d) = (d + 8q^2) + q*(2q-1); Dq(q,d) = d + 12q^2
Lam(q,d) = (sqrt(BigFloat(aq(q))) + sqrt(2*BigFloat(bq(q,d))))^2
# flat: 60 D_{q+1} < 72 Lambda_q for q >= 3, all d
ok = all(((q,d),) -> 60Dq(q+1,d) < 72Lam(q,d), ((q,d) for q in 3:200 for d in (1,2,3,5,10,100,10^4,10^8)))
println("60D_{q+1} < 72Lam_q (flat, q>=3): ", ok)
println("q=1,2 at d=1: ", [(q, 60Dq(q+1,1), Float64(72Lam(q,1))) for q in 1:2])
# Lambda >= 2d+22q^2-q (flat)
println("Lam >= 2d+22q^2-q: ", all(((q,d),)-> Lam(q,d) >= 2d+22q^2-q, ((q,d) for q in 1:100 for d in (1,10,1000))))
# path-local rows: M = ceil(5 Lambda_q(ell)/eps0^2), eps0 = eta/2
for eta in (1//4, 1//8), ell in (11, 51, 101), qq in (:half, :full)
    q = qq == :half ? cld(ell,2) : ell
    e0 = BigFloat(eta)/2
    M = ceil(Int, 5Lam(q, ell)/e0^2)
    println("eta=$eta ell=$ell q=$q M=$M  M/ell=", round(M/ell), "  M/ell^2=", round(M/ell^2, digits=1))
end
# delta/2 vs 8delta^2/9
println("rank-two: ", [(δ, δ/2 > 8δ^2/9) for δ in (0.5, 0.56, 0.5625, 0.57)])
