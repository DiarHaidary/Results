# Independent editor recheck of the prescription table (flat transforms)
setprecision(BigFloat, 256)
function lf(x::Integer)
    x < 2000 && return sum(log(BigFloat(i)) for i in 1:x; init=BigFloat(0))
    X = BigFloat(x)
    X*log(X) - X + log(2*BigFloat(pi)*X)/2 + 1/(12X) - 1/(360X^3) + 1/(1260X^5) - 1/(1680X^7)
end
logatom(n,M) = (rho = BigFloat(M)/n; lf(n-1)-lf(M)-lf(n-1-M) + M*log(rho) + (n-1-M)*log(1-rho))
Lam(q,d) = (a=q*(2q+1); b=d+8q^2+q*(2q-1); (sqrt(BigFloat(a))+sqrt(2*BigFloat(b)))^2)
Dq(q,d) = d+12q^2
function cert(q,M,n,d,eps)
    rho = BigFloat(M)/n; L = Lam(q,d)
    R = (2*sqrt((1-rho)*L/M) + abs(1-2rho)*L/M)/(1-exp(logatom(n,M)))
    d*(R/eps)^(2q)
end
n=2^30; eps=big"0.5"; del=big"0.01"
q4(x) = ceil(Int, log(BigFloat(x))/log(BigFloat(4)) - big"1e-30")
for (d,Mdir,qd,Mmin) in [(10,58406,4,58396),(100,125760,6,125750),(1000,346325,10,346315),(10000,1423405,18,1423395)]
    q = q4(4d/del); qs = q4(6d/(5del)); qf = q4(d/del)
    c216 = ceil(Int,216*Dq(q,d)/eps^2); c72 = ceil(Int,72*Lam(q,d)/eps^2)
    c60 = ceil(Int,60*Dq(qs,d)/eps^2); c20 = ceil(Int,20*Lam(qs,d)/eps^2)
    c72f = ceil(Int,72*Lam(qf,d)/eps^2)
    # check no q<=40 certifies Mmin-1
    anyfeas = any(cert(qq,Mmin-1,n,d,eps) <= del for qq in 1:40)
    println(d," q=",q," q*=",qs," qf=",qf," 216D=",c216," 72L=",c72," 60D=",c60," 20L=",c20, " 72L(qf)=",c72f,
      " cert(dir)=",Float64(cert(qd,Mdir,n,d,eps))," cert(min)=",Float64(cert(qd,Mmin,n,d,eps)),
      " anyq feasible at min-1: ",anyfeas, " cert20L=",Float64(cert(qs,c20,n,d,eps)))
end
# gate / rank-two tables at d=100
for G in [1, 252, factorial(10), binomial(big(1000),10)]
    qs = q4(6*100*G/(5del)); println("Gamma=",Float64(G)," q*=",qs," 20L=",ceil(Int,20*Lam(qs,100)/eps^2)," 60D=",ceil(Int,60*Dq(qs,100)/eps^2))
end
for d in [10,100,1000,10000]
    q1=q4(6d/(5*del/2)); q2=q4(6d/(5*8del^2/9))
    println("rank2 d=",d," ",q1," ",ceil(Int,20*Lam(q1,d)/eps^2)," ",q2," ",ceil(Int,20*Lam(q2,d)/eps^2))
end
println("kappa0=", exp(2-Base.MathConstants.eulergamma)/2)
println("sup check d=1: ", 20*Lam(q4(120),1))
