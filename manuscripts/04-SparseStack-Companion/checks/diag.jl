# Diagnostic (not a proof): build SS's scalar majorant T_q from alpha_nu, chi_nu with the
# new (C,B) prescriptions and check (T^{2q})_00^{1/2q} < eps/2 and row sums <= eps*Phi_{t,b}(u_k).
setprecision(200)
function Tq(d,q,s,b); d=big(d); q=big(q); s=big(s); b=big(b)
    m = s*b; h = b-1
    α(ν) = ν==0 ? sqrt(BigFloat(d+1)/m) :
        (sqrt(BigFloat((m+ν)*(d+ν+1))) + sqrt(BigFloat(h*ν*(d+ν-1))/2) + ν*sqrt(BigFloat(h)) + sqrt(BigFloat(b*h*ν*(ν-1))/2))/m
    χ(ν) = ν==0 ? big(0.0) :
        (d+2ν-1 + 2*sqrt(BigFloat(b*(ν-1)*(d+ν-1))) + (ν+1)*sqrt(BigFloat(b)) + h*ν)/BigFloat(m)
    T = zeros(BigFloat,Int(q)+1,Int(q)+1)
    for k in 0:Int(q)
        T[k+1,k+1] = χ(2k)
        if k<q; T[k+1,k+2]=α(2k); T[k+2,k+1]=α(2k); end
    end
    T
end
function Phi(t,b,u,C,B)
    θ = 1-t+t*u; G=(2sqrt(1-1/b)+1)/sqrt(b)
    2/sqrt(C)*sqrt((1+t*u/C)*θ) + (2+sqrt(2))/sqrt(C*B)*sqrt(u*θ) + θ/C + (1+sqrt(2)+G)*u/B
end
for (C,B) in [(69,8.0),(74,7.5),(79,7.0),(81,8.5)]
    worst = 0.0; worstpt = nothing; rsviol = 0
    for d in [1,2,3,5,10,30,100,1000,10^4,10^6], q in [2,3,4,6,10,20,40], ε in [1.0,0.9,0.5,0.1,0.013]
        v=2q+1; D=d+v
        s = ceil(Int, B*v/ε - 1e-12)
        b = ceil(Int, C*D/(s*ε^2) - 1e-12)
        T = Tq(d,q,s,b)
        x = zeros(BigFloat,q+1); x[1]=1
        y = copy(x); for _ in 1:q; y = T*y; end
        val = sum(y.^2)                      # (T^{2q})_00
        r = Float64(val^(1/(2q))/ε)
        if r>worst; worst=r; worstpt=(d,q,ε,s,b); end
        t = v/D
        for k in 1:q
            rs = Float64(sum(T[k+1,:]))/ε
            if rs > Phi(t,b,(2k+1)/v,C,B)*(1+1e-12); rsviol+=1; end
        end
    end
    println("(C,B)=($C,$B): max (T^{2q})_00^{1/2q}/eps = ", round(worst,digits=5), " at ", worstpt, "; row-sum>Phi violations: ", rsviol)
end
