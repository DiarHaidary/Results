# Diagnostic: large-q behaviour of SS's scalar majorant (tridiagonal recursion, BigFloat).
setprecision(128)
function ratio(d,q,ε,C,B)
    v=2q+1; D=d+v; s=big(ceil(Int,B*v/ε-1e-12)); b=big(ceil(Int,C*D/(Float64(s)*ε^2)-1e-12)); m=s*b; h=b-1; d=big(d)
    α(ν) = ν==0 ? sqrt(BigFloat(d+1)/m) : (sqrt(BigFloat((m+ν)*(d+ν+1))) + sqrt(BigFloat(h*ν*(d+ν-1))/2) + ν*sqrt(BigFloat(h)) + sqrt(BigFloat(b*h*ν*(ν-1))/2))/m
    χ(ν) = ν==0 ? big(0.0) : (d+2ν-1 + 2*sqrt(BigFloat(b*(ν-1)*(d+ν-1))) + (ν+1)*sqrt(BigFloat(b)) + h*ν)/BigFloat(m)
    a=[α(2k) for k in 0:q-1]; c=[χ(2k) for k in 0:q]
    y=zeros(BigFloat,q+1); y[1]=1
    for _ in 1:q
        z=similar(y)
        for k in 1:q+1
            z[k]=c[k]*y[k] + (k>1 ? a[k-1]*y[k-1] : 0) + (k<=q ? a[k]*y[k+1] : 0)
        end
        y=z
    end
    Float64(sum(y.^2)^(1/(2q))/ε)
end
for (C,B) in [(69,8.0),(35,8.0)], d in [10,1000,10^6,10^9], q in [200,1000,3000]
    println((C,B,d,q), " => ", round(ratio(d,q,1.0,C,B),digits=4))
end
