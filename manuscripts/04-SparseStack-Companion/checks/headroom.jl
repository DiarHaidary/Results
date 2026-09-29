# Diagnostic: headroom in SS's own scalar majorant T_q (no row-sum product/envelope loss).
include_string(Main, split(read("diag.jl",String),"for (C,B) in")[1])
function worst(C,B; ds=[1,10,100,10^4,10^6,10^8], qs=[2,5,10,40,100,200], es=[1.0,0.5,0.05])
    w=0.0; pt=nothing
    for d in ds, q in qs, ε in es
        v=2q+1; D=d+v; s=ceil(Int,B*v/ε-1e-12); b=ceil(Int,C*D/(s*ε^2)-1e-12)
        T=Tq(d,q,s,b); y=zeros(BigFloat,q+1); y[1]=1
        for _ in 1:q; y=T*y; end
        r=Float64(sum(y.^2)^(1/(2q))/ε); if r>w; w=r; pt=(d,q,ε,s,b); end
    end
    w,pt
end
for (C,B) in [(69,8.0),(40,8.0),(30,8.0),(25,8.0),(30,4.0),(20,4.0)]
    println((C,B), " => ", worst(C,B))
end
