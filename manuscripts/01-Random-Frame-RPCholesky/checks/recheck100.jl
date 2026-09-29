include("common.jl")
println(onestep([100.0,10,1,1,1], 2, 0.5)[1:2])
# finer quadrature
c,u,ρ = coeffs([100.0,10,1,1,1],2); println(lapint(t->exp(-0.5*sum(log1p.(t.*c))); h=0.001) , " vs ", lapint(t->exp(-0.5*sum(log1p.(t.*c)))))
