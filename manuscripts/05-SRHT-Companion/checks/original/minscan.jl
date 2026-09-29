# Exhaustive Float64 scan: is there any feasible M below the bisection value, for any q?
# p_{n,M} = P{Bin(n,M/n)=M} = sqrt(n/(2 pi M (n-M))) * exp(s(n)-s(M)-s(n-M)),  s = Stirling remainder.
s(k) = (x = Float64(k); 1/(12x) - 1/(360x^3) + 1/(1260x^5) - 1/(1680x^7))
function patom64(n, M)
    M < 30 && return Float64(binomial(big(n), M) * (big(M)/n)^M * (1 - big(M)/n)^(n-M))
    sqrt(n / (2pi * M * (n - M))) * exp(s(n) - s(M) - s(n - M))
end
Lam(q, d) = (sqrt(q*(2q+1)) + sqrt(2*(d + 8q^2 + q*(2q-1))))^2
function logcert(n, M, q, d, eps)
    rho = M/n; L = Lam(q, d)
    R = (2*sqrt((1-rho)*L/M) + abs(1-2rho)*L/M) / (1 - patom64(n, M))
    log(d) + 2q*log(R/eps)
end
n = 2^30; eps = 0.5; delta = 0.01
for (d, Mclaim) in ((10, 58396), (100, 125750), (1000, 346315), (10000, 1423395))
    best = (typemax(Int), 0)
    for q in 1:40
        for M in 1:Mclaim
            if logcert(n, M, q, d, eps) <= log(delta) + 1e-9
                if M < best[1]; best = (M, q); end
                break
            end
        end
    end
    println("d=$d: smallest feasible (Float64, tol 1e-9) M=", best[1], " at q=", best[2], "  (bisection value $Mclaim)")
end
