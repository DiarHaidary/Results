# Exact certificate for the interval envelopes Psi_{l,r}: every radical is
# rounded up to a rational with denominator K; everything else is exact.
const K = big(10)^12
function sqrt_up(x::Rational{BigInt})              # a rational >= sqrt(x)
    N, D = numerator(x) * K^2, denominator(x)
    k = isqrt(cld(N, D))
    while k * k * D < N; k += 1; end
    return k // K
end
const S2 = sqrt_up(big(2) // 1)
G_up(n) = (2 * sqrt_up((n - 1) // n) + 1) * sqrt_up(1 // n)
function Psi_up(u, l, r, C, B, n)
    th = 1 - l + l * u
    2 * sqrt_up((1 + r * u / C) * th / C) + (2 + S2) * sqrt_up(u * th / (C * B)) +
        th / C + (1 + S2 + G_up(n)) * u / B
end
bucket(C, B, r) = max(floor(BigInt, C / B), floor(BigInt, C / (r * (B + 1 // 5)))) + 1
intervals = [[(big(j) // 40, big(j + 1) // 40) for j in 0:19]; (big(1) // 2, big(1) // 1)]
pairs = [(big(69) // 1, big(8) // 1), (big(74) // 1, big(15) // 2), (big(79) // 1, big(7) // 1)]
nodes = [big(2j + 1) // 16 for j in 0:7]
ok = true
for (C, B) in pairs, (l, r) in intervals
    n = bucket(C, B, r)
    P = prod(Psi_up(u, l, r, C, B, n) for u in nodes)   # >= prod of Psi(u_j)
    lo, hi = big(0), big(500000)
    while hi - lo > 1
        mid = (lo + hi) >> 1
        if mid^8 >= P * big(10)^48; hi = mid; else; lo = mid; end
    end
    T = hi                    # least T with (T/1e6)^8 >= P, if P < 1/256
    global ok &= (P < big(1) // 256) && (T < 500000) && (9C < 800)
    println("C=$C B=$B [$l,$r] n=$n T=$T")
end
println("certificate valid: ", ok)
