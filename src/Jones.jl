using Polynomials

P = LaurentPolynomial{Int, :q}

function jonesPolynomial(l::Link; normalized=false) :: P
    if normalized && isEmpty(l)
        throw(Exception)
    end

    n = crossingNum(l)
    (n₊, n₋) = signedCrossingNums(l)
    s = n₊ - 2n₋

    q = variable(P)
    qinv = P([1], -1)
    
    e = (-1)^isodd(n₋)
    a = (s ≥ 0) ? q^s : qinv^(-s)

    res = zero(P)
    for u in 0 : 2^n - 1
        s = digits(u, base=2, pad=n)
        w = sum(+, s, init=0)
        d = resolve(l, s)
        r = normalized ? length(components(d)) - 1 : length(components(d))
        p = (-q)^w * (q + qinv)^r
        res += p
    end

    J = e * a * res
end