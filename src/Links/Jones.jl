using Polynomials

const JonesPolynomial = LaurentPolynomial{Int, :q}

function jonesPolynomial(l::Link; normalized=false) :: JonesPolynomial
    if normalized && isEmpty(l)
        throw(Exception)
    end

    P = JonesPolynomial
    n = crossingNum(l)
    (n₊, n₋) = signedCrossingNums(l)
    s = n₊ - 2n₋

    q = variable(P)
    qinv = P([1], -1)
    
    e = (-1)^isodd(n₋)
    a = (s ≥ 0) ? q^s : qinv^(-s)

    res = reduce(0 : 2^n - 1; init=zero(P)) do res, u
        s = digits(u, base=2, pad=n)
        w = sum(s)
        d = resolve(l, s)
        r = normalized ? length(components(d)) - 1 : length(components(d))
        p = (-q)^w * (q + qinv)^r
        res += p
    end

    e * a * res
end