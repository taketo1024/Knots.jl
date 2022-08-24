import Base: gcdx
using Polynomials

function leadcoeff(a::T) :: R where {R, T<:Polynomial{R}}
    coeffs(a)[end]
end

function isunit(a::T) :: Bool where {T<:Polynomial}
    length(a) == 1 && isunit(leadcoeff(a))
end

function normalizing_unit(a::T) :: Tuple{T, T} where {T<:Polynomial}
    c = leadcoeff(a)
    c == 0 ? (one(T), one(T)) : (T(inv(c)), T(c))
end

function Base.gcdx(f::T, g::T) :: Tuple{T, T, T} where {T<:Polynomial}
    s0, s1 = oneunit(T), zero(T)
    t0, t1 = s1, s0

    # The loop invariant is: s0*a0 + t0*b0 == a
    x = f
    y = g
    while y != 0
        q, r = divrem(x, y)
        x, y = y, r
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end

    (x, s0, t0)
end

function symbol(::Type{T}) :: String where {R, X, T<:Polynomial{R, X}}
    "$(symbol(R))[$X]"
end