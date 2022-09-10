const GaussInteger{T <: Integer} = Complex{T}
const GaussInt = GaussInteger{Int}

function isunit(a::T) :: Bool where {T<:GaussInteger}
    isone(a) || isone(-a) || a == T(im) || a == -T(im)
end

function normalizing_unit(a::T) :: Tuple{T, T} where {T<:GaussInteger}
    if iszero(a)
        (one(T), one(T))
    elseif iszero(a.re)
        imag(a) >= 0 ? (-T(im), T(im)) : (T(im), -T(im))
    else
        (u, uinv) = normalizing_unit(a.im)
        (T(u), T(uinv))
    end
end

function Base.div(x::T, y::T) :: T where {T<:GaussInteger}
    divrem(x, y)[1]
end

function Base.rem(x::T, y::T) :: T where {T<:GaussInteger}
    divrem(x, y)[2]
end

function Base.divrem(a::T, b::T) :: Tuple{T, T} where {I, T<:GaussInteger{I}}
    if iszero(b)
        throw(DivideError())
    end

    z = a / b # float-valued
    q = T(I(round(z.re)), I(round(z.im)))
    r = a - q * b

    return (q, r)
end

function Base.gcdx(a::T, b::T) :: Tuple{T, T, T} where {T<:GaussInteger}
    s0, s1 = oneunit(T), zero(T)
    t0, t1 = s1, s0

    x = a
    y = b
    while y != 0
        q, r = divrem(x, y)
        x, y = y, r
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end

    (x, s0, t0)
end

function computational_weight(a::T) :: Float64 where {I, T<:GaussInteger{I}}
    computational_weight(a.re) + computational_weight(a.im)
end

function symbol(::Type{T}) :: String where {I, T<:GaussInteger{I}}
    "$(symbol(I))[i]"
end