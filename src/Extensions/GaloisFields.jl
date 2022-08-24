using GaloisFields
using GaloisFields: AbstractGaloisField
using ..Utils: subscript

function isunit(a::T) :: Bool where {T<:AbstractGaloisField}
    true
end

function normalizing_unit(a::T) :: Tuple{T, T} where {T<:AbstractGaloisField}
    iszero(a) ? (one(T), one(T)) : (inv(a), a)
end

# workaround for Polynomial divrem
function Base.isapprox(a::T, b::Int) :: Bool where {T<:AbstractGaloisField}
    a == T(b)
end

function Base.gcdx(a::T, b::T) :: Tuple{T, T, T} where {T<:AbstractGaloisField}
    if iszero(a) 
        iszero(b) ? (zero(T), one(T), zero(T)) : (one(T), zero(T), inv(b))
    else
        (one(T), inv(a), zero(T))
    end
end

function symbol(::Type{T}) :: String where {T<:AbstractGaloisField}
    p = char(T)
    "F$(subscript(p))"
end