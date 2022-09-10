function isunit(a::T) :: Bool where {T<:Rational}
    !iszero(a)
end

function normalizing_unit(a::T) :: Tuple{T, T} where {T<:Rational}
    a == 0 ? (one(T), one(T)) : (inv(a), a)
end

function computational_weight(a::T) :: Float64 where {T<:Rational}
    computational_weight(numerator(a)) + computational_weight(denominator(a))
end

function symbol(::Type{T}) :: String where {T<:Rational}
    "Q"
end