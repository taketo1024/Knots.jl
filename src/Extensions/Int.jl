function isunit(a::T) :: Bool where {T<:Integer}
    isone(a) || isone(-a)
end

function normalizing_unit(a::T) :: Tuple{T, T} where {T<:Integer}
    a >= 0 ? (one(T), one(T)) : (-one(T), -one(T))
end