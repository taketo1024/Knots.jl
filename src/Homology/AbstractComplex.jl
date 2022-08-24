abstract type AbstractComplex{R} end

function hDegRange(C::AbstractComplex) :: UnitRange{Int}
    throw(MethodError(hDegRange, (C,)))
end

function generators(C::AbstractComplex, k::Int) :: Vector{Any}
    throw(MethodError(generators, (C, k)))
end

function differentialDegree(C::AbstractComplex) :: Int
    +1
end

function differential(C::AbstractComplex{R}, k::Int) :: AbstractMatrix{R} where {R}
    throw(MethodError(differential, (C, k)))
end