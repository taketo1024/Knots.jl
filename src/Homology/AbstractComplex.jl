using AbstractAlgebra: RingElement, Ring

abstract type AbstractComplex{R<:RingElement, RR<:Ring} end

function baseRing(C::AbstractComplex{R, RR}) :: RR where {R, RR <: Ring}
    throw(MethodError(baseRomg, (C,)))
end

function hDegRange(C::AbstractComplex) :: UnitRange{Int}
    throw(MethodError(hDegRange, (C,)))
end

function generators(C::AbstractComplex, k::Int) :: Vector{Any}
    throw(MethodError(generators, (C, k)))
end

function differentialDegree(C::AbstractComplex) :: Int
    +1
end

function differential(C::AbstractComplex{R}, k::Int) :: AbstractMatrix{R} where {R<:RingElement}
    throw(MethodError(differential, (C, k)))
end