module Env

using AbstractAlgebra
using AbstractAlgebra: Ring, RingElement, RingElem

export set_base_ring

_base_ring = nothing

function init()
end

function current_base_ring()
    _base_ring
end

function set_base_ring(base_ring::Union{Ring, Nothing})
    global _base_ring
    _base_ring = base_ring
end

function check_base_ring(type::Type{R}; warn=true) :: Bool where {R<:RingElement}
    if type <: Union{Integer, Rational}
        true
    elseif !isnothing(_base_ring) && type <: elem_type(_base_ring)
        true
    elseif isnothing(_base_ring)
        warn && @warn "You must explicitly call `Env.set_base_ring()` to enable calculations for AbstractAlgebra.jl types."
        false
    else
        warn && @warn "Base ring `$(_base_ring)` is incompatible with element type `$(R)`."
        false
    end
end

function get_base_ring(type::Type{R}) :: Ring where {R<:RingElement}
    if type <: Union{Integer, Rational}
        parent(zero(R))
    elseif check_base_ring(R)
        _base_ring
    else
        error()
    end
end

function Base.zero(::Type{R}) :: R where {R<:RingElem}
    AbstractAlgebra.zero(_base_ring)
end

function Base.one(::Type{R}) :: R where {R<:RingElem}
    AbstractAlgebra.one(_base_ring)
end

end #module

Env.init()