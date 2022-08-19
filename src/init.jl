module Env

using AbstractAlgebra
using AbstractAlgebra: Ring, RingElement, RingElem

export set_base_ring

_base_ring = nothing

function init()
end

function set_base_ring(base_ring::Union{Ring, Nothing})
    global _base_ring
    _base_ring = base_ring
end

function check_base_ring(type::Type{R}) :: Bool where {R<:RingElement}
    if type <: Union{Integer, Rational}
        true
    elseif !isnothing(_base_ring) && type <: elem_type(_base_ring)
        true
    elseif isnothing(_base_ring)
        @warn "You must explicitly call `Env.set_base_ring()` to enable calculations for AbstractAlgebra.jl types."
        false
    else
        @warn "Base ring `$(_base_ring)` is incompatible with element type `$(R)`."
        false
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