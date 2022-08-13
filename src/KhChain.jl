using ComputationalHomology: AbstractChain

import Base: length, copy, keys, values, getindex, setindex!, push!, iterate

# KhEnhancedState
# Generators of Cube(D) and CKh(D)

struct KhEnhancedState
    state::BitArray
    label::Vector{KhAlgGenerator}
end

asString(x::KhEnhancedState) = "$(join(x.label, ""))$(join(map(b -> b ? "₀" : "₁", x.state), ""))"
Base.show(io::IO, x::KhEnhancedState) = print(io, asString(x))

# KhChain
# Elements of Cube(D) and CKh(D)

mutable struct KhChain{R} <: AbstractChain{KhEnhancedState,R}
    dim::Int
    elements::Dict{KhEnhancedState,R}
end

KhChain(d::Int, elems::AbstractVector, coefs::AbstractVector) = KhChain(d, Dict(zip(elems, coefs)))
KhChain(elems::AbstractVector, coefs::AbstractVector) = KhChain(0, elems, coefs)

# implement interface
dim(ch::KhChain{R}) where {R} = ch.dim
length(ch::KhChain{R}) where {R} = length(ch.elements)
copy(ch::KhChain{R}) where {R} = KhChain{R}(ch.dim, copy(ch.elements))
keys(ch::KhChain{R}) where {R} = keys(ch.elements)
values(ch::KhChain{R}) where {R} = values(ch.elements)
getindex(ch::KhChain{R}, k::KhAlgGenerator) where {R} = get(ch.elements, k, zero(R))
setindex!(ch::KhChain{R}, c::R, k::KhAlgGenerator) where {R} = setindex!(ch.elements, c, k)

function push!(ch::KhChain{R}, e::Pair{KhAlgGenerator,R}) where {R}
    (k,v) = e
    if k ∈ ch
        ch.elements[k] += v
    else
        push!(ch.elements, e)
    end
    ch
end

function iterate(ch::KhChain{R}, state...) where {R}
    y = iterate(ch.elements, state...)
    y === nothing && return nothing
    return (y[1], y[2])
end

function simplify(ch::KhChain)
    for (k,v) in ch.elements
        v == zero(R) && delete!(ch.elements, k)
    end
    return ch
end