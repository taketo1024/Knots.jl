# KhEnhancedState
# Generators of Cube(D) and CKh(D)

struct KhChainGenerator
    state::State
    label::Vector{KhAlgGenerator}
end

function asString(x::KhChainGenerator) 
    isempty(x.label) ? "1" : "$(join(x.label, ""))$(join(map(b -> (b == 0) ? "₀" : "₁", x.state), ""))"
end

Base.show(io::IO, x::KhChainGenerator) = print(io, asString(x))
Base.hash(x::KhChainGenerator) = Base.hash(x.state, Base.hash(x.label))
Base.:(==)(x1::KhChainGenerator, x2::KhChainGenerator) = begin 
    (x1.state, x1.label) == (x2.state, x2.label)
end

# KhChain
# elements of Cube(D) and CKh(D)

using ComputationalHomology

mutable struct KhChain{R} <: ComputationalHomology.AbstractChain{KhChainGenerator,R}
    dim::Int
    elements::Dict{KhChainGenerator,R}
end

KhChain(d::Int, elems::AbstractVector, coefs::AbstractVector) = KhChain(d, Dict(zip(elems, coefs)))
KhChain(elems::AbstractVector, coefs::AbstractVector) = KhChain(0, elems, coefs)

asString(ch::KhChain{R}) where {R} = begin
    if isempty(ch.elements)
        "0"
    else
        join( map( p -> "$(p[2])$(p[1])", collect(ch.elements) ), " + ")
    end
end

# implement interface
Base.zero(::Type{KhChain{R}}) where {R} = KhChain(0, Dict{KhChainGenerator, R}())
Base.length(ch::KhChain{R}) where {R} = length(ch.elements)
Base.copy(ch::KhChain{R}) where {R} = KhChain{R}(ch.dim, copy(ch.elements))
Base.keys(ch::KhChain{R}) where {R} = keys(ch.elements)
Base.values(ch::KhChain{R}) where {R} = values(ch.elements)
Base.getindex(ch::KhChain{R}, k::KhChainGenerator) where {R} = get(ch.elements, k, zero(R))
Base.setindex!(ch::KhChain{R}, c::R, k::KhChainGenerator) where {R} = setindex!(ch.elements, c, k)

Base.push!(ch::KhChain{R}, e::Pair{KhChainGenerator, R}) where {R} = begin
    (k,v) = e
    if k ∈ ch
        ch.elements[k] += v
    else
        push!(ch.elements, e)
    end
end

Base.iterate(ch::KhChain{R}, state...) where {R} = begin
    y = iterate(ch.elements, state...)
    y === nothing && return nothing
    return (y[1], y[2])
end

Base.show(io::IO, ch::KhChain{R}) where {R} = print(io, asString(ch))

Base.:(==)(c1::KhChain{R}, c2::KhChain{R}) where {R} = begin 
    (c1.dim, c1.elements) == (c2.dim, c2.elements)
end

ComputationalHomology.dim(ch::KhChain{R}) where {R} = ch.dim
ComputationalHomology.simplify(ch::KhChain{R}) where {R} = begin
    for (k,v) in ch.elements
        v == zero(R) && delete!(ch.elements, k)
    end
    return ch
end