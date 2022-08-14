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

mutable struct KhChain{R}
    elements::Dict{KhChainGenerator,R}
end

KhChain(elems::AbstractVector, coefs::AbstractVector) = KhChain(Dict(zip(elems, coefs)))

function mapCoeffs(f, c::KhChain{R}) :: KhChain{R} where {R}
    KhChain( Dict( k => f(v) for (k, v) in c.elements ) )
end

function asString(ch::KhChain{R}) where {R}
    if isempty(ch.elements)
        "0"
    else
        join( map( p -> "$(p[2])$(p[1])", collect(ch.elements) ), " + ")
    end
end

function simplify!(ch::KhChain{R}) where {R}
    for (k, v) in ch.elements
        v == zero(R) && delete!(ch.elements, k)
    end
    ch
end

Base.zero(::Type{KhChain{R}}) where {R} = KhChain(Dict{KhChainGenerator, R}())
Base.length(ch::KhChain{R}) where {R} = length(ch.elements)
Base.copy(ch::KhChain{R}) where {R} = KhChain{R}(ch.dim, copy(ch.elements))
Base.keys(ch::KhChain{R}) where {R} = keys(ch.elements)
Base.values(ch::KhChain{R}) where {R} = values(ch.elements)
Base.getindex(ch::KhChain{R}, k::KhChainGenerator) where {R} = get(ch.elements, k, zero(R))
Base.setindex!(ch::KhChain{R}, c::R, k::KhChainGenerator) where {R} = setindex!(ch.elements, c, k)
Base.iterate(ch::KhChain{R}, state...) where {R} = begin
    y = iterate(ch.elements, state...)
    y === nothing && return nothing
    return (y[1], y[2])
end
Base.show(io::IO, ch::KhChain{R}) where {R} = print(io, asString(ch))

Base.:(+)(c1::KhChain{R}, c2::KhChain{R}) where {R} = KhChain(mergewith(+, c1.elements, c2.elements))
Base.:(-)(c::KhChain{R}) where {R} = mapCoeffs(-, c)
Base.:(-)(c1::KhChain{R}, c2::KhChain{R}) where {R} = c1 + (-c2)
Base.:(*)(r::R, c::KhChain{R}) where {R} = mapCoeffs( v -> r * v, c )
Base.:(*)(c::KhChain{R}, r::R) where {R} = mapCoeffs( v -> v * r, c )
Base.:(==)(c1::KhChain{R}, c2::KhChain{R}) where {R} = c1.elements == c2.elements