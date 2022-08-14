# KhCubeVertex

struct KhCubeVertex
    state::State
    circles::Vector{Component}
    generators::Vector{KhChainGenerator}
end

function KhCubeVertex(l::Link, s::State) :: KhCubeVertex
    circles = components(resolve(l, s))
    r = length(circles)
    generators = map(0 : 2^r - 1) do i in
        bits = digits(i, base=2, pad=r)
        label = map(b -> (b == 0) ? X : I, bits)
        KhChainGenerator(s, label)
    end
    KhCubeVertex(s, circles, generators)
end

# KhCubeEdge

abstract type KhCubeEdgeType end

struct mergeEdge <: KhCubeEdgeType
    sign::Int
    from::Tuple{Int, Int}
    to::Int
end

struct splitEdge <: KhCubeEdgeType
    sign::Int
    from::Int
    to::Tuple{Int, Int}
end

# KhCube

mutable struct KhCube{R} 
    structure::KhAlgStructure{R}
    link::Link
    vertices::Dict{State, KhCubeVertex}
    edges::Dict{Tuple{State, State}, Union{mergeEdge, splitEdge}}
end

function KhCube(str::KhAlgStructure{R}, l::Link) where {R} 
     KhCube(str, l, Dict{State, KhCubeVertex}(), Dict{Tuple{State, State}, Union{mergeEdge, splitEdge}}())
end

function dim(cube::KhCube{R}) where {R} 
    crossingNum(cube.link)
end

function vertex(cube::KhCube{R}, u::State) :: KhCubeVertex where {R}
    @assert length(u) == dim(cube)

    get!(cube.vertices, u) do 
        KhCubeVertex(cube.link, u)
    end
end

function nextVertices(cube::KhCube{R}, u::State) :: Vector{State} where {R}
    @assert length(u) == dim(cube)

    n = length(u)
    indices = filter(i -> u[i] == 0, 1 : n)
    map(indices) do i in
        v = copy(u)
        v[i] = 1
        v
    end
end

Base.findfirst(arr::AbstractArray{T, N}, elm::T) where {T, N} = findfirst( x -> x == elm, arr )

function edgeSign(cube::KhCube{R}, u::State, v::State) :: Int where {R}
    @assert length(u) == length(v) == dim(cube)
    @assert sum(u) + 1 == sum(v)

    n = dim(cube)
    i = findfirst(i -> u[i] != v[i], 1 : n)
    k = count(1 : i) do i
        u[i] == v[i] == 1
    end
    (-1)^isodd(k)
end

function edge(cube::KhCube{R}, u::State, v::State) :: Union{mergeEdge, splitEdge} where {R}
    get!(cube.edges, (u, v)) do 
        _edge(cube, u, v)
    end
end

function _edge(cube::KhCube{R}, u::State, v::State) :: Union{mergeEdge, splitEdge} where {R}
    @assert length(u) == length(v) == dim(cube)
    @assert sum(u) + 1 == sum(v)

    e = edgeSign(cube, u, v)
    Cᵤ = vertex(cube, u).circles
    Cᵥ = vertex(cube, v).circles
    cᵤ = filter(c -> c ∉ Cᵥ, Cᵤ)
    cᵥ = filter(c -> c ∉ Cᵤ, Cᵥ)

    return if (length(cᵤ), length(cᵥ)) == (2, 1)
        i₁ = findfirst(Cᵤ, cᵤ[1])
        i₂ = findfirst(Cᵤ, cᵤ[2])
        j  = findfirst(Cᵥ, cᵥ[1])
        if i₁ > i₂ 
            (i₁, i₂) = (i₂, i₁)
        end
        mergeEdge(e, (i₁, i₂), j)

    elseif (length(cᵤ), length(cᵥ)) == (1, 2)
        i  = findfirst(Cᵤ, cᵤ[1])
        j₁ = findfirst(Cᵥ, cᵥ[1])
        j₂ = findfirst(Cᵥ, cᵥ[2])
        if j₁ > j₂ 
            (j₁, j₂) = (j₂, j₁)
        end
        splitEdge(e, i, (j₁, j₂))

    else
        throw(Exception)
    end
end

function edgeMap(cube::KhCube{R}, u::State, v::State, x::KhChainGenerator) :: Vector{Tuple{KhChainGenerator, R}} where {R}
    @assert length(u) == length(v) == dim(cube)
    @assert sum(u) + 1 == sum(v)

    edg = edge(cube, u, v)
    if isa(edg, mergeEdge)
        _mergeEdgeMap(cube, edg, v, x)
    else
        _splitEdgeMap(cube, edg, v, x)
    end
end

function _mergeEdgeMap(cube::KhCube{R}, edg::mergeEdge, v::State, x::KhChainGenerator) :: Vector{Tuple{KhChainGenerator, R}} where {R}
    m = product(cube.structure) 
        
    (e, (i, j), k) = (edg.sign, edg.from, edg.to)
    (xᵢ, xⱼ) = (x.label[i], x.label[j])

    res = map(m(xᵢ, xⱼ)) do (yₖ, r) 
        label = copy(x.label)
        deleteat!(label, j)
        deleteat!(label, i)
        insert!(label, k, yₖ)
        y = KhChainGenerator(v, label)
        (y, e * r)
    end
end

function _splitEdgeMap(cube::KhCube{R}, edg::splitEdge, v::State, x::KhChainGenerator) :: Vector{Tuple{KhChainGenerator, R}} where {R}
    Δ = coproduct(cube.structure)

    (e, i, (j, k)) = (edg.sign, edg.from, edg.to)
    xᵢ = x.label[i]
    
    res = map(Δ(xᵢ)) do (yⱼ, yₖ, r) 
        label = copy(x.label)
        deleteat!(label, i)
        insert!(label, j, yⱼ)
        insert!(label, k, yₖ)
        y = KhChainGenerator(v, label)
        (y, e * r)
    end
end