using ..Links: Link, State, Component, Edge, resolve, components, crossingNum, edges as l_edges, isEmpty
using ..Utils: findfirst_elm, delete_where!

# KhCubeVertex

struct KhCubeVertex
    state::State
    circles::Vector{Component}
    generators::Vector{KhChainGenerator}
end

function KhCubeVertex(l::Link, s::State, reduce_at::Union{Edge, Nothing}=nothing) :: KhCubeVertex
    circles = components(resolve(l, s))
    r = length(circles)
    generators = map(0 : 2^r - 1) do i in
        bits = digits(i, base=2, pad=r)
        label = map(b -> (b == 0) ? X : I, bits)
        KhChainGenerator(s, label)
    end

    if !isnothing(reduce_at)
        e = reduce_at
        i = findfirst(c -> e in c.edges, circles)
        delete_where!(generators) do x 
            x.label[i] == I
        end
    end

    KhCubeVertex(s, circles, generators)
end

# KhCubeEdge

struct KhCubeMergeEdge
    from::State
    to::State
    sign::Int
    transition::Tuple{Tuple{Int, Int}, Int}
end

struct KhCubeSplitEdge
    from::State
    to::State
    sign::Int
    transition::Tuple{Int, Tuple{Int, Int}}
end

const KhCubeEdge = Union{KhCubeMergeEdge, KhCubeSplitEdge}

function mergeEdge(from::State, to::State, sign::Int, transition::Tuple{Tuple{Int, Int}, Int}) :: KhCubeMergeEdge
    KhCubeMergeEdge(from, to, sign, transition)
end

function splitEdge(from::State, to::State, sign::Int, transition::Tuple{Int, Tuple{Int, Int}}) :: KhCubeSplitEdge
    KhCubeSplitEdge(from, to, sign, transition)
end

Base.:(==)(e1::KhCubeEdge, e2::KhCubeEdge) :: Bool = 
    (e1.from, e1.to, e1.sign, e1.transition) == (e2.from, e2.to, e2.sign, e2.transition)

# KhCube

struct KhCube{R} 
    structure::KhAlgStructure{R}
    link::Link
    reduced::Bool
    reduce_at::Union{Edge, Nothing}
    vertices::Dict{State, KhCubeVertex}    # cache
    edges::Dict{State, Vector{KhCubeEdge}} # cache

    KhCube(str::KhAlgStructure{R}, l::Link; reduced=false) where {R} = begin
        @assert !reduced || iszero(str.t)
        vertices = Dict{State, KhCubeVertex}()
        edges = Dict{State, Vector{KhCubeEdge}}()
        reduce_at = reduced && !isEmpty(l) ? minimum(l_edges(l)) : nothing
        new{R}(str, l, reduced, reduce_at, vertices, edges)
    end
end

function dim(cube::KhCube)
    crossingNum(cube.link)
end

function vertex(cube::KhCube, u::State) :: KhCubeVertex
    @assert length(u) == dim(cube)

    get!(cube.vertices, u) do 
        KhCubeVertex(cube.link, u, cube.reduce_at)
    end
end

function vertices(cube::KhCube, degree::Int) :: Vector{KhCubeVertex}
    n = dim(cube)

    return if degree ∉ 0 : n
        []
    elseif n == degree == 0
        u = Int[]
        [vertex(cube, u)]
    else 
        bits = Utils.bitseq(n, degree)
        map(bits) do b
            u = digits(b, base=2, pad=n)
            vertex(cube, u)
        end
    end
end

function nextVertices(cube::KhCube, u::State) :: Vector{State}
    @assert length(u) == dim(cube)

    n = length(u)
    indices = filter(i -> u[i] == 0, 1 : n)
    map(indices) do i in
        v = copy(u)
        v[i] = 1
        v
    end
end

function edgeSign(cube::KhCube, u::State, v::State) :: Int
    @assert length(u) == length(v) == dim(cube)
    @assert sum(u) + 1 == sum(v)

    n = dim(cube)
    i = findfirst(i -> u[i] != v[i], 1 : n)
    k = count(1 : i) do i
        u[i] == v[i] == 1
    end
    (-1)^isodd(k)
end

function edge(cube::KhCube, u::State, v::State) :: KhCubeEdge
    @assert length(u) == length(v) == dim(cube)
    @assert sum(u) + 1 == sum(v)

    e = edgeSign(cube, u, v)
    Cᵤ = vertex(cube, u).circles
    Cᵥ = vertex(cube, v).circles
    cᵤ = filter(c -> c ∉ Cᵥ, Cᵤ)
    cᵥ = filter(c -> c ∉ Cᵤ, Cᵥ)

    return if (length(cᵤ), length(cᵥ)) == (2, 1)
        i₁ = findfirst_elm(Cᵤ, cᵤ[1])
        i₂ = findfirst_elm(Cᵤ, cᵤ[2])
        j  = findfirst_elm(Cᵥ, cᵥ[1])
        if i₁ > i₂ 
            (i₁, i₂) = (i₂, i₁)
        end
        KhCubeMergeEdge(u, v, e, ((i₁, i₂), j))

    elseif (length(cᵤ), length(cᵥ)) == (1, 2)
        i  = findfirst_elm(Cᵤ, cᵤ[1])
        j₁ = findfirst_elm(Cᵥ, cᵥ[1])
        j₂ = findfirst_elm(Cᵥ, cᵥ[2])
        if j₁ > j₂ 
            (j₁, j₂) = (j₂, j₁)
        end
        KhCubeSplitEdge(u, v, e, (i, (j₁, j₂)))

    else
        throw(Exception)
    end
end

function edges(cube::KhCube, u::State) :: Vector{KhCubeEdge}
    get!(cube.edges, u) do 
        vs = nextVertices(cube, u)
        map(vs) do v
            edge(cube, u, v)
        end
    end
end

function apply(cube::KhCube{R}, edg::KhCubeMergeEdge, x::KhChainGenerator) :: Vector{Pair{KhChainGenerator, R}} where {R}
    m = product(cube.structure) 

    (v, e, ((i, j), k)) = (edg.to, edg.sign, edg.transition)
    (xᵢ, xⱼ) = (x.label[i], x.label[j])

    map(m(xᵢ, xⱼ)) do (yₖ, r) 
        label = copy(x.label)
        deleteat!(label, j)
        deleteat!(label, i)
        insert!(label, k, yₖ)
        y = KhChainGenerator(v, label)
        (y => e * r)
    end
end

function apply(cube::KhCube{R}, edg::KhCubeSplitEdge, x::KhChainGenerator) :: Vector{Pair{KhChainGenerator, R}} where {R}
    Δ = coproduct(cube.structure)

    (v, e, (i, (j, k))) = (edg.to, edg.sign, edg.transition)
    xᵢ = x.label[i]
    
    map(Δ(xᵢ)) do (yⱼ, yₖ, r) 
        label = copy(x.label)
        deleteat!(label, i)
        insert!(label, j, yⱼ)
        insert!(label, k, yₖ)
        y = KhChainGenerator(v, label)
        (y => e * r)
    end
end

# KhComplex IF

function chain_generators(cube::KhCube, degree::Int) :: Vector{KhChainGenerator}
    Vs = vertices(cube, degree)
    reduce(Vs; init=KhChainGenerator[]) do res, V
        append!(res, V.generators)
    end
end

function differentiate(cube::KhCube{R}, x::KhChainGenerator) :: Vector{Pair{KhChainGenerator, R}} where {R}
    u = x.state
    es = edges(cube, u)

    reduce(es; init=[]) do res, e
        y = apply(cube, e, x)
        append!(res, y)
    end
end