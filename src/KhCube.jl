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

# KhCube

mutable struct KhCube{R} 
    structure::KhAlgStructure{R}
    link::Link
    vertices::Dict{State, KhCubeVertex}
end

KhCube(str::KhAlgStructure{R}, l::Link) where {R} = KhCube(str, l, Dict{State, KhCubeVertex}())

function dim(cube::KhCube{R}) where {R} 
    crossingNum(cube.link)
end

function vertex(cube::KhCube{R}, u::State) :: KhCubeVertex where {R}
    get!(cube.vertices, u, KhCubeVertex(cube.link, u))
end

abstract type EdgeType end

struct mergeEdge <: EdgeType
    from::Tuple{Int, Int}
    to::Int
end

struct splitEdge <: EdgeType
    from::Int
    to::Tuple{Int, Int}
end

Base.findfirst(arr::AbstractArray{T, N}, elm::T) where {T, N} = findfirst( x -> x == elm, arr )

function edge(cube::KhCube{R}, u::State, v::State) :: Union{mergeEdge, splitEdge} where {R}
    @assert sum(u) + 1 == sum(v)
    Cᵤ = vertex(cube, u).circles
    Cᵥ = vertex(cube, v).circles
    cᵤ = filter(c -> c ∉ Cᵥ, Cᵤ)
    cᵥ = filter(c -> c ∉ Cᵤ, Cᵥ)

    return if (length(cᵤ), length(cᵥ)) == (2, 1)
        println(Cᵤ, Cᵥ, cᵤ, cᵥ)
        i₁ = findfirst(Cᵤ, cᵤ[1])
        i₂ = findfirst(Cᵤ, cᵤ[2])
        j  = findfirst(Cᵥ, cᵥ[1])
        mergeEdge((i₁, i₂), j)

    elseif (length(cᵤ), length(cᵥ)) == (1, 2)
        i  = findfirst(Cᵤ, cᵤ[1])
        j₁ = findfirst(Cᵥ, cᵥ[1])
        j₂ = findfirst(Cᵥ, cᵥ[2])
        splitEdge(i, (j₁, j₂))

    else
        throw(Exception)
    end
end
