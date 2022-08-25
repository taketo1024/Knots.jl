using ..Links: Link, crossingNum, signedCrossingNums
using ..Homology
using ..Utils

struct KhComplex{R} <: AbstractComplex{R, KhChainGenerator}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}
    _generatorsCache::Dict{Int, Vector{KhChainGenerator}} # cache

    KhComplex(l::Link, cube::KhCube{R}, degShift::Tuple{Int, Int}) where {R} = begin
        gCache = Dict{Int, Vector{KhChainGenerator}}()
        new{R}(l, cube, degShift, gCache)
    end
end

KhComplex(str::KhAlgStructure{R}, l::Link; shift=true) where {R} = begin 
    cube = KhCube(str, l)
    degShift = if shift
        (n₊, n₋) = signedCrossingNums(l)
        (-n₋, n₊ - 2n₋)
    else
        (0, 0)
    end
    KhComplex(l, cube, degShift)
end

function Homology.hDegRange(C::KhComplex) :: UnitRange{Int}
    n = crossingNum(C.link)
    base = C.degShift[1] # <= 0
    base : base + n
end

function Homology.generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    get!(C._generatorsCache, k) do 
        _generators(C, k)
    end
end

function Homology.differentiate(C::KhComplex{R}, ::Int, x::KhChainGenerator) :: Vector{Tuple{KhChainGenerator, R}} where {R}
    differentiate(C.cube, x)
end

function _generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    base = C.degShift[1] # <= 0
    vs = vertices(C.cube, k - base)
    reduce(vs; init=KhChainGenerator[]) do res, v
        append!(res, v.generators)
    end
end