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

KhComplex(str::KhAlgStructure{R}, l::Link; chain_reduction=true, shift=true) where {R} = begin 
    cube = KhCube(str, l)
    chain_reduction && reduce!(cube)

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
        base = C.degShift[1] # <= 0
        generators(C.cube, k - base)
    end
end

function Homology.differentiate(C::KhComplex{R}, ::Int, x::KhChainGenerator) :: Vector{Pair{KhChainGenerator, R}} where {R}
    differentiate(C.cube, x)
end