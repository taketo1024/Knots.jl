using SparseArrays
using ..Links: Link, crossingNum, signedCrossingNums
using ..Homology
using ..Utils

const KhComplexMatrix = SparseMatrixCSC

struct KhComplex{R} <: AbstractComplex{R}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}
    _generatorsCache::Dict{Int, Vector{KhChainGenerator}}
    _differentialCache::Dict{Int, KhComplexMatrix{R}}

    KhComplex(l::Link, cube::KhCube{R}, degShift::Tuple{Int, Int}) where {R} = begin
        gCache = Dict{Int, Vector{KhChainGenerator}}()
        mCache = Dict{Int, KhComplexMatrix{R}}()
        new{R}(l, cube, degShift, gCache, mCache)
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

function Homology.differential(C::KhComplex{R}, k::Int) :: KhComplexMatrix{R} where {R}
    get!(C._differentialCache, k) do 
        _differential(C, k)
    end
end

function _generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    base = C.degShift[1] # <= 0
    vs = vertices(C.cube, k - base)
    reduce(vs; init=KhChainGenerator[]) do res, v
        append!(res, v.generators)
    end
end

function _differential(C::KhComplex{R}, k::Int) :: KhComplexMatrix{R} where {R}
    Gₖ   = generators(C, k)
    Gₖ₊₁ = generators(C, k + 1)

    n = length(Gₖ)
    m = length(Gₖ₊₁)

    gDict = _generatorsDict(Gₖ₊₁)
    vDict = Dict()
    
    Is = Vector{Int}()
    Js = Vector{Int}()
    Vs = Vector{R}()

    cube = C.cube

    for j in 1 : n 
        x = Gₖ[j]
        u = x.state
        vs = get!(vDict, u) do 
            nextVertices(cube, u)
        end

        for v in vs
            ys = edgeMap(cube, u, v, x)
            for (y, r) in ys
                i = gDict[y]
                push!(Is, i)
                push!(Js, j)
                push!(Vs, r)
            end
        end
    end
    
    sparse(Is, Js, Vs, m, n)
end

function _generatorsDict(gens::Vector{KhChainGenerator}) :: Dict{KhChainGenerator, Int} 
    Dict( gens[i] => i for i in 1 : length(gens) )
end