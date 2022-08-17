using AbstractAlgebra: RingElement, Ring
using SparseArrays
using .Utils

const SpMatrix = SparseMatrixCSC

struct KhComplex{R <: RingElement, RR <: Ring}
    link::Link
    cube::KhCube{R, RR}
    degShift::Tuple{Int, Int}
    _generatorsCache::Dict{Int, Vector{KhChainGenerator}}
    _differentialCache::Dict{Int, SpMatrix{R}}
end

KhComplex(str::KhAlgStructure{R}, l::Link; shift=true) where {R <: RingElement} = begin 
    cube = KhCube(str, l)
    degShift = if shift
        (n₊, n₋) = signedCrossingNums(l)
        (-n₋, n₊ - 2n₋)
    else
        (0, 0)
    end
    gCache = Dict{Int, Vector{KhChainGenerator}}()
    mCache = Dict{Int, SpMatrix{R}}()
    KhComplex(l, cube, degShift, gCache, mCache)
end

function hDegRange(C::KhComplex) :: UnitRange{Int}
    n = crossingNum(C.link)
    base = C.degShift[1] # <= 0
    base : base + n
end

function generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    get!(C._generatorsCache, k) do 
        base = C.degShift[1] # <= 0
        _generators(C.cube, k - base)
    end
end

function _generators(cube::KhCube, degree::Int) :: Vector{KhChainGenerator}
    n = dim(cube)

    if degree ∉ 0 : n
        []
    elseif n == degree == 0
        u = Int[]
        vertex(cube, u).generators
    else 
        bits = Utils.bitseq(n, degree)
        reduce(1 : length(bits); init=[]) do res, i
            u = digits(bits[i], base=2, pad=n)
            gens = vertex(cube, u).generators
            append!(res, gens)
        end
    end
end

function differential(C::KhComplex{R}, k::Int) :: SpMatrix{R} where {R <: RingElement}
    get!(C._differentialCache, k) do 
        base = C.degShift[1] # <= 0
        _differential(C.cube, k - base)
    end
end

function _differential(cube::KhCube{R}, degree::Int) :: SpMatrix{R} where {R <: RingElement}
    k = degree

    Gₖ   = _generators(cube, k)
    Gₖ₊₁ = _generators(cube, k + 1)

    n = length(Gₖ)
    m = length(Gₖ₊₁)

    gDict = _generatorsDict(Gₖ₊₁)
    vDict = Dict()
    
    Is = Vector{Int}()
    Js = Vector{Int}()
    Vs = Vector{R}()

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