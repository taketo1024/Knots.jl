using AbstractAlgebra
using SparseArrays
using .Utils

Matrix = SparseMatrixCSC

struct KhComplex{R <: RingElement}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}
    _generatorsCache::Dict{Int, Vector{KhChainGenerator}}
    _matrixCache::Dict{Int, Matrix{R}}
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
    mCache = Dict{Int, Matrix{R}}()
    KhComplex(l, cube, degShift, gCache, mCache)
end

function hDegRange(C::KhComplex{R}) :: UnitRange{Int} where {R <: RingElement}
    n = crossingNum(C.link)
    base = C.degShift[1] # <= 0
    base : base + n
end

function chainGenerators(C::KhComplex{R}, k::Int) :: Vector{KhChainGenerator} where {R <: RingElement} 
    get!(C._generatorsCache, k) do 
        base = C.degShift[1] # <= 0
        gens = _chainGenerators(C.cube, k - base)
        reduce(gens; init=[]) do res, entry
            append!(res, entry[2])
        end
    end
end

function matrix(C::KhComplex{R}, k::Int) :: Matrix{R} where {R <: RingElement} 
    get!(C._matrixCache, k) do 
        base = C.degShift[1] # <= 0
        _matrix(C.cube, k - base)
    end
end

function _chainGenerators(cube::KhCube{R}, degree::Int) :: Vector{ Tuple{ State, Vector{KhChainGenerator} } } where {R <: RingElement} 
    n = dim(cube)

    if degree ∉ 0 : n
        []
    elseif n == degree == 0
        u = Int[]
        [(u, vertex(cube, u).generators)]
    else 
        vertices = Utils.bitseq(n, degree)
        reduce(1 : length(vertices); init=[]) do res, i
            u = digits(vertices[i], base=2, pad=n)
            gens = vertex(cube, u).generators
            push!(res, (u, gens))
        end
    end
end

function _indexDict(gens::Vector) :: Dict{KhChainGenerator, Int} 
    res = Dict()
    for (_, g) in gens
        N = length(res)
        merge!(res, Dict( g[i] => N + i for i in 1:length(g)) )
    end
    res
end

function _matrix(cube::KhCube{R}, degree::Int) :: Matrix{R} where {R <: RingElement}
    k = degree
    Gₖ   = _chainGenerators(cube, k)
    Gₖ₊₁ = _chainGenerators(cube, k + 1)
    
    Is = Vector{Int}()
    Js = Vector{Int}()
    Vs = Vector{R}()

    m = reduce( (res, e) -> res + length(e[2]), Gₖ₊₁, init=0)
    n = reduce( (res, e) -> res + length(e[2]), Gₖ,   init=0)

    dict = _indexDict(Gₖ₊₁)
    j = 1
    
    for (u, xs) in Gₖ 
        vs = nextVertices(cube, u)
        for x in xs
            for v in vs
                ys = edgeMap(cube, u, v, x)
                for (y, r) in ys
                    i = dict[y]
                    push!(Is, i)
                    push!(Js, j)
                    push!(Vs, r)
                end
            end
            j += 1
        end
    end
    
    sparse(Is, Js, Vs, m, n)
end