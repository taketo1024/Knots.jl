struct KhComplex{R}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}
end

KhComplex(str::KhAlgStructure, l::Link; shift=true) where {R} = begin 
    cube = KhCube(str, l)
    degShift = if shift
        (n₊, n₋) = signedCrossingNums(l)
        (-n₋, n₊ - 2n₋)
    else
        (0, 0)
    end
    KhComplex(l, cube, degShift)
end

function hDegRange(C::KhComplex{R}) :: UnitRange{Int} where {R}
    n = crossingNum(C.link)
    base = C.degShift[1] # <= 0
    base : base + n
end

function _bitseq(length::Int, degree::Int) :: Vector{Int} 
    if length <= 0 || degree < 0 || length < degree 
        []
    elseif length > 1
        s₀ = map( b -> (b << 1) | 1, _bitseq(length - 1, degree - 1) )
        s₁ = map( b -> b << 1, _bitseq(length - 1, degree) )
        append!(s₀, s₁)
    else # 0 ≤ degree ≤ length == 1
        [degree] 
    end
end

function _collectGenerators(cube::KhCube{R}, degree::Int) :: Vector{ Tuple{ State, Vector{KhChainGenerator} } } where {R} 
    n = dim(cube)
    if degree ∉ 0 : n
        return []
    end

    vertices = _bitseq(n, degree)
    reduce(1 : length(vertices); init=[]) do res, i
        u = digits(vertices[i], base=2, pad=n)
        gens = vertex(cube, u).generators
        push!(res, (u, gens))
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

using SparseArrays

function _matrix(cube::KhCube{R}, degree::Int) :: SparseMatrixCSC{R} where {R}
    k = degree
    Gₖ   = _collectGenerators(cube, k)
    Gₖ₊₁ = _collectGenerators(cube, k + 1)
    
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