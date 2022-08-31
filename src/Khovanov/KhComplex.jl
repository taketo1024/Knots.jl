using Permutations
using ..Links: Link, crossingNum, signedCrossingNums
using ..Matrix
using ..Homology
using ..Utils

struct KhComplex{R} <: AbstractComplex{R, KhChainGenerator}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}

    generators::Dict{Int, Vector{KhChainGenerator}} # cache
    differentials::Dict{Int, SparseMatrix{R}} # cache

    KhComplex(l::Link, cube::KhCube{R}, degShift::Tuple{Int, Int}) where {R} = begin
        generators = Dict{Int, Vector{KhChainGenerator}}()
        differentials = Dict{Int, SparseMatrix{R}}()
        new{R}(l, cube, degShift, generators, differentials)
    end
end

function KhComplex(str::KhAlgStructure{R}, l::Link; reduced=false, shifted=true, perform_reduction=true) where {R}
    cube = KhCube(str, l; reduced=reduced)
    degShift = if shifted
        (n₊, n₋) = signedCrossingNums(l)
        (-n₋, n₊ - 2n₋ + (reduced ? 1 : 0))
    else
        (0, 0)
    end

    C = KhComplex(l, cube, degShift)

    if perform_reduction
        reduce_all!(C)
    end

    C
end

function Homology.hDegRange(C::KhComplex) :: UnitRange{Int}
    n = crossingNum(C.link)
    base = C.degShift[1] # <= 0
    base : base + n
end

function Homology.generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    get!(C.generators, k) do 
        base = C.degShift[1] # <= 0
        chain_generators(C.cube, k - base)
    end
end

function Homology.drop_generators!(C::KhComplex, k::Int, r::Int, p::Permutation)
    gens = generators(C, k)
    n = length(gens)

    C.generators[k] = [ gens[p(i)] for i in r+1:n ]

    k₋₁ = k - differentialDegree(C)
    if k₋₁ in keys(C.differentials)
        A = C.differentials[k₋₁]
        C.differentials[k₋₁] = permute_row(A, p)[r+1:n, :]
    end
end

function Homology.differentiate(C::KhComplex{R}, ::Int, x::KhChainGenerator) :: Vector{Pair{KhChainGenerator, R}} where {R}
    differentiate(C.cube, x)
end

function Homology.differential(C::KhComplex{R}, k::Int) :: SparseMatrix{R} where {R}
    get!(C.differentials, k) do 
        invoke(Homology.differential, Tuple{AbstractComplex{R, KhChainGenerator}, Int}, C, k)
    end
end

function Homology.set_differential!(C::KhComplex{R}, k::Int, A::SparseMatrix{R}) where {R}
    C.differentials[k] = A
end