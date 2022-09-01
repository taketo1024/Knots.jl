using Permutations
using ..Links: Link, crossingNum, signedCrossingNums
using ..Matrix
using ..Homology
using ..Utils

struct KhComplex{R} <: AbstractComplex{R, KhChainGenerator}
    link::Link
    cube::KhCube{R}
    degShift::Tuple{Int, Int}

    perform_reduction::Bool
    with_tranform::Bool

    generators::Dict{Int, Vector{KhChainGenerator}} # cache
    differentials::Dict{Int, SparseMatrix{R}} # cache
    transforms::Dict{Int, SparseMatrix{R}} # cache

    KhComplex(l::Link, cube::KhCube{R}, degShift::Tuple{Int, Int}; perform_reduction=true, with_transform=false) where {R} = begin
        generators = Dict{Int, Vector{KhChainGenerator}}()
        differentials = Dict{Int, SparseMatrix{R}}()
        transforms = Dict{Int, SparseMatrix{R}}()

        new{R}(l, cube, degShift, perform_reduction, with_transform, generators, differentials, transforms)
    end
end

function KhComplex(str::KhAlgStructure{R}, l::Link; reduced=false, shifted=true, perform_reduction=true, with_transform=false) where {R}
    cube = KhCube(str, l; reduced=reduced)
    degShift = if shifted
        (n₊, n₋) = signedCrossingNums(l)
        (-n₋, n₊ - 2n₋ + (reduced ? 1 : 0))
    else
        (0, 0)
    end

    C = KhComplex(l, cube, degShift; perform_reduction=perform_reduction, with_transform=with_transform)

    if perform_reduction
        reduce_all!(C; with_transform=with_transform)
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
        Homology.original_generators(C, k)
    end
end

function Homology.original_generators(C::KhComplex, k::Int) :: Vector{KhChainGenerator}
    base = C.degShift[1] # <= 0
    chain_generators(C.cube, k - base)
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
        Homology.generate_differential(C, k)
    end
end

function Homology.set_differential!(C::KhComplex{R}, k::Int, A::SparseMatrix{R}) where {R}
    C.differentials[k] = A
end

function Homology.transform(C::KhComplex{R}, k::Int) :: Union{SparseMatrix{R}, Nothing} where {R}
    get(C.transforms, k, nothing)
end

function Homology.set_transform!(C::KhComplex{R}, k::Int, T::SparseMatrix{R}) where {R}
    if haskey(C.transforms, k)
        T = T * C.transforms[k]
    end
    C.transforms[k] = T
end

function Homology.vectorize(C::KhComplex{R}, k::Int, z::KhChain{R}) :: Vector{R} where {R}
    Homology.vectorize(C, k, z.elements)
end