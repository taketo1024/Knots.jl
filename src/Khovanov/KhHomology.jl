using SparseArrays
using ..Links: Link
using ..Homology
using ..Utils

# KhHomology

struct KhHomology{R} <: AbstractHomology{R} 
    complex::KhComplex{R}
end

function KhHomology(l::Link, str::KhAlgStructure{R}; reduced=false, shifted=true, perform_reduction=true, with_transform=false) where {R}
    KhHomology(
        KhComplex(l, str; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction, with_transform=with_transform)
    )
end

function KhHomology(l::Link, name::String; reduced=false, shifted=true, perform_reduction=true, with_transform=false)
    A = KhAlgStructure(name)
    KhHomology(l, A; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction, with_transform=with_transform)
end

function KhHomology(l::Link, h::R, t::R; reduced=false, shifted=true, perform_reduction=true, with_transform=false) where {R} 
    A = KhAlgStructure(h, t)
    KhHomology(l, A; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction, with_transform=with_transform)
end

function KhHomology(l::Link, ::Type{R}; reduced=false, shifted=true, perform_reduction=true, with_transform=false) where {R} 
    A = KhAlgStructure(zero(R), zero(R))
    KhHomology(l, A; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction, with_transform=with_transform)
end

function KhHomology(l::Link; reduced=false, shifted=true, perform_reduction=true, with_transform=false) where {R} 
    KhHomology(l, Int; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction, with_transform=with_transform)
end

Homology.complex(H::KhHomology) = 
    H.complex

function Homology.compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R}
    C = H.complex
    p = C.perform_reduction
    t = C.with_tranform

    Homology.compute_single(H, k; preprocess=!p, with_transform=t)
end

function Homology.asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    str = invoke(Homology.asString, Tuple{AbstractHomology}, H)
    lines = ["L = $l", A, "---", str, "---"]
    join(lines, "\n")
end

# KhHomologySummand

struct KhHomologySummand{R} <: AbstractHomologySummand{R}
    parent::KhHomology{R}
    h_degree::Int
    rank::Int
    torsions::Vector{R}
    transform::SparseMatrix{R}
end

Homology.h_degree(s::KhHomologySummand) = s.h_degree
Homology.parent(s::KhHomologySummand) = s.parent
Homology.rank(s::KhHomologySummand) = s.rank
Homology.torsions(s::KhHomologySummand) = s.torsions
Homology.transform(s::KhHomologySummand) = s.transform

Homology.makeSummand(H::KhHomology{R}, h_degree::Int, rank::Int, torsions::Vector{R}, transform::SparseMatrix{R}) where {R} = 
    KhHomologySummand(H, h_degree, rank, torsions, transform)

function Homology.vectorize(s::KhHomologySummand{R}, z::KhChain{R}) :: Vector{R} where {R}
    Homology.vectorize(s, z.elements) 
end