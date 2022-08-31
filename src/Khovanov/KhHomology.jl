using SparseArrays
using ..Links: Link
using ..Homology
using ..Utils

struct KhHomologySummand{R} <: AbstractHomologySummand{R}
    rank::Int
    torsions::Vector{R}
end

Homology.rank(s::KhHomologySummand) = s.rank
Homology.torsions(s::KhHomologySummand) = s.torsions

# KhHomology

struct KhHomology{R} <: AbstractHomology{R} 
    link::Link
    complex::KhComplex{R}
end

KhHomology(str::KhAlgStructure{R}, l::Link; reduced=false, shifted=true, perform_reduction=true) where {R} = begin 
    C = KhComplex(str, l; reduced=reduced, shifted=shifted, perform_reduction=perform_reduction)
    KhHomology{R}(l, C)
end

Homology.complex(H::KhHomology) = 
    H.complex

Homology.makeSummand(H::KhHomology, rank::Int, torsions::Vector) = 
    KhHomologySummand(rank, torsions)

function Homology.compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R}
    Homology.compute_single(H, k; preprocess=false)
end

function Homology.asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    str = invoke(Homology.asString, Tuple{AbstractHomology}, H)
    lines = ["L = $l", A, "---", str, "---"]
    join(lines, "\n")
end