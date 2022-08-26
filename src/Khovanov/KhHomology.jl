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
    chain_reduction::Bool
end

KhHomology(str::KhAlgStructure{R}, l::Link; chain_reduction=true, shift=true) where {R} = begin 
    C = KhComplex(str, l; chain_reduction=chain_reduction, shift=shift)
    KhHomology{R}(l, C, chain_reduction)
end

Homology.complex(H::KhHomology) = 
    H.complex

Homology.makeSummand(H::KhHomology, rank::Int, torsions::Vector) = 
    KhHomologySummand(rank, torsions)

function Homology.compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R}
    prep = !H.chain_reduction
    Homology.compute_single(H, k; preprocess=prep)
end

function Homology.asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    str = invoke(Homology.asString, Tuple{AbstractHomology}, H)
    lines = ["L = $l", A, "---", str, "---"]
    join(lines, "\n")
end