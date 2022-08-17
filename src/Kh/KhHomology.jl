using AbstractAlgebra: RingElement, Ring
using SparseArrays
using ..Links: Link
using ..Matrix: SNF
using ..Homology
using ..Utils

struct KhHomologySummand{R<:RingElement, RR<:Ring} <: AbstractHomologySummand{R, RR}
    baseRing::RR
    rank::Int
    torsions::Vector{R}
end

Homology.baseRing(s::KhHomologySummand) = s.baseRing
Homology.rank(s::KhHomologySummand) = s.rank
Homology.torsions(s::KhHomologySummand) = s.torsions

# KhHomology

struct KhHomology{R<:RingElement, RR<:Ring} <: AbstractHomology{R, RR} 
    link::Link
    complex::KhComplex{R, RR}
    _SNFCache::Dict{Int, SNF{R}}
end

KhHomology(str::KhAlgStructure{R, RR}, l::Link; shift=true) where {R, RR} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, SNF{R}}()
    KhHomology{R, RR}(l, C, sCache)
end

Homology.complex(H::KhHomology) = 
    H.complex

Homology.makeSummand(H::KhHomology, rank::Int, torsions::Vector) = 
    KhHomologySummand(baseRing(H), rank, torsions)

function Homology.compute(H::KhHomology{R, RR}, k::Int) :: KhHomologySummand{R, RR} where {R, RR}
    Homology.compute_single(H, k)
end

function Homology.asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    str = invoke(Homology.asString, Tuple{AbstractHomology}, H)
    lines = ["L = $l", A, "---", str]
    join(lines, "\n")
end