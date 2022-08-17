using AbstractAlgebra: RingElement, Ring
using SparseArrays
using .Utils

struct KhHomologySummand{R<:RingElement, RR<:Ring} <: AbstractHomologySummand{R, RR}
    baseRing::RR
    rank::Int
    torsions::Vector{R}
end

baseRing(s::KhHomologySummand) = s.baseRing
rank(s::KhHomologySummand) = s.rank
torsions(s::KhHomologySummand) = s.torsions

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

function complex(H::KhHomology{R, RR}) :: KhComplex{R, RR} where {R, RR}
    H.complex
end

function makeSummand(H::KhHomology{R, RR}, rank::Int, torsions::Vector{R}) :: KhHomologySummand{R, RR} where {R, RR}
    KhHomologySummand(baseRing(H), rank, torsions)
end

function compute(H::KhHomology{R, RR}, k::Int) :: KhHomologySummand{R, RR} where {R, RR}
    compute_single(H, k)
end

function asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    str = invoke(asString, Tuple{AbstractHomology}, H)
    lines = ["L = $l", A, "---", str]
    join(lines, "\n")
end