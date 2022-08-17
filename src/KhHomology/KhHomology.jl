using AbstractAlgebra
using SparseArrays
using .Utils

struct KhHomologySummand{R<:RingElement, RR<:AbstractAlgebra.Ring} <: AbstractHomologySummand{R, RR}
    baseRing::RR
    rank::Int
    torsions::Vector{R}
end

baseRing(s::KhHomologySummand) = s.baseRing
rank(s::KhHomologySummand) = s.rank
torsions(s::KhHomologySummand) = s.torsions

# KhHomology

struct KhHomology{R <: RingElement} 
    link::Link
    complex::KhComplex{R}
    _SNFCache::Dict{Int, SNF{R}}
end

KhHomology(str::KhAlgStructure{R}, l::Link; shift=true) where {R <: RingElement} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, SNF{R}}()
    KhHomology(l, C, sCache)
end

function hDegRange(H::KhHomology) :: UnitRange{Int}
    hDegRange(H.complex)
end

function compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R <: RingElement}
    
    #           dₖ₋₁          dₖ
    #     Cₖ₋₁ -------> Cₖ --------> Cₖ₊₁
    #     :             :            :
    #     :     Aₖ₋₁    :      Aₖ     :
    #     Rⁿ ---------> Rᵐ --------> Rᵖ
    #     |             ^            | Pₖ
    #     |          Qₖ |      Sₖ     V
    #     |             Rʳ --------> Rʳ ⊕ ..
    #     |             ⊕ 
    #     |             Rᶠ \ 
    #     |     Sₖ₋₁    ⊕   | Zₖ
    #     Rᵗ ---------> Rᵗ /
    #     ⊕
    #     :
    #
    #    Hₖ = Ker(dₖ) / Im(dₖ₋₁)
    #       ≅ Rᶠ ⊕ Rᵗ/Im(Sₖ₋₁)

    baseRing = H.complex.cube.structure.baseRing

    nₖ = length(generators(H.complex, k))
    nₖ == 0 && return KhHomologySummand(baseRing, 0, R[])

    Fₖ₋₁ = _snf(H, k - 1)
    Fₖ   = _snf(H, k)

    rₖ₋₁ = length(Fₖ₋₁.S)
    rₖ   = length(Fₖ.S)
    
    fₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), Fₖ₋₁.S)

    KhHomologySummand(baseRing, fₖ, tors)
end

function _snf(H::KhHomology{R}, k::Int) :: SNF{R} where {R <: RingElement}
    get!(H._SNFCache, k) do 
        Aₖ = differential(H.complex, k)
        baseRing = H.complex.cube.structure.baseRing
        snf(Aₖ, baseRing)
    end
end

function asString(H::KhHomology) :: String
    l = H.complex.cube.link
    A = H.complex.cube.structure
    lines = ["L = $l", A, "---"]
    for i in hDegRange(H)
        Hi = asString(H[i])
        push!(lines, "H[$i] = $Hi")
    end
    push!(lines, "---")
    join(lines, "\n")
end

Base.getindex(H::KhHomology, k::Int) = 
    compute(H, k)

Base.show(io::IO, H::KhHomology) = 
    print(io, asString(H))