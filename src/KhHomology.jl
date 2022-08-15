using AbstractAlgebra
using SparseArrays
using .Utils

struct KhHomologySummand{R <: RingElement}
    rank::Int
    torsions::Vector{R}
end

function asString(s::KhHomologySummand{R}; R_symbol="R") :: String where {R <: RingElement}
    iszero(s) && return "⋅"
    res = (s.rank > 0) ? ["$R_symbol$( s.rank > 1 ? Utils.superscript(s.rank) : "" )"] : []
    for t in s.torsions
        push!(res, "$R_symbol/$t")
    end
    join(res, " ⊕ ")
end

Base.zero(::Type{KhHomologySummand{R}}) where {R <: RingElement} = 
    KhHomologySummand{R}(0, R[])

Base.iszero(s::KhHomologySummand{R}) where {R <: RingElement} = 
    (s.rank == 0 && isempty(s.torsions))

Base.show(io::IO, s::KhHomologySummand{R}) where {R <: RingElement} = 
    print(io, asString(s))

# KhHomology

struct KhHomology{R <: RingElement} 
    link::Link
    complex::KhComplex{R}
    _SNFCache::Dict{Int, SNF{R}}
    _withGenerators::Bool
end

KhHomology(str::KhAlgStructure{R}, l::Link; shift=true, withGenerators=false) where {R <: RingElement} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, SNF{R}}()
    KhHomology(l, C, sCache, withGenerators)
end

function hDegRange(H::KhHomology{R}) :: UnitRange{Int} where {R <: RingElement}
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

    nₖ = length(chainGenerators(H.complex, k))
    nₖ == 0 && return zero(KhHomologySummand{R})

    Fₖ₋₁ = _snf(H, k - 1)
    Fₖ   = _snf(H, k)

    rₖ₋₁ = length(Fₖ₋₁.S)
    rₖ   = length(Fₖ.S)
    
    fₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), Fₖ₋₁.S)

    KhHomologySummand(fₖ, tors)
end

function _snf(H::KhHomology{R}, k::Int) :: SNF{R} where {R <: RingElement}
    get!(H._SNFCache, k) do 
        Aₖ = differential(H.complex, k)
        baseRing = H.complex.cube.structure.baseRing
        withTrans = H._withGenerators
        snf(Aₖ, baseRing; withTrans)
    end
end

function asString(H::KhHomology{R}) :: String where {R <: RingElement}
    L = H.complex.cube.link
    A = H.complex.cube.structure
    lines = ["L = $L", A, "---"]
    for i in hDegRange(H)
        push!(lines, "H[$i] = $(asString(H[i]; A.R_symbol))")
    end
    push!(lines, "---")
    join(lines, "\n")
end

Base.getindex(H::KhHomology{R}, k::Int) where {R <: RingElement} = 
    compute(H, k)

Base.show(io::IO, H::KhHomology{R}) where {R <: RingElement} = 
    print(io, asString(H))