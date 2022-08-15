using AbstractAlgebra
using SparseArrays
using .Utils

struct KhHomologySummand{R <: RingElement}
    rank::Int
    torsions::Vector{R}
end

function asString(s::KhHomologySummand{R}, R_symbol="R") :: String where {R <: RingElement}
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
end

KhHomology(str::KhAlgStructure{R}, l::Link; shift=true) where {R <: RingElement} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, SNF{R}}()
    KhHomology(l, C, sCache)
end

function hDegRange(H::KhHomology{R}) :: UnitRange{Int} where {R <: RingElement}
    hDegRange(H.complex)
end

function compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R <: RingElement}
    #          Aₖ₋₁        Aₖ
    #     Cₖ₋₁ -----> Cₖ ------> Cₖ₊₁
    #     ^           | 
    #   Q |           | P 
    #     |    Dₖ₋₁   V
    #     Cₖ₋₁ -----> Cₖ' 
    #                 ⊕    Bₖ
    #                 Cₖ''-----> Cₖ₊₁
    #
    #   Hₖ = Ker(dₖ) / Im(dₖ₋₁)
    #      ≅ Ker(Dₖ) ⊕ Coker(Dₖ₋₁)
    #         ^ free    ^ tor

    str = H.complex.cube.structure
    ref = str.h

    Aₖ = differential(H.complex, k)
    nₖ = size(Aₖ)[2]
    nₖ == 0 && return zero(KhHomologySummand{R})

    # TODO: use cache
    Aₖ₋₁ = differential(H.complex, k - 1)
    Fₖ₋₁ = snf(Aₖ₋₁, ref)

    # non-zero diagonal entries of SNF(Dₖ₋₁)
    eₖ₋₁ = Fₖ₋₁.diag
    rₖ₋₁ = length(eₖ₋₁)

    P⁻¹ = Fₖ₋₁.P⁻¹
    Bₖ = Aₖ * view(P⁻¹, :, rₖ₋₁ + 1 : nₖ)
    Fₖ = snf(Bₖ, ref)
    rₖ = length(filter(r -> !iszero(r), Fₖ.diag))

    zₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), eₖ₋₁)

    KhHomologySummand(zₖ, tors)
end

function asString(H::KhHomology{R}) :: String where {R <: RingElement}
    L = H.complex.cube.link
    A = H.complex.cube.structure
    lines = ["L = $L", A, "---"]
    for i in hDegRange(H)
        push!(lines, "H[$i] = $(asString(H[i], A.R_symbol))")
    end
    push!(lines, "---")
    join(lines, "\n")
end

Base.getindex(H::KhHomology{R}, k::Int) where {R <: RingElement} = 
    compute(H, k)

Base.show(io::IO, H::KhHomology{R}) where {R <: RingElement} = 
    print(io, asString(H))