using AbstractAlgebra
using SparseArrays
using SmithNormalForm: Smith, smith
using .Utils

struct KhHomologySummand{R <: RingElement}
    rank::Int
    torsions::Vector{R}
end

function asString(s::KhHomologySummand{R}) :: String where {R <: RingElement}
    iszero(s) && return "0"
    
    symbol = (R == Int) ? "Z" : "R" # TODO
    res = (s.rank > 0) ? ["$symbol$( s.rank > 1 ? Utils.superscript(s.rank) : "" )"] : []
    for t in s.torsions
        push!(res, "$symbol/$t")
    end
    join(res, "⊕")
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
    _SNFCache::Dict{Int, Smith}
end

KhHomology(str::KhAlgStructure{R}, l::Link; shift=true) where {R <: RingElement} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, Smith}()
    KhHomology(l, C, sCache)
end

Base.getindex(H::KhHomology{R}, k::Int) where {R <: RingElement} = 
    compute(H, k)

function compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R <: RingElement}
    #          Aₖ₋₁        Aₖ
    #     Cₖ₋₁ -----> Cₖ ------> Cₖ₊₁
    #     |           ^ 
    #   T |           | S 
    #     V    Dₖ₋₁   |
    #     Cₖ₋₁ -----> Cₖ' 
    #                 ⊕    Dₖ
    #                 Cₖ''-----> Cₖ₊₁
    #
    #   Hₖ = Ker(dₖ) / Im(dₖ₋₁)
    #      ≅ Ker(Dₖ) ⊕ Coker(Dₖ₋₁)
    #         ^ free    ^ tor

    str = H.complex.cube.structure
    Aₖ = differential(H.complex, k)
    nₖ = size(Aₖ)[2]

    nₖ == 0 && return zero(KhHomologySummand{R})

    # TODO: use cache
    Aₖ₋₁ = differential(H.complex, k - 1)
    Fₖ₋₁ = smith(Aₖ₋₁)

    # non-zero diagonal entries of SNF(Dₖ₋₁)
    eₖ₋₁ = Vector(filter(r -> !iszero(r), Fₖ₋₁.SNF)) 
    rₖ₋₁ = length(eₖ₋₁)

    # Cₖ'  = S[:, 1 : rₖ₋₁],      rk = rₖ₋₁,
    # Cₖ'' = S[:, rₖ₋₁ + 1 : nₖ], rk = nₖ - rₖ₋₁.

    S = Fₖ₋₁.S
    Dₖ = Aₖ * view(S, :, rₖ₋₁ + 1 : nₖ)
    Fₖ = smith(Dₖ)

    rₖ = length(filter(r -> !iszero(r), Fₖ.SNF))
    zₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> r ∉ (str.one, -str.one), eₖ₋₁)

    KhHomologySummand(zₖ, tors)
end