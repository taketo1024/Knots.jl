using .Utils
using SparseArrays
using SmithNormalForm: Smith, smith

struct KhHomologySummand{R}
    rank::Int
    torsions::Vector{R}
end

function asString(s::KhHomologySummand{R}) :: String where {R}
    iszero(s) && return "0"
    
    symbol = (R == Int) ? "Z" : "R" # TODO
    res = (s.rank > 0) ? ["$symbol$( s.rank > 1 ? Utils.superscript(s.rank) : "" )"] : []
    for t in s.torsions
        push!(res, "$symbol/$t")
    end
    join(res, "⊕")
end

Base.zero(::Type{KhHomologySummand{R}}) where {R} = KhHomologySummand{R}(0, R[])
Base.iszero(s::KhHomologySummand{R}) where {R} = (s.rank == 0 && isempty(s.torsions))
Base.show(io::IO, s::KhHomologySummand{R}) where {R} = print(io, asString(s))

# KhHomology

struct KhHomology{R} 
    link::Link
    complex::KhComplex{R}
    _SNFCache::Dict{Int, Smith}
end

KhHomology(str::KhAlgStructure{R}, l::Link; shift=true) where {R} = begin 
    C = KhComplex(str, l; shift=shift)
    sCache = Dict{Int, Smith}()
    KhHomology(l, C, sCache)
end

KhHomology(l::Link; shift=true) where {R} = begin 
    A = KhAlgStructure(Kh)
    KhHomology(A, l, shift=shift)
end

Base.getindex(H::KhHomology{R}, k::Int) where {R} = compute(H, k)

function compute(H::KhHomology{R}, k::Int) :: KhHomologySummand{R} where {R}
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
    Aₖ = matrix(H.complex, k)
    nₖ = size(Aₖ)[2]

    nₖ == 0 && return zero(KhHomologySummand{R})

    # TODO: use cache
    Aₖ₋₁ = matrix(H.complex, k - 1)
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