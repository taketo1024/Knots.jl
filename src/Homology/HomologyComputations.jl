using AbstractAlgebra: RingElement, Ring, is_unit
using ..Matrix: SNF, snf

#  Homology computations:
#
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

function compute_single(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    C = complex(H)
    deg = differentialDegree(C)

    nₖ = length(generators(C, k))
    nₖ == 0 && return makeSummand(H, 0, R[])

    Aₖ₋₁ = differential(C, k - deg)
    Aₖ   = differential(C, k)

    Fₖ₋₁ = snf(Aₖ₋₁)
    Fₖ   = snf(Aₖ)

    rₖ₋₁ = length(Fₖ₋₁.factors)
    rₖ   = length(Fₖ.factors)
    
    fₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), Fₖ₋₁.factors)

    makeSummand(H, fₖ, tors)
end

function compute_incremental(H::AbstractHomology{R, RR}, k::Int; previous=nothing) :: Tuple{AbstractHomologySummand{R, RR}, Any} where {R, RR}
    C = complex(H)
    deg = differentialDegree(C)
    flags = (false, true, false, false) # only need P⁻¹

    nₖ = length(generators(C, k))
    nₖ == 0 && return (makeSummand(H, 0, R[]), nothing)

    nₖ₋₁ = length(generators(C, k - deg))
    if nₖ₋₁ > 0
        Aₖ₋₁ = differential(C, k - deg)
        Fₖ₋₁ = isnothing(previous) ? 
            snf(Aₖ₋₁; flags=flags) :
            previous

        Pₖ₋₁⁻¹ = Fₖ₋₁.T.P⁻¹
        eₖ₋₁ = Fₖ₋₁.factors
        rₖ₋₁ = length(eₖ₋₁)
        tors = filter(r -> !is_unit(r), eₖ₋₁)
    else
        Pₖ₋₁⁻¹ = nothing # not to be acccessed
        rₖ₋₁ = 0
        tors = R[]
    end

    Aₖ = differential(C, k)
    Bₖ = (rₖ₋₁ > 0) ? 
        Aₖ * Pₖ₋₁⁻¹[:, rₖ₋₁ + 1 : nₖ] : # restrict Aₖ to the complement of R^{rₖ₋₁}
        Aₖ

    Fₖ = snf(Bₖ; flags=flags)
    eₖ = Fₖ.factors
    rₖ = length(eₖ)

    fₖ = nₖ - rₖ₋₁ - rₖ

    (makeSummand(H, fₖ, tors), Fₖ)
end

function compute_reverse_incremental(H::AbstractHomology{R, RR}, k::Int; previous=nothing) :: Tuple{AbstractHomologySummand{R, RR}, Any} where {R, RR}
    C = complex(H)
    deg = differentialDegree(C)
    flags=(false, false, false, true) # only need Q⁻¹

    nₖ = length(generators(C, k))
    nₖ == 0 && return (makeSummand(H, 0, R[]), nothing)

    Fₖ = if isnothing(previous)
        Aₖ = differential(H.complex, k)
        snf(Aₖ; flags=flags) 
    else 
        previous
    end

    rₖ = length(Fₖ.factors)
    Qₖ⁻¹ = Fₖ.T.Q⁻¹

    Aₖ₋₁ = differential(H.complex, k - deg)
    Bₖ₋₁ = (rₖ < nₖ) ? 
        Qₖ⁻¹[rₖ + 1 : nₖ, : ] * Aₖ₋₁ : # compose with projection to Zₖ
        Aₖ₋₁
    
    Fₖ₋₁ = snf(Bₖ₋₁; flags=flags)
    eₖ₋₁ = Fₖ₋₁.factors
    rₖ₋₁ = length(eₖ₋₁)

    fₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), eₖ₋₁)

    (makeSummand(H, fₖ, tors), Fₖ₋₁)
end