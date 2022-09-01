using SparseArrays: sparse_vcat, spzeros
using ..Extensions: isunit
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

function compute_single(H::AbstractHomology{R}, k::Int; preprocess=true, with_transform=false) :: AbstractHomologySummand{R} where {R}
    @debug "compute H[$k]" preprocess=preprocess with_transform=with_transform

    C = complex(H)
    deg = differentialDegree(C)
    p = preprocess

    nₖ = length(generators(C, k))
    nₖ == 0 && return makeSummand(H, k, 0, R[], spzeros(R, 0, 0))

    Aₖ₋₁ = differential(C, k - deg)

    flgₖ₋₁ = (with_transform, true, false, false) # need Pₖ₋₁, Pₖ₋₁⁻¹
    Fₖ₋₁ = snf(Aₖ₋₁; preprocess=p, flags=flgₖ₋₁)
    rₖ₋₁ = length(Fₖ₋₁.factors)

    Aₖ = differential(C, k)
    Tₖ = Fₖ₋₁.T.P⁻¹[:, rₖ₋₁ + 1 : nₖ] # restricts Aₖ, size = (nₖ, nₖ - rₖ₋₁)
    Bₖ = (rₖ₋₁ > 0) ? Aₖ * Tₖ : Aₖ

    flgₖ = (false, false, false, with_transform) # only need Qₖ⁻¹
    Fₖ = snf(Bₖ; preprocess=p, flags=flgₖ)
    rₖ = length(Fₖ.factors)
    
    fₖ = nₖ - rₖ₋₁ - rₖ

    tors = filter(r -> !isunit(r), Fₖ₋₁.factors)
    tₖ = length(tors)

    trans = if with_transform
        Pₖ₋₁ = Fₖ₋₁.T.P                      # size = (nₖ, nₖ)
        T₁ = Pₖ₋₁[rₖ₋₁ - tₖ + 1 : rₖ₋₁, :]    # tor part, size = (tₖ, nₖ)
        T₂ = Pₖ₋₁[rₖ₋₁ + 1 : nₖ, :]           # domain of Bₖ, size = (nₖ - rₖ₋₁, nₖ)
        Z  = Fₖ.T.Q⁻¹[rₖ + 1 : nₖ - rₖ₋₁, :]  # kernel of Bₖ, size = (fₖ, nₖ - rₖ₋₁)
        T₃ = Z * T₂                          # free part, size = (fₖ, nₖ)

        sparse_vcat(T₃, T₁)
    else
        spzeros(R, fₖ + tₖ, nₖ)
    end

    makeSummand(H, k, fₖ, tors, trans)
end