using AbstractAlgebra: RingElement, Ring

abstract type AbstractHomology{R<:RingElement, RR<:Ring} end

# must override 
function complex(H::AbstractHomology{R, RR}) :: AbstractComplex{R, RR} where {R, RR}
    throw(MethodError(complex, (H,)))
end

# must override 
function makeSummand(H::AbstractHomology{R, RR}, rank::Int, torsions::Vector{R}) :: AbstractHomologySummand{R, RR} where {R, RR}
    throw(MethodError(_makeSummand, (H, rank, torsions)))
end

# can override (e.g. use different algorithm / cache result etc)
function compute(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    compute_simple(H, k)
end

function baseRing(H::AbstractHomology{R, RR}) :: RR where {R, RR <: Ring}
    baseRing(complex(H))
end

function hDegRange(H::AbstractHomology) :: UnitRange{Int}
    hDegRange(complex(H))
end

function Base.getindex(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    compute(H, k)
end

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

function compute_simple(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    r = baseRing(H)
    C = complex(H)
    deg = differentialDegree(C)

    nₖ = length(generators(C, k))
    nₖ == 0 && return makeSummand(H, 0, R[])

    Aₖ₋₁ = differential(C, k - deg)
    Aₖ   = differential(C, k)

    Fₖ₋₁ = snf(Aₖ₋₁, r)
    Fₖ   = snf(Aₖ, r)

    rₖ₋₁ = length(Fₖ₋₁.S)
    rₖ   = length(Fₖ.S)
    
    fₖ = nₖ - rₖ₋₁ - rₖ
    tors = filter(r -> !is_unit(r), Fₖ₋₁.S)

    makeSummand(H, fₖ, tors)
end

function asString(H::AbstractHomology) :: String
    lines = String[]
    for i in hDegRange(H)
        Hi = asString(H[i])
        push!(lines, "H[$i] = $Hi")
    end
    join(lines, "\n")
end

function Base.show(io::IO, H::AbstractHomology) 
    print(io, asString(H))
end