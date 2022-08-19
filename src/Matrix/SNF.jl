using AbstractAlgebra
using AbstractAlgebra: Ring

export SNF, snf

# SNF

mutable struct SNF{R<:RingElement}
    # PAQ = S
    S  ::Vector{R}
    P  ::Union{SparseMatrix{R}, Nothing}
    P⁻¹::Union{SparseMatrix{R}, Nothing}
    Q  ::Union{SparseMatrix{R}, Nothing}
    Q⁻¹::Union{SparseMatrix{R}, Nothing}
end

function snf(A::SparseMatrix{R}, baseRing::RR; preprocess=true, flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement, RR<:Ring}
    if preprocess
        _snf_preprocess(A, baseRing; flags=flags)
    else
        _snf(A, baseRing; flags=flags)
    end
end

function _snf(A::SparseMatrix{R}, baseRing::RR; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement, RR<:Ring}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A, baseRing)
    (Sᵈ, Pᵈ, Pinvᵈ, Qᵈ, Qinvᵈ) = snf_x(Aᵈ; flags=flags) # PAQ = S

    S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    isempty(S) && (S = R[])

    P    = isnothing(Pᵈ)    ?  nothing : _toSparse(Pᵈ)
    Pinv = isnothing(Pinvᵈ) ?  nothing : _toSparse(Pinvᵈ)
    Q    = isnothing(Qᵈ)    ?  nothing : _toSparse(Qᵈ)
    Qinv = isnothing(Qinvᵈ) ?  nothing : _toSparse(Qinvᵈ)

    SNF(S, P, Pinv, Q, Qinv)
end

function _snf_preprocess(A::SparseMatrix{R}, baseRing::RR; flags=(false, false, false, false)) where {R<:RingElement, RR<:Ring}
    (p, q) = pivotPermutations(A)

    B = permute(A, p, q)     # B = p⁻¹ A q
    F = _snf(B, baseRing, flags=flags) # S = PBQ = (Pp⁻¹) A (qQ)

    if any(flags)
        isnothing(F.P)   || (F.P   = permute_col(F.P, inv(p)))    # Pp⁻¹
        isnothing(F.P⁻¹) || (F.P⁻¹ = permute_row(F.P⁻¹, inv(p)))  # pP⁻¹
        isnothing(F.Q)   || (F.Q   = permute_row(F.Q, inv(q)))    # qQ
        isnothing(F.Q⁻¹) || (F.Q⁻¹ = permute_col(F.Q⁻¹, inv(q)))  # Q⁻¹q⁻¹
    end

    return F
end