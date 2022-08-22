using AbstractAlgebra
using SparseArrays: blockdiag

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

function snf(A::SparseMatrix{R}; preprocess=true, flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    if preprocess && density(A) < 0.5
        _snf_preprocess(A; flags=flags)
    else
        _snf(A; flags=flags)
    end
end

function _snf(A::SparseMatrix{R}; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A)
    (Sᵈ, Pᵈ, Pinvᵈ, Qᵈ, Qinvᵈ) = snf_x(Aᵈ; flags=flags) # PAQ = S

    S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    isempty(S) && (S = R[])

    P    = isnothing(Pᵈ)    ?  nothing : _toSparse(Pᵈ)
    Pinv = isnothing(Pinvᵈ) ?  nothing : _toSparse(Pinvᵈ)
    Q    = isnothing(Qᵈ)    ?  nothing : _toSparse(Qᵈ)
    Qinv = isnothing(Qinvᵈ) ?  nothing : _toSparse(Qinvᵈ)

    SNF(S, P, Pinv, Q, Qinv)
end

function _snf_preprocess(A::SparseMatrix{R}; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    piv = pivot(A)
    npivots(piv) == 0 && return _snf(A; flags=flags)

    (S, r, P, Pinv, Q, Qinv) = schur_complement(A, piv, flags=flags)

    if min(size(S)...) > 0 
        next = snf(S; flags=flags)
        _snf_compose(r, P, Pinv, Q, Qinv, next)
    else
        SNF(fill(one(R), r), P, Pinv, Q, Qinv)
    end
end

function _snf_compose(r, P, Pinv, Q, Qinv, next::SNF{R}) where {R} 
    I(k) = sparse_identity_matrix(R, k)

    d = vcat(fill(one(R), r), next.S)
    
    if !isnothing(P)
        P = blockdiag(I(r), next.P) * P
    end

    if !isnothing(Pinv)
        Pinv = Pinv * blockdiag(I(r), next.P⁻¹)
    end

    if !isnothing(Q)
        Q = Q * blockdiag(I(r), next.Q)
    end

    if !isnothing(Qinv)
        Qinv = blockdiag(I(r), next.Q⁻¹) * Qinv
    end

    if !isnothing(P) && !isnothing(Pinv)
        @assert is_one(P * Pinv)
    end

    if !isnothing(Q) && !isnothing(Qinv)
        @assert is_one(Q * Qinv)
    end

    SNF(d, P, Pinv, Q, Qinv)
end