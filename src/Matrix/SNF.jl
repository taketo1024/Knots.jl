using AbstractAlgebra
using AbstractAlgebra: Ring
using SparseArrays: SparseMatrixCSC, findnz, permute
using Permutations: Permutation, inv

# SNF

mutable struct SNF{R<:RingElement}
    # PAQ = S
    S  ::Vector{R}
    P  ::Union{SparseMatrixCSC{R}, Nothing}
    P⁻¹::Union{SparseMatrixCSC{R}, Nothing}
    Q  ::Union{SparseMatrixCSC{R}, Nothing}
    Q⁻¹::Union{SparseMatrixCSC{R}, Nothing}
end

function snf(A::SparseMatrixCSC{R}, baseRing::RR; preprocess=true, flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement, RR<:Ring}
    if preprocess
        _snf_preprocess(A, baseRing; flags=flags)
    else
        _snf(A, baseRing; flags=flags)
    end
end

function _snf(A::SparseMatrixCSC{R}, baseRing::RR; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement, RR<:Ring}
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

function _snf_preprocess(A::SparseMatrixCSC{R}, baseRing::RR; flags=(false, false, false, false)) where {R<:RingElement, RR<:Ring}
    (p, q) = pivotPermutations(A)

    B = permute(A, p.data, q.data)     # B = p⁻¹ A q
    F = _snf(B, baseRing, flags=flags) # S = PBQ = (Pp⁻¹) A (qQ)

    if any(flags)
        (m, n) = size(A)
        (p⁻¹, q⁻¹) = (inv(p).data, inv(q).data)
        (p, q) = (p.data, q.data)
        (idₘ, idₙ) = (collect(1:m), collect(1:n))
            
        isnothing(F.P)   || (F.P   = permute(F.P, idₘ, p⁻¹))    # Pp⁻¹
        isnothing(F.P⁻¹) || (F.P⁻¹ = permute(F.P⁻¹, p⁻¹, idₘ))  # pP⁻¹
        isnothing(F.Q)   || (F.Q   = permute(F.Q, q⁻¹, idₙ))    # qQ
        isnothing(F.Q⁻¹) || (F.Q⁻¹ = permute(F.Q⁻¹, idₙ, q⁻¹))  # Q⁻¹q⁻¹
    end

    return F
end

function inverse(P::SparseMatrixCSC{R}, baseRing::RR) :: SparseMatrixCSC{R} where {R<:RingElement, RR<:Ring}
    is_one(P) && return P

    Pᵈ = _toDense(P, baseRing)
    Xᵈ = hnf_with_transform(Pᵈ)[2]
    _toSparse(Xᵈ)
end

function _toDense(A::SparseMatrixCSC{R}, baseRing::RR) :: MatrixElem where {R<:RingElement, RR<:Ring}
    (m, n) = size(A)

    if R <: Union{Integer, Rational}
        Aᵈ = matrix(baseRing, A)
    else
        Aᵈ = zero_matrix(baseRing, m, n)
        for (i, j, a) in zip(findnz(A)...)
            Aᵈ[i, j] = a
        end
    end
    Aᵈ
end

function _toSparse(Aᵈ::MatrixElem{R}) :: SparseMatrixCSC{R} where {R<:RingElement}
    (m, n) = size(Aᵈ)

    Is = Int[]
    Js = Int[]
    Vs = R[]

    ij = (0, 1)
    while (it = iterate(Aᵈ, ij)) !== nothing
        (r, ij) = it
        if !iszero(r) 
            push!(Is, ij[1])
            push!(Js, ij[2])
            push!(Vs, r)
        end
    end

    sparse(Is, Js, Vs, m, n)
end