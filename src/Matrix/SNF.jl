using AbstractAlgebra
using AbstractAlgebra: Ring
using SparseArrays: SparseMatrixCSC, findnz

# SNF

struct SNF{R <: RingElement}
    # PAQ = S
    S  ::Vector{R}
    P  ::Union{SparseMatrixCSC{R}, Nothing}
    P⁻¹::Union{SparseMatrixCSC{R}, Nothing}
    Q  ::Union{SparseMatrixCSC{R}, Nothing}
    Q⁻¹::Union{SparseMatrixCSC{R}, Nothing}
end

function snf(A::SparseMatrixCSC{R}, baseRing::RR; flags=(false, false, false, false)) :: SNF{R} where {R <: RingElement, RR <: Ring}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A, baseRing)
    (Sᵈ, Pᵈ, Pinvᵈ, Qᵈ, Qinvᵈ) = snf_x(Aᵈ; flags=flags) # PAQ = S

    S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    isempty(S) && (S = R[])

    P = isnothing(Pᵈ) ? nothing : _toSparse(Pᵈ)
    Pinv = isnothing(Pinvᵈ) ? nothing : _toSparse(Pinvᵈ)
    Q = isnothing(Qᵈ) ? nothing : _toSparse(Qᵈ)
    Qinv = isnothing(Qinvᵈ) ? nothing : _toSparse(Qinvᵈ)

    SNF(S, P, Pinv, Q, Qinv)
end

function inverse(P::SparseMatrixCSC{R}, baseRing::RR) :: SparseMatrixCSC{R} where {R <: RingElement, RR <: Ring}
    is_one(P) && return P

    Pᵈ = _toDense(P, baseRing)
    Xᵈ = hnf_with_transform(Pᵈ)[2]
    _toSparse(Xᵈ)
end

function _toDense(A::SparseMatrixCSC{R}, baseRing::RR) :: MatrixElem where {R <: RingElement, RR <: Ring}
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

function _toSparse(Aᵈ::MatrixElem{R}) :: SparseMatrixCSC{R} where {R <: RingElement}
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