using AbstractAlgebra
using AbstractAlgebra: Ring
using SparseArrays: SparseMatrixCSC

# SNF

struct SNF{R <: RingElement}
    # PAQ = S
    S::Vector{R}
    P::SparseMatrixCSC{R}
    Q::SparseMatrixCSC{R}
end

function snf(A::SparseMatrixCSC{R}, baseRing::RR; withTrans=false) :: SNF{R} where {R <: RingElement, RR <: Ring}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A, baseRing)

    if withTrans 
        (Sᵈ, Pᵈ, Qᵈ) = AbstractAlgebra.snf_with_transform(Aᵈ) # PAQ = S

        S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
        isempty(S) && (S = R[])

        P = _toSparse(Pᵈ)
        Q = _toSparse(Qᵈ)

        SNF(S, P, Q)
    else
        Sᵈ = AbstractAlgebra.snf(Aᵈ)

        S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
        isempty(S) && (S = R[])

        P = spzeros(R, m, m)
        Q = spzeros(R, n, n)

        SNF(S, P, Q)
    end
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