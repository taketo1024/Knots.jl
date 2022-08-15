using AbstractAlgebra
using SparseArrays: SparseMatrixCSC

# SNF

struct SNF{R <: RingElement}
    diag::Vector{R}

    # PAQ = S
    P::SparseMatrixCSC{R}
    P⁻¹::SparseMatrixCSC{R}
    Q::SparseMatrixCSC{R}
    Q⁻¹::SparseMatrixCSC{R}
end

function snf(A::SparseMatrixCSC{R}, ref::R) :: SNF{R} where {R <: RingElement}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A, ref)
    (Sᵈ, Pᵈ, Qᵈ) = AbstractAlgebra.snf_with_transform(Aᵈ) # PAQ = S

    diag = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    P = _toSparse(Pᵈ)
    Q = _toSparse(Qᵈ)
    P⁻¹ = _toSparse(inv(Pᵈ))
    Q⁻¹ = _toSparse(inv(Qᵈ))

    SNF(diag, P, P⁻¹, Q, Q⁻¹)
end

function _toDense(A::SparseMatrixCSC{R}, ref::R) :: AbstractAlgebra.MatrixElem where {R <: RingElement}
    (m, n) = size(A)
    P = parent(ref)
    
    if R <: Union{Integer, Rational}
        Aᵈ = AbstractAlgebra.matrix(P, A)
    else
        Aᵈ = AbstractAlgebra.zero_matrix(P, m, n)
        for (i, j, a) in zip(findnz(A)...)
            Aᵈ[i, j] = a
        end
    end
    Aᵈ
end

function _toSparse(Aᵈ::AbstractAlgebra.MatrixElem{R}) :: SparseMatrixCSC{R} where {R <: RingElement}
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