using SparseArrays
using AbstractAlgebra
using AbstractAlgebra: RingElement, Ring

const SparseMatrix = SparseArrays.SparseMatrixCSC
const DenseMatrix = AbstractAlgebra.MatrixElem

function _toDense(A::SparseMatrix{R}, baseRing::RR) :: DenseMatrix{R} where {R<:RingElement, RR<:Ring}
    (m, n) = size(A)

    if R <: Union{Integer, Rational}
        Aᵈ = AbstractAlgebra.matrix(baseRing, A)
    else
        Aᵈ = AbstractAlgebra.zero_matrix(baseRing, m, n)
        for (i, j, a) in zip(SparseArrays.findnz(A)...)
            Aᵈ[i, j] = a
        end
    end
    Aᵈ
end

function _toSparse(Aᵈ::DenseMatrix{R}) :: SparseMatrix{R} where {R<:RingElement}
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

    SparseArrays.sparse(Is, Js, Vs, m, n)
end