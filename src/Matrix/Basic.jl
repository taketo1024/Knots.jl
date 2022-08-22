using SparseArrays
using AbstractAlgebra
using AbstractAlgebra: RingElement, Ring
using Permutations

export SparseMatrix, DenseMatrix
export print_matrix

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
    # WARN: must not call simply `sparse(Aᵈ)`, since `R` might be an AbstractAlgebra-type.

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

function permute(A::SparseMatrix{R}, p::Permutation, q::Permutation) :: SparseMatrix{R} where {R}
    SparseArrays.permute(A, p.data, q.data)
end

function permute_row(A::SparseMatrix{R}, p::Permutation) :: SparseMatrix{R} where {R}
    SparseArrays.permute(A, p.data, collect(1:size(A, 2)))
end

function permute_col(A::SparseMatrix{R}, q::Permutation) :: SparseMatrix{R} where {R}
    SparseArrays.permute(A, collect(1:size(A, 1)), q.data)
end

function print_matrix(A::AbstractMatrix)
    Base.print_matrix(stdout, A, "[", " ", "]\n")
end