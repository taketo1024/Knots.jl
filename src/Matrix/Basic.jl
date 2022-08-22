using SparseArrays
using AbstractAlgebra
using AbstractAlgebra: RingElement, Ring
using Permutations
using ..Env

export SparseMatrix, DenseMatrix
export print_matrix

const SparseMatrix = SparseArrays.SparseMatrixCSC
const DenseMatrix = AbstractAlgebra.MatrixElem

function _toDense(A::SparseMatrix{R}) :: DenseMatrix{R} where {R<:RingElement, RR<:Ring}
    baseRing = Env.get_base_ring(R)
    (m, n) = size(A)
    Aᵈ = AbstractAlgebra.zero_matrix(baseRing, m, n)
    for (i, j, a) in zip(SparseArrays.findnz(A)...)
        Aᵈ[i, j] = a
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


function sparse_identity_matrix(R, n::Int) :: SparseMatrix{R}
    Is = collect(1 : n)
    Vs = fill(one(R), n)
    SparseArrays.sparse(Is, Is, Vs, n, n)
end

function is_identity(A::SparseMatrix) :: Bool
    for (i, j, a) in zip(SparseArrays.findnz(A)...)
        if i == j && isone(a)
            continue
        elseif i != j && iszero(a)
            continue
        else
            return false
        end
    end
    true
end

function is_zero(A::SparseMatrix) :: Bool
    for a in SparseArrays.findnz(A)[3]
        if !iszero(a)
            return false
        end
    end
    true
end

function density(A::SparseMatrix) :: Float64
    (m, n) = size(A)
    N = m * n

    if N > 0
        nz = count(x -> !iszero(x), SparseArrays.findnz(A)[3])
        nz / N
    else 
        0
    end
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