using SparseArrays
using Permutations

export SparseMatrix, permute, permute_row, permute_col, print_matrix

const DenseMatrix{R} = Base.Matrix{R}
const SparseMatrix{R} = SparseArrays.SparseMatrixCSC{R, Int}

function sparse_identity_matrix(R, n::Int) :: SparseMatrix{R}
    Is = collect(1 : n)
    Vs = fill(one(R), n)
    SparseArrays.sparse(Is, Is, Vs, n, n)
end

function density(A::SparseMatrix) :: Float64
    (m, n) = size(A)
    N = m * n

    if N > 0
        nz = count(x -> !iszero(x), findnz(A)[3])
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

function diagonal_entries(A::AbstractSparseArray{R}) :: Vector{R} where {R}
    r = minimum(size(A))
    d = fill(zero(R), r)
    for (i, j, a) in zip(findnz(A)...)
        i == j && (d[i] = a)
    end
    d
end