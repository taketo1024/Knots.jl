using SparseArrays: blockdiag, sparse_hvcat, spdiagm
using LinearAlgebra: UnitUpperTriangular

function schur_complement(A::SparseMatrix{R}, piv::Pivot{R}; flags=(true, true, true, true)) :: Tuple{SparseMatrix{R}, Transform{SparseMatrix{R}}} where {R}
    (m, n) = size(A)

    r = npivots(piv)
    (p, q) = permutations(piv)

    B = permute(A, p, q) # B = p⁻¹ A q

    (S, T) = _schur_complement_U(B, r, flags)

    if any(flags)
        T = permute(T, p, q)
    end

    (S, T)
end

function _schur_complement_U(A::SparseMatrix{R}, r::Int, flags) where {R}
    # A = [U X] ~> [I  ]
    #     [Y Z]    [  S]
    #
    # by 
    #
    # [U⁻¹    ] [U X] [I -U⁻¹X] = [I  ]
    # [-YU⁻¹ I] [Y Z] [    I  ]   [  S]

    I(k) = sparse_identity_matrix(R, k)
    O(k, l) = spzeros(R, k, l)

    (m, n) = size(A)

    U = A[1:r, 1:r]
    Uinv = inv_upper_triangular(U)

    X = A[1 : r, r + 1 : n]
    Y = A[r + 1 : m, 1 : r]
    Z = A[r + 1 : m, r + 1 : n]
    S = Z - Y * Uinv * X

    P = flags[1] ? 
        sparse_hvcat(
            (2, 2), 
            Uinv, O(r, m - r),
            -Y * Uinv, I(m - r)
        ) : I(m)

    Pinv = flags[2] ? 
        sparse_hvcat(
            (2, 2), 
            U, O(r, m - r),
            Y, I(m - r)
        ) : I(m)

    Q = flags[3] ? 
        sparse_hvcat(
            (2, 2), 
            I(r), -Uinv * X,
            O(n - r, r), I(n - r)
        ) : I(n)

    Qinv = flags[4] ?
        sparse_hvcat(
            (2, 2), 
            I(r), Uinv * X,
            O(n - r, r), I(n - r)
        ) : I(n)

    T = Transform(P, Pinv, Q, Qinv)
    
    (S, T)
end