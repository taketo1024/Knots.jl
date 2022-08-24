using SparseArrays: blockdiag, sparse_hvcat, spdiagm
using LinearAlgebra: UnitUpperTriangular

function schur_complement(A::SparseMatrix{R}, piv::Pivot{R}; flags=(false, false, false, false)) where {R}
    I(k) = sparse_identity_matrix(R, k)

    (m, n) = size(A)
    (p, q, r) = permutations(piv)

    B = permute(A, p, q) # B = p⁻¹ A q

    # make left-upper of B, unittriangular

    d = fill(one(R), n)

    for i in 1 : r 
        u = B[i, i]
        isone(u) && continue

        B[:, i] .*= u # col-ops are faster
        d[i] = u
    end

    T0 = Transform{SparseMatrix{R}}(
        flags[1] ? I(m) : nothing,
        flags[2] ? I(m) : nothing,
        flags[3] ? spdiagm(d) : nothing,
        flags[4] ? spdiagm(d) : nothing
    )
    (S, T1) = _schur_complement_U(B, r, flags)
    T = permute( compose(T0, T1), p, q)

    (r, S, T)
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

    U = UnitUpperTriangular(A[1:r, 1:r])
    Uinv = sparse(inv(U))

    X = A[1 : r, r + 1 : n]
    Y = A[r + 1 : m, 1 : r]
    Z = A[r + 1 : m, r + 1 : n]
    S = Z - Y * Uinv * X

    P = flags[1] ? 
        sparse_hvcat(
            (2, 2), 
            Uinv, O(r, m - r),
            -Y * Uinv, I(m - r)
        ) : nothing

    Pinv = flags[2] ? 
        sparse_hvcat(
            (2, 2), 
            U, O(r, m - r),
            Y, I(m - r)
        ) : nothing

    Q = flags[3] ? 
        sparse_hvcat(
            (2, 2), 
            I(r), -Uinv * X,
            O(n - r, r), I(n - r)
        ) : nothing

    Qinv = flags[4] ?
        sparse_hvcat(
            (2, 2), 
            I(r), Uinv * X,
            O(n - r, r), I(n - r)
        ) : nothing

    T = Transform{SparseMatrix{R}}(P, Pinv, Q, Qinv)
    (S, T)
end