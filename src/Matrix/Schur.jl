using SparseArrays: blockdiag, sparse_hvcat, spdiagm
using LinearAlgebra: UnitUpperTriangular

function schur_complement(A::SparseMatrix{R}, piv::Pivot{R}; flags=(false, false, false, false)) where {R}
    I(k) = sparse_identity_matrix(R, k)
    O(k, l) = spzeros(R, k, l)

    (m, n) = size(A)
    (p, q, r) = permutations(piv)
    B = permute(A, p, q) # B = p⁻¹ A q

    # make left-upper of B, unittriangular

    d = fill(one(R), n)

    for i in 1 : r 
        u = B[i, i]
        isone(u) && continue

        B[:, i] *= u # col-ops are faster
        d[i] = u
    end

    # B = [U X] ~> [I  ]
    #     [Y Z]    [  S]
    #
    # by 
    #
    # [U⁻¹    ] [U X] [I -U⁻¹X] = [I  ]
    # [-YU⁻¹ I] [Y Z] [    I  ]   [  S]

    U = UnitUpperTriangular(B[1:r, 1:r])
    Uinv = sparse(inv(U))

    X = B[1 : r, r + 1 : n]
    Y = B[r + 1 : m, 1 : r]
    Z = B[r + 1 : m, r + 1 : n]
    S = Z - Y * Uinv * X

    if flags[1]
        P = sparse_hvcat(
            (2, 2), 
            Uinv, O(r, m - r),
            -Y * Uinv, I(m - r)
        )
        P = permute_col(P, inv(p))
    else
        P = nothing
    end

    if flags[2]
        Pinv = sparse_hvcat(
            (2, 2), 
            U, O(r, m - r),
            Y, I(m - r)
        )
        Pinv = permute_row(Pinv, inv(p))
    else
        Pinv = nothing
    end

    if flags[3]
        Q0 = spdiagm(d)
        Q1 = sparse_hvcat(
            (2, 2), 
            I(r), -Uinv * X,
            O(n - r, r), I(n - r)
        )
        Q = permute_row(Q0 * Q1, inv(q))
    else
        Q = nothing
    end

    if flags[4]
        Q0inv = spdiagm(d)
        Q1inv = sparse_hvcat(
            (2, 2), 
            I(r), Uinv * X,
            O(n - r, r), I(n - r)
        )
        Qinv = permute_col(Q1inv * Q0inv, inv(q))
    else
        Qinv = nothing
    end

    (S, r, P, Pinv, Q, Qinv)
end
