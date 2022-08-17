# The following codes are copied & modified from AbstractAlgebra.
#
# hnf: https://github.com/Nemocas/AbstractAlgebra.jl/blob/0a2868547f71d2cf8b48b8da347276a183804b9e/src/Matrix.jl#L4480
# snf: https://github.com/Nemocas/AbstractAlgebra.jl/blob/0a2868547f71d2cf8b48b8da347276a183804b9e/src/Matrix.jl#L5123
#
# The individual .jl files in the AbstractAlgebra package are
# licensed under the Simplified "2-clause" BSD License.
# https://github.com/Nemocas/AbstractAlgebra.jl/blob/master/LICENSE.md

using AbstractAlgebra
using AbstractAlgebra: RingElement, MatrixElem, kb_search_first_pivot, swap_rows!, div

function hnf_x(A::MatrixElem{T}; flags=(true, true)) where {T <: RingElement}
    H = deepcopy(A)
    m = nrows(H)
    P =    flags[1] ? identity_matrix(A, m) : nothing
    Pinv = flags[2] ? identity_matrix(A, m) : nothing

    hnf_kb_x!(H, P, Pinv)
    
    return H, P, Pinv
end 

function snf_x(A::MatrixElem{T}; flags=(true, true, true, true)) where {T <: RingElement}
    S = deepcopy(A)
    m = nrows(S)
    n = ncols(S)

    P =    flags[1] ? identity_matrix(A, m) : nothing
    Pinv = flags[2] ? identity_matrix(A, m) : nothing
    Q =    flags[3] ? identity_matrix(A, n) : nothing
    Qinv = flags[4] ? identity_matrix(A, n) : nothing

    snf_kb_x!(S, P, Pinv, Q, Qinv)
    
    return S, P, Pinv, Q, Qinv
end

# private (HNF)

function hnf_kb_x!(H, P, Pinv, start_element::Int = 1)
    m = nrows(H)
    n = ncols(H)
    pivot = zeros(Int, n) # pivot[j] == i if the pivot of column j is in row i

    # Find the first non-zero entry of H
    row1, col1 = kb_search_first_pivot(H, start_element)
    if row1 == 0
        return nothing
    end
    pivot[col1] = row1
    kb_canonical_row_x!(H, P, Pinv, row1, col1)
    pivot_max = col1
    t = base_ring(H)()
    t1 = base_ring(H)()
    t2 = base_ring(H)()
    for i = row1 + 1:m
        new_pivot = false
        for j = start_element:n
            if iszero(H[i, j])
                continue
            end
            if pivot[j] == 0
                # We found a non-zero entry in a column without a pivot: This is a
                # new pivot
                pivot[j] = i
                pivot_max = max(pivot_max, j)
                new_pivot = true
            else
                # We have a pivot for this column: Use it to write 0 in H[i, j]

                # MEMO: (x, y) = (H[p, j], H[i, j])
                # ux + vy = d, a = x/d, b = -y/d.
                #
                # [u v][x] = [d]
                # [b a][y]   [0]
                #   ^ 
                #  det = 1
                
                p = pivot[j]
                d, u, v = gcdx(H[p, j], H[i, j])
                a = divexact(H[p, j], d)
                b = -divexact(H[i, j], d)

                # updte H
                for c = j:n
                    t = deepcopy(H[i, c])
                    t1 = mul_red!(t1, a, H[i, c], false)
                    t2 = mul_red!(t2, b, H[p, c], false)
                    H[i, c] = reduce!(t1 + t2)
                    t1 = mul_red!(t1, u, H[p, c], false)
                    t2 = mul_red!(t2, v, t, false)
                    H[p, c] = reduce!(t1 + t2)
                end

                # update P
                !isnothing(P) && for c = 1:m
                    t = deepcopy(P[i, c])
                    t1 = mul_red!(t1, a, P[i, c], false)
                    t2 = mul_red!(t2, b, P[p, c], false)
                    P[i, c] = reduce!(t1 + t2)
                    t1 = mul_red!(t1, u, P[p, c], false)
                    t2 = mul_red!(t2, v, t, false)
                    P[p, c] = reduce!(t1 + t2)
                end
                
                # update P⁻¹
                # multiply 
                #
                # [u v]⁻¹ = [a  -v]
                # [b a]     [-b  u]
                #
                # from right.

                !isnothing(Pinv) && for c = 1:m
                    t = deepcopy(Pinv[c, i])
                    t1 = mul_red!(t1, u, Pinv[c, i], false)
                    t2 = mul_red!(t2, v, Pinv[c, p], false)
                    Pinv[c, i] = reduce!(t1 - t2)
                    t1 = mul_red!(t1, a, Pinv[c, p], false)
                    t2 = mul_red!(t2, b, t, false)
                    Pinv[c, p] = reduce!(t1 - t2)
                end
            end

            # We changed the pivot of column j (or found a new one).
            # We have do reduce the entries marked with # in
            # ( 0 0 0 . * )
            # ( . # # * * )
            # ( 0 0 . * * )
            # ( 0 . # * * )
            # ( * * * * * )
            # where . are pivots and i = 4, j = 2. (This example is for the
            # "new pivot" case.)
            kb_canonical_row_x!(H, P, Pinv, pivot[j], j)
            for c = j:pivot_max
                if pivot[c] == 0
                    continue
                end
                kb_reduce_column_x!(H, P, Pinv, pivot, c, start_element)
            end
            if new_pivot
                break
            end
        end
    end
    kb_sort_rows_x!(H, P, Pinv, pivot, start_element)
    return nothing
end

function kb_canonical_row_x!(H, P, Pinv, r::Int, c::Int)
    m = nrows(H)
    n = ncols(H)

    cu = canonical_unit(H[r, c])
    if !isone(cu)
        # update H
        for j = c:n
            H[r, j] = divexact(H[r, j], cu)
        end

        # update P
        !isnothing(P) && for j = 1:m
            P[r, j] = divexact(P[r, j], cu)
        end

        # update P⁻¹.
        # D(r, r; u)⁻¹ = D(r, r; u⁻¹).
        # (multiply D(r, r; u⁻¹) from right) == (col[r] *= u⁻¹)
        !isnothing(Pinv) && for j = 1:m
            Pinv[j, r] = Pinv[j, r] * cu
        end
    end
    return nothing
end

function kb_reduce_column_x!(H, P, Pinv, pivot::Vector{Int}, c::Int, start_element::Int = 1)
    # Let c = 4 and pivot[c] = 4. H could look like this:
    # ( 0 . * # * )
    # ( . * * # * ) < p = pivot[i] (i = 1)
    # ( 0 0 0 0 . )
    # ( 0 0 0 . * ) < r = pivot[c]
    # ( * * * * * )
    #         ^ c
    # (. are pivots, we want to reduce the entries marked with #)
    # The #'s are in rows whose pivot is in a column left of column c.

    m = nrows(H)
    n = ncols(H)
    r = pivot[c]
    t = base_ring(H)()

    for i = start_element:c - 1
        p = pivot[i]
        if p == 0
            continue
        end
        
        # So, the pivot in row p is in a column left of c.
        if iszero(H[p, c])
            continue
        end

        q = -div(H[p, c], H[r, c])

        # update H
        for j = c:n
            t = mul!(t, q, H[r, j])
            H[p, j] += t
        end

        # update P
        !isnothing(P) && for j = 1:m
            t = mul!(t, q, P[r, j])
            P[p, j] += t
        end

        # update P⁻¹
        # A(p, r; q)⁻¹ = A(p, r; -q).
        # (multiply A(p, r; -q) from right) == (col[r] += -q * col[p])
        !isnothing(Pinv) && for j = 1:m
            t = mul!(t, q, Pinv[j, p])
            Pinv[j, r] -= t
        end
    end
    return nothing
end

function kb_sort_rows_x!(H, P, Pinv, pivot::Vector{Int}, start_element::Int = 1)

    m = nrows(H)
    n = ncols(H)
    pivot2 = zeros(Int, m)
    for i = 1:n
        if pivot[i] == 0
            continue
        end
        pivot2[pivot[i]] = i
    end

    r1 = start_element
    for i = start_element:n
        r2 = pivot[i]
        if r2 == 0
            continue
        end
        if r1 != r2
            # update H
            swap_rows!(H, r1, r2)

            # update P
            !isnothing(P) && swap_rows!(P, r1, r2)

            # update P⁻¹
            !isnothing(Pinv) && swap_cols!(Pinv, r1, r2)

            p = pivot2[r1]
            pivot[i] = r1
            if p != 0
                pivot[p] = r2
            end
            pivot2[r1] = i
            pivot2[r2] = p
        end
        r1 += 1
        if r1 == m
            break
        end
    end
    return nothing
end

# private (SNF)

function snf_kb_x!(S, P, Pinv, Q, Qinv)
    snf_kb_x_step1!(S, P, Pinv, Q, Qinv)
    snf_kb_x_step2!(S, P, Pinv, Q, Qinv)
    return nothing
end

function snf_kb_x_step1!(S, P, Pinv, Q, Qinv)
    m = nrows(S)
    n = ncols(S)
    l = min(m, n)
    i = 1
    while i <= l
        kb_clear_row_x!(S, Q, Qinv, i)
        hnf_kb_x!(S, P, Pinv, i)
        c = i + 1
        while c <= n && iszero(S[i, c])
            c += 1
        end
        if c != n + 1
            continue
        end
        i+=1
    end
end

function snf_kb_x_step2!(S, P, Pinv, Q, Qinv)
    m = nrows(S)
    n = ncols(S)
    l = min(m, n)
    t = base_ring(S)()
    t1 = base_ring(S)()
    t2 = base_ring(S)()
    for i = 1:l-1
        for j = i + 1:l
            if isone(S[i, i])
                break
            end
            if iszero(S[i, i]) && iszero(S[j, j])
                continue
            end
            
            # MEMO: 
            # (x, y) = (S[i, i], S[j, j]).
            # ux + vy = d, a = x/d, b = -y/d.
            # ua - vb = 1
            #
            # [1  1 ][x  ][u b] = [d      ]
            # [vb ua][  y][v a]   [   xy/d]

            d, u, v = gcdx(S[i, i], S[j, j])
            a = divexact(S[i, i], d)
            b = -divexact(S[j, j], d)

            # update S
            S[j, j] = divexact(S[i, i]*S[j, j], d)
            S[i, i] = d
            
            # update P
            t1 = mul!(t1, b, v)
            !isnothing(P) && for c = 1:m
                t = deepcopy(P[i, c])
                P[i, c] += P[j, c]
                t2 = mul_red!(t2, t1, P[j, c], false)
                P[j, c] += t2
                t2 = mul_red!(t2, t1, t, false)
                P[j, c] = reduce!(P[j, c] + t2)
            end

            # update P⁻¹
            # multiply
            #
            # [1  1 ]⁻¹ = [ ua -1]
            # [vb ua]     [-vb  1]
            #
            # from right.
            !isnothing(Pinv) && for c = 1:m
                t = deepcopy(Pinv[c, i])
                t2 = mul_red!(t2, t1, t, false)
                Pinv[c, i] += t2
                t2 = mul_red!(t2, t1, Pinv[c, j], false)
                Pinv[c, i] = reduce!(Pinv[c, i] - t2)
                Pinv[c, j] -= t
            end

            # update Q
            !isnothing(Q) && for r = 1:n
                t = deepcopy(Q[r, i])
                t1 = mul_red!(t1, Q[r, i], u, false)
                t2 = mul_red!(t2, Q[r, j], v, false)
                Q[r, i] = reduce!(t1 + t2)
                t1 = mul_red!(t1, t, b, false)
                t2 = mul_red!(t2, Q[r, j], a, false)
                Q[r, j] = reduce!(t1 + t2)
            end

            # update Q⁻¹
            # multiply
            #
            # [u b]⁻¹ = [a  -b]
            # [v a]     [-v  u]
            #
            # from left.
            !isnothing(Qinv) && for r = 1:n
                t = deepcopy(Qinv[i, r])
                t1 = mul_red!(t1, Qinv[i, r], a, false)
                t2 = mul_red!(t2, Qinv[j, r], b, false)
                Qinv[i, r] = reduce!(t1 - t2)
                t1 = mul_red!(t1, t, v, false)
                t2 = mul_red!(t2, Qinv[j, r], u, false)
                Qinv[j, r] = reduce!(t2 - t1)
            end
        end
    end
end

function kb_clear_row_x!(S, Q, Qinv, i::Int)
    m = nrows(S)
    n = ncols(S)
    t = base_ring(S)()
    t1 = base_ring(S)()
    t2 = base_ring(S)()
    for j = i+1:n
        if iszero(S[i, j])
            continue
        end

        # MEMO: (x, y) = (S[i, i], S[i, j])
        # ux + vy = d, a = x/d, b = -y/d.
        #
        # [x y][u b] = [d]
        #      [v a]   [0]
        #        ^ 
        #       det = 1

        d, u, v = gcdx(S[i, i], S[i, j])
        a = divexact(S[i ,i], d)
        b = -divexact(S[i, j], d)

        # update S
        for r = i:m
            t = deepcopy(S[r, j])
            t1 = mul_red!(t1, a, S[r, j], false)
            t2 = mul_red!(t2, b, S[r, i], false)
            S[r, j] = reduce!(t1 + t2)
            t1 = mul_red!(t1, u, S[r, i], false)
            t2 = mul_red!(t2, v, t, false)
            S[r, i] = reduce!(t1 + t2)
        end

        # update Q
        !isnothing(Q) && for r = 1:n
            t = deepcopy(Q[r,j])
            t1 = mul_red!(t1, a, Q[r, j], false)
            t2 = mul_red!(t2, b, Q[r, i], false)
            Q[r, j] = reduce!(t1 + t2)
            t1 = mul_red!(t1, u, Q[r, i], false)
            t2 = mul_red!(t2, v, t, false)
            Q[r, i] = reduce!(t1 + t2)
        end

        # update Q⁻¹.
        # mutiply
        #
        # [u b]⁻¹ = [a  -b]
        # [v a]     [-v  u]
        #
        # from left.
        !isnothing(Qinv) && for r = 1:n
            t = deepcopy(Qinv[j, r])
            t1 = mul_red!(t1, u, Qinv[j, r], false)
            t2 = mul_red!(t2, -v, Qinv[i, r], false)
            Qinv[j, r] = reduce!(t1 + t2)
            t1 = mul_red!(t1, a, Qinv[i, r], false)
            t2 = mul_red!(t2, -b, t, false)
            Qinv[i, r] = reduce!(t1 + t2)
        end
    end
    return nothing
end