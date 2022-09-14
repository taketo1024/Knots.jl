using LinearAlgebra: I as id

function hnf_lll(G::AbstractMatrix{R}) where {R<:Integer}
    (m, n) = size(G)

    A = copy(G)
    B = Base.Matrix(one(R) * id, m, m)

    F = Float64
    D = Dict(i => one(F) for i in 0:m)
    λ = zeros(F, m, m)

    # hnf_lll_init(A, B)

    k = 2
    while k <= m 
        # println("k = $k")
        # print_matrix(A)
        # print_matrix(B)
        # println(D)
        # print_matrix(λ)

        should_swap = hnf_lll_reduce(A, B, D, λ, k, k - 1)
        if should_swap
            hnf_lll_swap(A, B, D, λ, k)
            if k > 2 
                k -= 1
            end
        else
            for i in reverse(1:k-2)
                hnf_lll_reduce(A, B, D, λ, k, i)
            end
            k += 1
        end
    end

    (A, B)
end

function hnf_lll_init(A, B)
    (m, n) = size(A)
    for j in 1:n
        e = count(!iszero, view(A, :, j))
        e == 0 && continue

        if count(!iszero, view(A, :, j)) == 1 && A[m, j] < 0
            A[m, :] *= -1
            B[m, m]  = -1
        end
        break
    end
end

function hnf_lll_reduce(A, B, D, λ, k, i)
    (m, n) = size(A)

    function _findnz_inrow(i)
        for j in 1:n 
            !iszero(A[i, j]) && return j
        end
        return n + 1
    end

    col1 = _findnz_inrow(i)
    if col1 <= n && A[i, col1] < 0 
        hnf_lll_minus(λ, m, i)
        A[i, :] *= -1
        B[i, :] *= -1
    end

    col2 = _findnz_inrow(k)

    q = if col1 <= n
        floor(Int, A[k, col1] / A[i, col1])
    elseif 2 * abs(λ[k, i]) > D[i]
        round(Int, λ[k, i] / D[i])
    else
        0
    end

    if !iszero(q)
        A[k, :] -= q * A[i, :]
        B[k, :] -= q * B[i, :]
        λ[k, i] -= q * D[i]
        for j in 1 : i - 1
            λ[k, j] -= q * λ[i, j]
        end
    end

    return col1 <= min(col2, n) || 
        (col1 == col2 == n + 1 && 4(D[k - 2] * D[k] + λ[k, k-1]^2) < 3 * D[k - 1]^2)
end

function hnf_lll_minus(λ, m, j)
    for r in 2:m
        for s in 1:r-1
            if r == j || s == j 
                λ[r, s] *= -1
            end
        end
    end
end

function hnf_lll_swap(A, B, D, λ, k)
    (m, n) = size(A)

    (A[k, :], A[k - 1, :]) = (A[k - 1, :], A[k, :])
    (B[k, :], B[k - 1, :]) = (B[k - 1, :], B[k, :])
    
    for j in 1 : k - 2 
        (λ[k, j], λ[k - 1, j]) = (λ[k - 1, j], λ[k, j])
    end

    for i in k + 1 : m
        t = λ[i, k - 1] * D[k] - λ[i, k] * λ[k, k - 1] 
        λ[i, k - 1] = (λ[i, k - 1] * λ[k, k - 1] + λ[i, k] * D[k - 2]) / D[k - 1]
        λ[i, k] = t / D[k - 1]
    end

    D[k - 1] = (D[k - 2] * D[k] + λ[k, k - 1]^2 / D[k - 1])
end
