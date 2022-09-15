using LinearAlgebra: I as id

function hnf_lll(A::AbstractMatrix{R}; flags=(true, true)) where {R<:Integer}
    hnf_lll!(copy(A); flags=flags)
end

function hnf_lll!(A::AbstractMatrix{R}; flags=(true, true)) where {R<:Integer}
    (m, n) = size(A)

    P    = flags[1] ? Base.Matrix(one(R) * id, m, m) : nothing
    Pinv = flags[2] ? Base.Matrix(one(R) * id, m, m) : nothing

    F = Float64
    D = Dict(i => one(F) for i in 0:m)
    λ = zeros(F, m, m)

    k = 2
    while k <= m 
        should_swap = hnf_lll_reduce(A, P, Pinv, D, λ, k, k - 1; flags=flags)
        if should_swap
            hnf_lll_swap(A, P, Pinv, D, λ, k; flags=flags)
            if k > 2 
                k -= 1
            end
        else
            for i in reverse(1:k-2)
                hnf_lll_reduce(A, P, Pinv, D, λ, k, i; flags=flags)
            end
            k += 1
        end
    end

    hnf_lll_finalize(A, P, Pinv; flags=flags)

    (A, P, Pinv)
end

function hnf_lll_reduce(A, P, Pinv, D, λ, k, i; flags::NTuple{2, Bool})
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
        flags[1] && (P[i, :] *= -1)
        flags[2] && (Pinv[:, i] *= -1)
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
        flags[1] && (P[k, :] -= q * P[i, :])
        flags[2] && (Pinv[:, i] += q * Pinv[:, k])

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

function hnf_lll_swap(A, P, Pinv, D, λ, k; flags::NTuple{2, Bool})
    (m, n) = size(A)

    (A[k, :], A[k - 1, :]) = (A[k - 1, :], A[k, :])
    flags[1] && ((P[k, :], P[k - 1, :]) = (P[k - 1, :], P[k, :]))
    flags[2] && ((Pinv[:, k], Pinv[:, k - 1]) = (Pinv[:, k - 1], Pinv[:, k]))
    
    for j in 1 : k - 2 
        (λ[k, j], λ[k - 1, j]) = (λ[k - 1, j], λ[k, j])
    end

    for i in k + 1 : m
        t = λ[i, k - 1] * D[k] - λ[i, k] * λ[k, k - 1] 
        λ[i, k - 1] = (λ[i, k - 1] * λ[k, k - 1] + λ[i, k] * D[k - 2]) / D[k - 1]
        λ[i, k] = t / D[k - 1]
    end

    D[k - 1] = (D[k - 2] * D[k] + λ[k, k - 1]^2) / D[k - 1]
end

function hnf_lll_finalize(A, P, Pinv; flags::NTuple{2, Bool})
    m = size(A, 1)
    Base.permutecols!!(A', reverse(collect(1:m)))
    flags[1] && Base.permutecols!!(P', reverse(collect(1:m)))
    flags[2] && Base.permutecols!!(Pinv, reverse(collect(1:m)))
end
