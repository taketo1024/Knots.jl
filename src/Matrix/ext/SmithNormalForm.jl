# `SmithNormalForm.jl` (modified, see history)
# https://github.com/wildart/SmithNormalForm.jl

module SmithNormalForm_x

using LinearAlgebra
using SparseArrays
using Base.CoreLogging
using ...Extensions: isunit, normalizing_unit

function snf(A::AbstractMatrix{R}; inverse=true) where {R}
    (Pinv, Qinv, B, P, Q) = init(A, inverse=inverse)

    _snf_step1(Pinv, Qinv, B, P, Q; inverse=inverse)
    _snf_step2(Pinv, Qinv, B, P, Q; inverse=inverse)
    _snf_step3(Pinv, Qinv, B, P, Q; inverse=inverse)

    (Pinv, Qinv, B, P, Q)
end

function bezout(a::R, b::R) where {R}
    (g, s, t) = gcdx(a, b)

    if g == a
        s = one(R)
        t = zero(R)
    elseif g == -a
        s = -one(R)
        t = zero(R)
    end
    
    (s, t, g)
end

function divisable(a::R, b::R) where {R}
    b == zero(R) && return a == zero(R)
    return iszero(a % b)
end

function divide(a::R, b::R) where {R}
    if b != -one(R)
        return div(a, b)
    else
        return a * b
    end
end

function rcountnz(A::AbstractMatrix{R}, j::Int) where {R}
    c = 0
    z = zero(R)
    @inbounds for row in eachrow(A)
        if row[j] != z
            c += 1
        end
    end
    return c
end

function ccountnz(A::AbstractMatrix{R}, i::Int) where {R}
    c = 0
    z = zero(R)
    @inbounds for col in eachcol(A)
        if col[i] != z
            c += 1
        end
    end
    return c
end

function rswap!(A::AbstractMatrix, i1::Int, i2::Int)
    i1 == i2 && return A
    @inbounds for col in eachcol(A)
        col[i1], col[i2] = col[i2], col[i1]
    end
    return A
end

function cswap!(A::AbstractMatrix, j1::Int, j2::Int)
    j1 == j2 && return A
    @inbounds for row in eachrow(A)
        row[j1], row[j2] = row[j2], row[j1]
    end
    return A
end

function rowelimination(A::AbstractMatrix{R}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    @inbounds for col in eachcol(A)
        t = col[i]
        s = col[j]
        col[i] = a * t + b * s
        col[j] = c * t + d * s
    end
    return A
end

function colelimination(A::AbstractMatrix{R}, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R}
    @inbounds for row in eachrow(A)
        t = row[i]
        s = row[j]
        row[i] = a * t + b * s
        row[j] = c * t + d * s
    end
    return A
end

function select_pivot(A::AbstractMatrix{R}, t::Int, j::Int) :: Int where {R}
    # Good pivot row for j-th column is the one
    # that have a smallest number of elements
    rows = size(A)[1]
    prow = 0
    rsize = typemax(Int)
    for i in t:rows
        iszero(A[i, j]) && continue
        c = count(!iszero, view(A, i, :))
        if c < rsize
            rsize = c
            prow = i
        end
    end
    (prow > 0) ? prow : error()
end

function rmul(A::AbstractMatrix{R}, i::Int, a::R) where {R}
    @views A[i, :] .*= a
end

function cmul(A::AbstractMatrix{R}, j::Int, a::R) where {R}
    @views A[:, j] .*= a
end

function rowpivot(
    Pinv::AbstractMatrix{R},
    P::AbstractMatrix{R},
    A::AbstractMatrix{R},
    i, j; inverse=true) where {R}

    for k in reverse!(findall(!iszero, view(A, :, j)))
        a = A[i, j]
        b = A[k, j]

        i == k && continue

        s, t, g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        rowelimination(A, s, t, -y, x, i, k)
        inverse && rowelimination(P, s, t, -y, x, i, k)
        colelimination(Pinv, x, y, -t, s, i, k)
    end
end

function colpivot(
    Qinv::AbstractMatrix{R},
    Q::AbstractMatrix{R},
    A::AbstractMatrix{R},
    i, j; inverse=true) where {R}

    for k in reverse!(findall(!iszero, view(A, i, :)))
        a = A[i, j]
        b = A[i, k]

        j == k && continue

        s, t, g = bezout(a, b)
        x = divide(a, g)
        y = divide(b, g)

        colelimination(A, s, t, -y, x, j, k)
        inverse && colelimination(Q, s, t, -y, x, j, k)
        rowelimination(Qinv, x, y, -t, s, j, k)
    end
end

function smithpivot(
    Pinv::AbstractMatrix{R},
    P::AbstractMatrix{R},
    Qinv::AbstractMatrix{R},
    Q::AbstractMatrix{R},
    A::AbstractMatrix{R},
    i, j; inverse=true) where {R}

    pivot = A[i, j]
    @assert pivot != zero(R) "Pivot cannot be zero"
    while ccountnz(A, i) > 1 || rcountnz(A, j) > 1
        colpivot(Qinv, Q, A, i, j, inverse=inverse)
        rowpivot(Pinv, P, A, i, j, inverse=inverse)
    end
end

function init(A::AbstractSparseMatrix{R,Ti}; inverse=true) where {R,Ti}
    B = copy(A)
    rows, cols = size(A)

    Pinv = spzeros(R, rows, rows)
    for i in 1:rows
        Pinv[i, i] = one(R)
    end
    P = inverse ? copy(Pinv) : spzeros(R, 0, 0)

    Qinv = spzeros(R, cols, cols)
    for i in 1:cols
        Qinv[i, i] = one(R)
    end
    Q = inverse ? copy(Qinv) : spzeros(R, 0, 0)

    return Pinv, Qinv, B, P, Q
end

function init(A::AbstractMatrix{R}; inverse=true) where {R}
    B = copy(A)
    rows, cols = size(A)

    Pinv = zeros(R, rows, rows)
    for i in 1:rows
        Pinv[i, i] = one(R)
    end
    P = inverse ? copy(Pinv) : zeros(R, 0, 0)

    Qinv = zeros(R, cols, cols)
    for i in 1:cols
        Qinv[i, i] = one(R)
    end
    Q = inverse ? copy(Qinv) : zeros(R, 0, 0)

    return Pinv, Qinv, B, P, Q
end

formatmtx(M) = size(M, 1) == 0 ? "[]" : repr(collect(M); context=IOContext(stdout, :compact => true))

function _snf_step1(
    Pinv::AbstractMatrix{R},
    Qinv::AbstractMatrix{R},
    A::AbstractMatrix{R},
    P::AbstractMatrix{R},
    Q::AbstractMatrix{R}
    ; inverse=true) where {R}

    cols = size(A)[2]
    t = 1

    for j in 1:cols
        # @debug "Working on column $j out of $cols" D = formatmtx(D)

        rcountnz(A, j) == 0 && continue

        prow = select_pivot(A, t, j)

        @debug "Pivot Row selected: t = $t, pivot = ($prow, $j): $(A[prow, :])"

        # swap rows
        rswap!(A, t, prow)
        inverse && rswap!(P, t, prow)
        cswap!(Pinv, t, prow)

        # swap cols
        cswap!(A, t, j)
        inverse && cswap!(Q, t, j)
        rswap!(Qinv, t, j)

        # normalize
        (u, uinv) = normalizing_unit(A[t, t])
        if !isone(u)
            cmul(A, t, u)
            rmul(Qinv, t, uinv)
            inverse && cmul(Q, t, u)
        end

        # @debug "Performing the pivot step at (i=$t, j=$t)" D = formatmtx(D)
        smithpivot(Pinv, P, Qinv, Q, A, t, t, inverse=inverse)

        if issparse(A)
            dropzeros!(A)
        end

        t += 1

        # @logmsg (Base.CoreLogging.Debug - 1) "Factorization" D = formatmtx(D) U = formatmtx(U) V = formatmtx(V) U⁻¹ = formatmtx(Uinv) V⁻¹ = formatmtx(Vinv)
    end
end

function _snf_step2(
    Pinv::AbstractMatrix{R},
    Qinv::AbstractMatrix{R},
    A::AbstractMatrix{R},
    P::AbstractMatrix{R},
    Q::AbstractMatrix{R}
    ; inverse=true) where {R}

    # Make sure that d_i is divisible be d_{i+1}.
    r = minimum(size(A))
    pass = true
    while pass
        pass = false
        for i in 1:r-1
            divisable(A[i+1, i+1], A[i, i]) && continue
            pass = true
            A[i+1, i] = A[i+1, i+1]

            colelimination(Q, one(R), one(R), zero(R), one(R), i, i + 1)
            rowelimination(Qinv, one(R), zero(R), -one(R), one(R), i, i + 1)

            smithpivot(Pinv, P, Qinv, Q, A, i, i, inverse=inverse)
        end
    end
end

function _snf_step3(
    Pinv::AbstractMatrix{R},
    Qinv::AbstractMatrix{R},
    A::AbstractMatrix{R},
    P::AbstractMatrix{R},
    Q::AbstractMatrix{R}
    ; inverse=true) where {R}

    rows, cols = size(A)

    # To guarantee SNFⱼ = Λⱼ ≥ 0 we absorb the sign of Λ into T and T⁻¹, s.t.
    #    Λ′ = Λ*sign(Λ),   T′ = sign(Λ)*T,    and    T⁻¹′ = T⁻¹*sign(Λ),
    # with the convention that sign(0) = 1. Then we still have that X = SΛT = SΛ′T′
    # and also that Λ = S⁻¹XT⁻¹ ⇒ Λ′ = S⁻¹XT⁻¹′.
    for j in 1:rows
        j > cols && break
        d = A[j, j]

        (u, uinv) = normalizing_unit(d)
        isone(u) && continue

        A[j, j] *= u
        rmul(Qinv, j, uinv)
        inverse && cmul(Q, j, u)
    end
    # @logmsg (Base.CoreLogging.Debug - 1) "Factorization" D = formatmtx(D) U = formatmtx(U) V = formatmtx(V) U⁻¹ = formatmtx(Uinv) V⁻¹ = formatmtx(Vinv)

    if issparse(A)
        return dropzeros!(Pinv), dropzeros!(Qinv), dropzeros!(A), dropzeros!(P), dropzeros!(Q)
    else
        return Pinv, Qinv, A, P, Q
    end
end

end