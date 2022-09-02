# `SmithNormalForm.jl` (modified, see history)
# https://github.com/wildart/SmithNormalForm.jl

module SmithNormalForm_x

using LinearAlgebra: diagm
using SparseArrays
using Base.CoreLogging
using ...Extensions: isunit, normalizing_unit

const Option{T} = Union{T, Nothing}

function snf(A::M; flags=(true, true, true, true)) :: 
    Tuple{M, Option{M}, Option{M}, Option{M}, Option{M}} where {R, M <: AbstractMatrix{R}}

    (B, P, Pinv, Q, Qinv) = init(A, flags)

    _snf_step1(B, P, Pinv, Q, Qinv, flags)
    _snf_step2(B, P, Pinv, Q, Qinv, flags)
    _snf_step3(B, P, Pinv, Q, Qinv, flags)

    (B, P, Pinv, Q, Qinv)
end

function bezout(a::R, b::R) where {R}
    (d, s, t) = gcdx(a, b) # d = sa + tb

    e = normalizing_unit(d)[1]
    if !isone(e)
        d *= e
        s *= e
        t *= e
    end

    u = div(a, d)
    isunit(u) ? (d, u, zero(R)) : (d, s, t)
end

function divisible(a::R, b::R) where {R}
    iszero(b) ? iszero(a) : iszero(a % b)
end

function rcountnz(A::M, j::Int) where {R, M <: AbstractMatrix{R}}
    c = 0
    z = zero(R)
    @inbounds for row in eachrow(A)
        if row[j] != z
            c += 1
        end
    end
    return c
end

function ccountnz(A::M, i::Int) where {R, M <: AbstractMatrix{R}}
    c = 0
    z = zero(R)
    @inbounds for col in eachcol(A)
        if col[i] != z
            c += 1
        end
    end
    return c
end

function rswap!(A::M, i1::Int, i2::Int) where {M <: AbstractMatrix}
    i1 == i2 && return A
    @inbounds for col in eachcol(A)
        col[i1], col[i2] = col[i2], col[i1]
    end
    return A
end

function cswap!(A::M, j1::Int, j2::Int) where {M <: AbstractMatrix}
    j1 == j2 && return A
    @inbounds for row in eachrow(A)
        row[j1], row[j2] = row[j2], row[j1]
    end
    return A
end

# Multiply from left: [a b] 
#                     [c d]
function rowelimination(A::M, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R, M <: AbstractMatrix{R}}
    @inbounds for col in eachcol(A)
        t = col[i]
        s = col[j]

        iszero(t) && iszero(s) && continue

        col[i] = a * t + b * s
        col[j] = c * t + d * s
    end
    return A
end

# Multiply from right: [a c] 
#                      [b d]
function colelimination(A::M, a::R, b::R, c::R, d::R, i::Int, j::Int) where {R, M <: AbstractMatrix{R}}
    @inbounds for row in eachrow(A)
        t = row[i]
        s = row[j]

        iszero(t) && iszero(s) && continue

        row[i] = a * t + b * s
        row[j] = c * t + d * s
    end
    return A
end

function select_pivot(A::M, t::Int, j::Int) :: Int where {R, M <: AbstractMatrix{R}}
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

function rmul(A::M, i::Int, a::R) where {R, M <: AbstractMatrix{R}}
    @views A[i, :] .*= a
end

function cmul(A::M, j::Int, a::R) where {R, M <: AbstractMatrix{R}}
    @views A[:, j] .*= a
end

function rowpivot(
    A::M,
    P::Option{M},
    Pinv::Option{M},
    i, j, flags) where {R, M <: AbstractMatrix{R}}

    for k in reverse!(findall(!iszero, view(A, :, j)))
        i == k && continue

        # sa + tb = d, x = a/d, y = b/d.
        #
        # [ s t][a] = [d]
        # [-y x][b]   [0]
        #
        # [ s t]⁻¹ = [x -t]
        # [-y x]     [y  s]

        a = A[i, j]
        b = A[k, j]

        d, s, t = bezout(a, b)
        x = div(a, d)
        y = div(b, d)

        rowelimination(A, s, t, -y, x, i, k)
        flags[1] && rowelimination(P, s, t, -y, x, i, k)
        flags[2] && colelimination(Pinv, x, y, -t, s, i, k)
    end
end

function colpivot(
    A::M,
    Q::Option{M},
    Qinv::Option{M},
    i, j, flags) where {R, M <: AbstractMatrix{R}}

    for k in reverse!(findall(!iszero, view(A, i, :)))
        j == k && continue

        # sa + tb = d, x = a/d, y = b/d.
        #
        # [a b][s -y] = [d 0]
        #      [t  x]
        #
        # [s -y]⁻¹ = [ x  y]
        # [t  x]     [-t  s]

        a = A[i, j]
        b = A[i, k]

        d, s, t = bezout(a, b)
        x = div(a, d)
        y = div(b, d)

        colelimination(A, s, t, -y, x, j, k)
        flags[3] && colelimination(Q, s, t, -y, x, j, k)
        flags[4] && rowelimination(Qinv, x, y, -t, s, j, k)
    end
end

function smithpivot(
    A::M,
    P::Option{M},
    Pinv::Option{M},
    Q::Option{M},
    Qinv::Option{M},
    i, j, flags) where {R, M <: AbstractMatrix{R}}

    pivot = A[i, j]
    @assert pivot != zero(R) "Pivot cannot be zero"

    while ccountnz(A, i) > 1 || rcountnz(A, j) > 1
        colpivot(A, Q, Qinv, i, j, flags)
        rowpivot(A, P, Pinv, i, j, flags)

        if pivot == A[i, j] && (ccountnz(A, i) > 1 || rcountnz(A, j) > 1)
            error("Detect endless process! pivot = ($i, $j), $(A[i, j])")
        end
    end
end

function init(A::M, flags::NTuple{4, Bool}) where {R, M <: AbstractSparseMatrix{R}}
    m, n = size(A)

    I(k) = sparse((1:k), (1:k), fill(one(R), k), k, k)

    B = copy(A)

    P    = flags[1] ? I(m) : nothing
    Pinv = flags[2] ? I(m) : nothing
    Q    = flags[3] ? I(n) : nothing
    Qinv = flags[4] ? I(n) : nothing

    return B, P, Pinv, Q, Qinv
end

function init(A::M, flags::NTuple{4, Bool}) where {R, M <: AbstractMatrix{R}}
    m, n = size(A)
    
    I(k) = diagm(fill(one(R), k))

    B = copy(A)

    P    = flags[1] ? I(m) : nothing
    Pinv = flags[2] ? I(m) : nothing
    Q    = flags[3] ? I(n) : nothing
    Qinv = flags[4] ? I(n) : nothing

    return B, P, Pinv, Q, Qinv
end

formatmtx(M) = size(M, 1) == 0 ? "[]" : repr(collect(M); context=IOContext(stdout, :compact => true))

function _snf_step1(
    A::M,
    P::Option{M},
    Pinv::Option{M},
    Q::Option{M},
    Qinv::Option{M}, 
    flags) where {R, M <: AbstractMatrix{R}}

    cols = size(A)[2]
    t = 1

    for j in 1:cols
        # @debug "Working on column $j out of $cols" D = formatmtx(D)

        rcountnz(A, j) == 0 && continue

        prow = select_pivot(A, t, j)

        # @debug "Pivot Row selected: t = $t, pivot = ($prow, $j): $(A[prow, :])"

        # swap rows
        rswap!(A, t, prow)
        flags[1] && rswap!(P, t, prow)
        flags[2] && cswap!(Pinv, t, prow)

        # swap cols
        cswap!(A, t, j)
        flags[3] && cswap!(Q, t, j)
        flags[4] && rswap!(Qinv, t, j)

        # normalize
        (u, uinv) = normalizing_unit(A[t, t])
        if !isone(u)
            cmul(A, t, u)
            flags[3] && cmul(Q, t, u)
            flags[4] && rmul(Qinv, t, uinv)
        end

        # @debug "Performing the pivot step at (i=$t, j=$t)" D = formatmtx(D)
        smithpivot(A, P, Pinv, Q, Qinv, t, t, flags)

        if issparse(A)
            dropzeros!(A)
        end

        t += 1

        # @logmsg (Base.CoreLogging.Debug - 1) "Factorization" D = formatmtx(D) U = formatmtx(U) V = formatmtx(V) U⁻¹ = formatmtx(Uinv) V⁻¹ = formatmtx(Vinv)
    end
end

function _snf_step2(
    A::M,
    P::Option{M},
    Pinv::Option{M},
    Q::Option{M},
    Qinv::Option{M}, 
    flags) where {R, M <: AbstractMatrix{R}}

    # Make sure that d_i is divisible be d_{i+1}.
    r = minimum(size(A))
    
    while true
        done = true
        for i in 1:r-1
            a = A[i, i]
            b = A[i + 1, i + 1]

            divisible(b, a) && continue

            # sa + tb = d, x = a/d, y = b/d.
            #
            # [1   1 ][a   ][s  -y] = [d      ]
            # [-ty sx][   b][t  x ]   [   ab/d]

            (d, s, t) = bezout(a, b)

            x = div(a, d)
            y = div(b, d)
            sx = s * x
            ty = t * y

            A[i, i] = d
            A[i + 1, i + 1] = x * b

            flags[1] && rowelimination(P, one(R), one(R), -ty, sx, i, i + 1)
            flags[2] && colelimination(Pinv, sx, ty, -one(R), one(R), i, i + 1)
            flags[3] && colelimination(Q, s, t, -y, x, i, i + 1)
            flags[4] && rowelimination(Qinv, x, y, -t, s, i, i + 1)

            done = false
        end
        done && break
    end
end

function _snf_step3(
    A::M,
    P::Option{M},
    Pinv::Option{M},
    Q::Option{M},
    Qinv::Option{M}, 
    flags) where {R, M <: AbstractMatrix{R}}

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
        flags[3] && cmul(Q, j, u)
        flags[4] && rmul(Qinv, j, uinv)
    end
    # @logmsg (Base.CoreLogging.Debug - 1) "Factorization" D = formatmtx(D) U = formatmtx(U) V = formatmtx(V) U⁻¹ = formatmtx(Uinv) V⁻¹ = formatmtx(Vinv)

    if issparse(A)
        dropzeros!(A)
        flags[1] && dropzeros!(P)
        flags[2] && dropzeros!(Pinv)
        flags[3] && dropzeros!(Q)
        flags[4] && dropzeros!(Qinv)
    end
end

end