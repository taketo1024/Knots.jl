function inv_upper_triangular(U::SparseMatrix{R}) :: SparseMatrix{R} where {R}
    n = size(U, 1)

    zero = Base.zero(R)
    one = Base.one(R)

    Is = Int[]
    Js = Int[]
    Vs = R[]

    x = fill(zero, n)
    e = fill(zero, n)
    d = diagonal_entries(U)

    for j in 1:n
        e[j] = one

        _solve_upper_trianguler!(U, d, e, x)

        for (i, a) in enumerate(x)
            iszero(a) && continue

            push!(Is, i)
            push!(Js, j)
            push!(Vs, a)
        end

        # it is guaranteed that e = 0.
        fill!(x, zero)
    end

    sparse(Is, Js, Vs, n, n)
end

function _solve_upper_trianguler!(U::SparseMatrix{R}, d::Vector{R}, b::Vector{R}, x::Vector{R}) :: Vector{R} where {R}
    n = size(U, 1)

    rows = rowvals(U)
    vals = nonzeros(U)

    for j in reverse(1:n)
        u = d[j] # must be unit
        v = x[j] = div(b[j], u)
        for k in nzrange(U, j)
            i = rows[k]
            a = vals[k]
            b[i] -= a * v
        end
    end

    x
end