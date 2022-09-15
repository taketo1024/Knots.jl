using Base.Threads

function inv_upper_triangular(U::SparseMatrix{R}) :: SparseMatrix{R} where {R}
    if Threads.nthreads() >= 8 && minimum(size(U)) >= 10000
        _inv_upper_triangular_m(U)
    else
        _inv_upper_triangular_s(U)
    end
end

function _inv_upper_triangular_s(U::SparseMatrix{R}) :: SparseMatrix{R} where {R}
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

        for i in 1:j
            a = x[i]
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

function _inv_upper_triangular_m(U::SparseMatrix{R}) :: SparseMatrix{R} where {R}
    n = size(U, 1)

    zero = Base.zero(R)
    one = Base.one(R)
    d = diagonal_entries(U)

    Is = Int[]
    Js = Int[]
    Vs = R[]

    l = ReentrantLock()
    resources = Vector(undef, Threads.nthreads())

    @threads for j in 1:n
        tid = Threads.threadid()
        if isassigned(resources, tid)
            (x, e) = resources[tid]
        else
            x = fill(zero, n)
            e = fill(zero, n)
            resources[tid] = (x, e)
        end
    
        e[j] = one

        _solve_upper_trianguler!(U, d, e, x)
        
        for i in 1:j
            a = x[i]
            iszero(a) && continue

            lock(l) do 
                push!(Is, i)
                push!(Js, j)
                push!(Vs, a)
            end
        end

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
            if !iszero(a)
                b[i] -= a * v
            end
        end
    end

    x
end