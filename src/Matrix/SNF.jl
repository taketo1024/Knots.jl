include("ext/SmithNormalForm.jl")

using SparseArrays: blockdiag

export SNF, snf

# SNF

struct SNF{R}
    factors::Vector{R}
    T::Transform{SparseMatrix{R}}
end

function Base.iterate(S::SNF{R}, i = 0) where {R}
    if i == 0
        S.factors, 1
    elseif i == 1
        S.T.P, 2
    elseif i == 2
        S.T.P⁻¹, 3
    elseif i == 3
        S.T.Q, 4
    elseif i == 4
        S.T.Q⁻¹, 5
    else
        nothing
    end
end

function snf(A::SparseMatrix{R}; preprocess=true, flags::Flags4=(true, true, true, true)) :: SNF{R} where {R}
    @debug "snf A: $(size(A)), density: $(density(A))"
    d_threshold = 0.5

    if iszero(A)
        SNF(R[], identity_transform(SparseMatrix{R}, size(A)))
    elseif preprocess
        _snf_preprocess(A, flags)
    elseif density(A) < d_threshold
        _snf_sparse(A, flags)
    else
        _snf_dense(A, flags)
    end
end

function _snf_preprocess(A::SparseMatrix{R}, flags::Flags4) :: SNF{R} where {R}
    @debug "snf-preprocess A: $(size(A)), density: $(density(A))"

    piv = pivot(A)
    npivots(piv) == 0 && return snf(A; preprocess=false, flags=flags)

    (r, S, T) = schur_complement(A, piv; flags=flags)

    if min(size(S)...) > 0 
        next = snf(S; flags=flags)
        _snf_compose(r, T, next)
    else
        d = fill(one(R), r)
        SNF(d, T)
    end
end

function _snf_sparse(A::SparseMatrix{R}, flags::Flags4) :: SNF{R} where {R}
    @debug "snf-sparse A: $(size(A)), density: $(density(A))"

    r = min(size(A)...)
    (S, P, Pinv, Q, Qinv) = SmithNormalForm_x.snf(A; flags=flags)

    d = filter(!iszero, map( i -> S[i, i], 1 : r ))
    isempty(d) && (d = R[])

    T = Transform(P, Pinv, Q, Qinv)

    SNF(d, T)
end

function _snf_dense(A::SparseMatrix{R}, flags::Flags4) :: SNF{R} where {R}
    @debug "snf-dense A: $(size(A)), density: $(density(A))"
    I(k) = sparse_identity_matrix(R, k)
    (m, n) = size(A)

    rows = Set{Int}()
    cols = Set{Int}()

    for (i, j) in zip(findnz(A)...)
        push!(rows, i)
        push!(cols, j)
    end

    (k, l) = (length(rows), length(cols))

    if (k, l) == (m, n)
        _snf_dense_sorted(A, flags)
    else
        p = permutation(collect(rows), m)
        q = permutation(collect(cols), n)
        B = permute(A, p, q)[1:k, 1:l]

        F = _snf_dense_sorted(B, flags)

        d = F.factors
        I = identity_transform(SparseMatrix{R}, (m - k, n - l))
        T = permute(F.T ⊕ I, p, q)
    
        SNF(d, T)
    end
end

function _snf_dense_sorted(A::SparseMatrix{R}, flags::Flags4) :: SNF{R} where {R}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = DenseMatrix{R}(A)
    (Sᵈ, Pᵈ, Pinvᵈ, Qᵈ, Qinvᵈ) = SmithNormalForm_x.snf(Aᵈ; flags=flags) # PAQ = S

    d = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    isempty(d) && (d = R[])

    T = Transform(
        SparseMatrix{R}(Pᵈ),
        SparseMatrix{R}(Pinvᵈ),
        SparseMatrix{R}(Qᵈ),
        SparseMatrix{R}(Qinvᵈ)
    )

    SNF(d, T)
end

function _snf_compose(r, T1, next::SNF{R}) where {R} 
    I(k) = sparse_identity_matrix(R, k)

    d = vcat(fill(one(R), r), next.factors)
    I = identity_transform(SparseMatrix{R}, (r, r))
    T = T1 * (I ⊕ next.T)

    SNF(d, T)
end