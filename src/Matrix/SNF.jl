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

function snf_zero(A::SparseMatrix{R}) :: SNF{R} where {R}
    SNF(R[], identity_transform(SparseMatrix{R}, size(A)))
end

function snf(A::SparseMatrix{R}; preprocess=true, flags::Flags4=(true, true, true, true)) :: SNF{R} where {R}
    d_threshold = 0.5

    if iszero(A)
        snf_zero(A)
    elseif preprocess
        snf_with_preprocess(A; flags=flags)
    elseif density(A) < d_threshold
        snf_sparse(A; flags=flags)
    else
        snf_dense(A; flags=flags)
    end
end

function snf_with_preprocess(A::SparseMatrix{R}; flags::Flags4) :: SNF{R} where {R}
    (F, S) = snf_preprocess(A; flags=flags)
    
    if !iszero(S)
        F₂ = snf(S; preprocess=false, flags=flags)
        snf_compose(F, F₂)
    else
        F
    end
end

function snf_preprocess(A::SparseMatrix{R}; flags::Flags4, itr=1) :: Tuple{SNF{R}, SparseMatrix{R}, Permutation, Permutation} where {R}
    @debug "snf-preprocess (step $itr)" A = size(A) density = density(A)

    piv = pivot(A)
    r = npivots(piv)
    (p, q) = permutations(piv)

    if r == 0 
        return (snf_zero(A), A, p, q)
    end

    (S, T) = schur_complement(A, piv; flags=flags)

    d = fill(one(R), r)
    F = SNF(d, T)

    if !iszero(S)
        (F₂, S₂, p₂, q₂) = snf_preprocess(S; flags=flags, itr=itr+1)

        F = snf_compose(F, F₂)
        S = S₂
        p = p * shift(p₂, r)
        q = q * shift(q₂, r)
    end

    if itr == 1
        @debug "snf-preprocess done, total_pivots = $(length(F.factors))"
    end
    
    (F, S, p, q)
end

function snf_sparse(A::SparseMatrix{R}; flags::Flags4) :: SNF{R} where {R}
    @debug "snf-sparse" A = size(A) density = density(A)

    r = min(size(A)...)
    (S, P, Pinv, Q, Qinv) = SmithNormalForm_x.snf(A; flags=flags)

    d = filter(!iszero, map( i -> S[i, i], 1 : r ))
    isempty(d) && (d = R[])

    T = Transform(P, Pinv, Q, Qinv)

    SNF(d, T)
end

function snf_dense(A::SparseMatrix{R}; flags::Flags4) :: SNF{R} where {R}
    @debug "snf-dense" A = size(A) density = density(A)

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
        _snf_dense_sorted(A; flags=flags)
    else
        p = permutation(collect(rows), m)
        q = permutation(collect(cols), n)
        B = permute(A, p, q)[1:k, 1:l]

        F = _snf_dense_sorted(B; flags=flags)

        d = F.factors
        I = identity_transform(SparseMatrix{R}, (m - k, n - l))
        T = permute(F.T ⊕ I, p, q)
    
        SNF(d, T)
    end
end

function _snf_dense_sorted(A::SparseMatrix{R}; flags::Flags4) :: SNF{R} where {R}
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

function snf_compose(F1::SNF{R}, F2::SNF{R}) :: SNF{R} where {R}
    if length(F2.factors) == 0
        return F1
    end

    I(k) = sparse_identity_matrix(R, k)

    r = length(F1.factors)
    d = vcat(F1.factors, F2.factors)
    I = identity_transform(SparseMatrix{R}, (r, r))
    T = F1.T * (I ⊕ F2.T)

    SNF(d, T)
end

function shift(p::Permutation, s::Int) :: Permutation
   indices = append!( Array(1:s), map(i -> i + s, p.data) )
   Permutation(indices)
end