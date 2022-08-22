using AbstractAlgebra
using SparseArrays: blockdiag

export SNF, snf

# SNF

mutable struct SNF{R<:RingElement}
    # PAQ = S
    S  ::Vector{R}
    P  ::Union{SparseMatrix{R}, Nothing}
    P⁻¹::Union{SparseMatrix{R}, Nothing}
    Q  ::Union{SparseMatrix{R}, Nothing}
    Q⁻¹::Union{SparseMatrix{R}, Nothing}
end

function snf(A::SparseMatrix{R}; preprocess=true, flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    d_threshold = 0.75

    if iszero(A)
        _zero(R, size(A); flags=flags)
    elseif preprocess && density(A) < d_threshold
        _snf_preprocess(A; flags=flags)
    else
        _snf_dense_sort(A; flags=flags)
    end
end

function _snf_preprocess(A::SparseMatrix{R}; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    piv = pivot(A)
    npivots(piv) == 0 && return _snf_dense_sort(A; flags=flags)

    (S, r, P, Pinv, Q, Qinv) = schur_complement(A, piv, flags=flags)

    if min(size(S)...) > 0 
        next = snf(S; flags=flags)
        _snf_compose(r, P, Pinv, Q, Qinv, next)
    else
        SNF(fill(one(R), r), P, Pinv, Q, Qinv)
    end
end

function _snf_compose(r, P, Pinv, Q, Qinv, next::SNF{R}) where {R} 
    I(k) = sparse_identity_matrix(R, k)

    d = vcat(fill(one(R), r), next.S)
    
    if !isnothing(P)
        P = blockdiag(I(r), next.P) * P
    end

    if !isnothing(Pinv)
        Pinv = Pinv * blockdiag(I(r), next.P⁻¹)
    end

    if !isnothing(Q)
        Q = Q * blockdiag(I(r), next.Q)
    end

    if !isnothing(Qinv)
        Qinv = blockdiag(I(r), next.Q⁻¹) * Qinv
    end

    SNF(d, P, Pinv, Q, Qinv)
end

function _snf_dense_sort(A::SparseMatrix{R}; flags=(false, false, false, false)) :: SNF{R} where {R}
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
        _snf_dense(A; flags=flags)
    else
        p = permutation(collect(rows), m)
        q = permutation(collect(cols), n)
        B = permute(A, p, q)[1:k, 1:l]

        F = _snf_dense(B; flags=flags)

        d = F.S
        P    = isnothing(F.P)   ?  nothing : permute_col(blockdiag(F.P, I(m - k)), inv(p))
        Pinv = isnothing(F.P⁻¹) ?  nothing : permute_row(blockdiag(F.P⁻¹, I(m - k)), inv(p))
        Q    = isnothing(F.Q)   ?  nothing : permute_row(blockdiag(F.Q, I(n - l)), inv(q))
        Qinv = isnothing(F.Q⁻¹) ?  nothing : permute_col(blockdiag(F.Q⁻¹, I(n - l)), inv(q))
    
        SNF(d, P, Pinv, Q, Qinv)
    end
end

function _snf_dense(A::SparseMatrix{R}; flags=(false, false, false, false)) :: SNF{R} where {R<:RingElement}
    (m, n) = size(A)
    r = min(m, n)

    Aᵈ = _toDense(A)
    (Sᵈ, Pᵈ, Pinvᵈ, Qᵈ, Qinvᵈ) = snf_x(Aᵈ; flags=flags) # PAQ = S

    S = filter(!iszero, map( i -> Sᵈ[i, i], 1 : r ))
    isempty(S) && (S = R[])

    P    = isnothing(Pᵈ)    ?  nothing : _toSparse(Pᵈ)
    Pinv = isnothing(Pinvᵈ) ?  nothing : _toSparse(Pinvᵈ)
    Q    = isnothing(Qᵈ)    ?  nothing : _toSparse(Qᵈ)
    Qinv = isnothing(Qinvᵈ) ?  nothing : _toSparse(Qinvᵈ)

    SNF(S, P, Pinv, Q, Qinv)
end



function _zero(R::Type{RR}, size::Tuple{Int, Int}; flags=(false, false, false, false)) :: SNF{R} where {RR}
    I(k) = sparse_identity_matrix(R, k)
    (m, n) = size

    d = R[]
    P    = flags[1] ? I(m) : nothing
    Pinv = flags[2] ? I(m) : nothing
    Q    = flags[3] ? I(n) : nothing
    Qinv = flags[4] ? I(n) : nothing

    SNF(d, P, Pinv, Q, Qinv)
end     