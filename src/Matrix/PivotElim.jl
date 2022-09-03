function pivotal_elim(A::M; flags::Flags4, itr=1) :: Tuple{M, Int, Transform{M}, Permutation, Permutation} where {R, M <: SparseMatrix{R}}
    @debug "pivotal-elim (step $itr)" A = size(A) density = density(A)

    (m, n) = size(A)
    I(k, l) = identity_transform(M, (k, l); flags=flags)

    piv = pivot(A)
    r = npivots(piv)
    (p, q) = permutations(piv)

    if r == 0 
        return (A, 0, I(m, n), p, q)
    end

    (S, T) = schur_complement(A, piv; flags=flags)

    if !iszero(S)
        (S₂, r₂, T₂, p₂, q₂) = pivotal_elim(S; flags=flags, itr=itr+1)

        @debug "compose results (step $itr)" S₂ = size(S₂) S₁ = size(S)
    
        S = S₂
        T = T * (I(r, r) ⊕ T₂)
        p = p * shift(p₂, r)
        q = q * shift(q₂, r)
        r += r₂
    end

    (itr == 1) && @debug "pivotal-elim done." A = size(A) S = size(S) total_pivots = r
    
    (S, r, T, p, q)
end

# Instead of constructing the transform matrix, 
# this method directly applies the transformation to the given matrix.

function pivotal_elim_with(A::SparseMatrix{R}, v::Vector{R}; left::Bool, itr=1) :: Tuple{SparseMatrix{R}, Vector{R}, Int, Permutation, Permutation} where {R}
    @debug "pivotal-elim-with (step $itr)" A = size(A) density = density(A)
    @assert left && length(v) == size(A, 1) || !left && length(v) == size(A, 2)

    (m, n) = size(A)

    piv = pivot(A)
    r = npivots(piv)
    (p, q) = permutations(piv)

    if r == 0 
        return (A, v, 0, p, q)
    end

    flags = left ? (true, false, false, false) : (false, false, false, true)

    (S, T) = schur_complement(A, piv; flags=flags)

    if left 
        v = (T.P * v)[r + 1 : m]
    else 
        v = (T.Q⁻¹ * v)[r + 1 : n]
    end

    if !iszero(S)
        (S, v, r₂, p₂, q₂) = pivotal_elim_with(S, v; left=left, itr=itr+1)

        p = p * shift(p₂, r)
        q = q * shift(q₂, r)
        r = r + r₂
    end

    (itr == 1) && @debug "pivotal-elim-with done." A = size(A) S = size(S) total_pivots = r
    
    (S, v, r, p, q)
end
