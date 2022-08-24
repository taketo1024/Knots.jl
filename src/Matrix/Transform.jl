using SparseArrays: blockdiag
using Permutations

const Flags4 = Tuple{Bool, Bool, Bool, Bool}

struct Transform{M<:AbstractMatrix}
    P  ::M
    P⁻¹::M
    Q  ::M
    Q⁻¹::M
    
    function Transform(P::M, P⁻¹::M, Q::M, Q⁻¹::M) where {M<:AbstractMatrix}
        @assert size(P)[1] == size(P)[2]
        @assert size(P⁻¹)[1] == size(P⁻¹)[2]
        @assert size(Q)[1] == size(Q)[2]
        @assert size(Q⁻¹)[1] == size(Q⁻¹)[2]
        new{M}(P, P⁻¹, Q, Q⁻¹)
    end
end

function identity_transform(::Type{M}, size::Tuple{Int, Int}) where {R, M<:SparseMatrix{R}}
    (m, n) = size
    I(k) = sparse_identity_matrix(R, k)
    Transform(I(m), I(m), I(n), I(n))
end

function Base.iterate(T::Transform, i = 1)
    if i == 1
        T.P, 2
    elseif i == 2
        T.P⁻¹, 3
    elseif i == 3
        T.Q, 4
    elseif i == 4
        T.Q⁻¹, 5
    else
        nothing
    end
end

function (⊕)(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    block_diagonal(T₁, T₂)
end

function Base.:*(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    compose(T₁, T₂)
end

function block_diagonal(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    P   = blockdiag(T₁.P,   T₂.P)
    P⁻¹ = blockdiag(T₁.P⁻¹, T₂.P⁻¹)
    Q   = blockdiag(T₁.Q,   T₂.Q)
    Q⁻¹ = blockdiag(T₁.Q⁻¹, T₂.Q⁻¹) 
    Transform(P, P⁻¹, Q, Q⁻¹)
end

# P₁ A Q₁ = B, 
# P₂ B Q₂ = C 
# -> P₂ P₁ A Q₁ Q₂ = C
function compose(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    P   = T₂.P   * T₁.P
    P⁻¹ = T₁.P⁻¹ * T₂.P⁻¹
    Q   = T₁.Q   * T₂.Q 
    Q⁻¹ = T₂.Q⁻¹ * T₁.Q⁻¹
    Transform(P, P⁻¹, Q, Q⁻¹)
end

# p⁻¹ B q = (p⁻¹P) A (Qq)
function permute(T::Transform{M}, p::Permutation, q::Permutation) :: Transform{M} where {M} 
    P   = permute_col(T.P, inv(p))
    P⁻¹ = permute_row(T.P⁻¹, inv(p))
    Q   = permute_row(T.Q, inv(q))
    Q⁻¹ = permute_col(T.Q⁻¹, inv(q))
    Transform(P, P⁻¹, Q, Q⁻¹)
end