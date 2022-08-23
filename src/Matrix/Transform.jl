using AbstractAlgebra: RingElement
using SparseArrays: blockdiag
using Permutations

# PAQ = B
struct Transform{M<:AbstractMatrix}
    P  ::Union{M, Nothing}
    P⁻¹::Union{M, Nothing}
    Q  ::Union{M, Nothing}
    Q⁻¹::Union{M, Nothing}
    
    function Transform{M}(P::Union{M, Nothing}, P⁻¹::Union{M, Nothing}, Q::Union{M, Nothing}, Q⁻¹::Union{M, Nothing}) where {M<:AbstractMatrix}
        @assert isnothing(P)   || size(P)[1] == size(P)[2]
        @assert isnothing(P⁻¹) || size(P⁻¹)[1] == size(P⁻¹)[2]
        @assert isnothing(Q)   || size(Q)[1] == size(Q)[2]
        @assert isnothing(Q⁻¹) || size(Q⁻¹)[1] == size(Q⁻¹)[2]
        new{M}(P, P⁻¹, Q, Q⁻¹)
    end
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

# P₁ A Q₁ = B, 
# P₂ B Q₂ = C 
# -> P₂ P₁ A Q₁ Q₂ = C
function compose(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    P = !isnothing(T₁.P) && !isnothing(T₂.P) ? 
        T₂.P * T₁.P : 
        nothing

    P⁻¹ = !isnothing(T₁.P⁻¹) && !isnothing(T₂.P⁻¹) ? 
        T₁.P⁻¹ * T₂.P⁻¹ : 
        nothing

    Q = !isnothing(T₁.Q) && !isnothing(T₂.Q) ?
        T₁.Q * T₂.Q : 
        nothing

    Q⁻¹ = !isnothing(T₁.Q⁻¹) && !isnothing(T₂.Q⁻¹) ?
        T₂.Q⁻¹ * T₁.Q⁻¹ : 
        nothing

    Transform{M}(P, P⁻¹, Q, Q⁻¹)
end

function identity_transform(::Type{SparseMatrix{R}}, size::Tuple{Int, Int}; flags=(true, true, true, true)) :: Transform{SparseMatrix{R}} where {R}
    I(k) = sparse_identity_matrix(R, k)
    (m, n) = size

    P   = flags[1] ? I(m) : nothing
    P⁻¹ = flags[2] ? I(m) : nothing
    Q   = flags[3] ? I(n) : nothing
    Q⁻¹ = flags[4] ? I(n) : nothing
    
    Transform{SparseMatrix{R}}(P, P⁻¹, Q, Q⁻¹)
end

function block_diagonal(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M} 
    P = !isnothing(T₁.P) && !isnothing(T₂.P) ? 
        blockdiag(T₁.P, T₂.P) : 
        nothing

    P⁻¹ = !isnothing(T₁.P⁻¹) && !isnothing(T₂.P⁻¹) ? 
        blockdiag(T₁.P⁻¹, T₂.P⁻¹) : 
        nothing

    Q = !isnothing(T₁.Q) && !isnothing(T₂.Q) ?
        blockdiag(T₁.Q, T₂.Q) : 
        nothing

    Q⁻¹ = !isnothing(T₁.Q⁻¹) && !isnothing(T₂.Q⁻¹) ?
        blockdiag(T₁.Q⁻¹, T₂.Q⁻¹) : 
        nothing

    Transform{M}(P, P⁻¹, Q, Q⁻¹)
end

# p⁻¹ B q = (p⁻¹P) A (Qq)
function permute(T::Transform{M}, p::Permutation, q::Permutation) :: Transform{M} where {M} 
    P = !isnothing(T.P) ? 
        permute_col(T.P, inv(p)) : 
        nothing

    P⁻¹ = !isnothing(T.P⁻¹) ? 
        permute_row(T.P⁻¹, inv(p)) : 
        nothing

    Q = !isnothing(T.Q) ?
        permute_row(T.Q, inv(q)) : 
        nothing

    Q⁻¹ = !isnothing(T.Q⁻¹) ?
        permute_col(T.Q⁻¹, inv(q)) : 
        nothing

    Transform{M}(P, P⁻¹, Q, Q⁻¹)
end