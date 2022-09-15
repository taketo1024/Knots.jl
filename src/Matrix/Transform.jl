using SparseArrays: blockdiag
using Permutations

const Flags4 = Tuple{Bool, Bool, Bool, Bool}

struct Transform{M<:AbstractMatrix}
    P  ::Union{M, Nothing}
    P⁻¹::Union{M, Nothing}
    Q  ::Union{M, Nothing}
    Q⁻¹::Union{M, Nothing}
    flags::Flags4
    
    function Transform(::Type{M}, P::Union{M, Nothing}, P⁻¹::Union{M, Nothing}, Q::Union{M, Nothing}, Q⁻¹::Union{M, Nothing}) where {M<:AbstractMatrix}
        flags = (!isnothing(P), !isnothing(P⁻¹), !isnothing(Q), !isnothing(Q⁻¹))

        @assert !flags[1] || size(P)[1]   == size(P)[2]
        @assert !flags[2] || size(P⁻¹)[1] == size(P⁻¹)[2]
        @assert !flags[3] || size(Q)[1]   == size(Q)[2]
        @assert !flags[4] || size(Q⁻¹)[1] == size(Q⁻¹)[2]

        @assert !(flags[1] && flags[2]) || size(P) == size(P⁻¹)
        @assert !(flags[3] && flags[4]) || size(Q) == size(Q⁻¹)

        new{M}(P, P⁻¹, Q, Q⁻¹, flags)
    end
end

function identity_transform(m_type::Type{M}, size::Tuple{Int, Int}; flags=(true, true, true, true)) where {R, M<:SparseMatrix{R}}
    (m, n) = size
    I(k) = sparse_identity_matrix(R, k)
    Transform(
        m_type,
        flags[1] ? I(m) : nothing, 
        flags[2] ? I(m) : nothing, 
        flags[3] ? I(n) : nothing, 
        flags[4] ? I(n) : nothing 
    )
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

function (⊕)(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M <: SparseMatrix} 
    block_diagonal(T₁, T₂)
end

function Base.:*(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M <: AbstractMatrix} 
    compose(T₁, T₂)
end

function block_diagonal(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M <: SparseMatrix} 
    P   = (T₁.flags[1] && T₂.flags[1]) ? 
        blockdiag(T₁.P, T₂.P) : 
        nothing
    P⁻¹ = (T₁.flags[2] && T₂.flags[2]) ? 
        blockdiag(T₁.P⁻¹, T₂.P⁻¹) : 
        nothing
    Q   = (T₁.flags[3] && T₂.flags[3]) ? 
        blockdiag(T₁.Q, T₂.Q) : 
        nothing
    Q⁻¹ = (T₁.flags[4] && T₂.flags[4]) ? 
        blockdiag(T₁.Q⁻¹, T₂.Q⁻¹) : 
        nothing
    Transform(M, P, P⁻¹, Q, Q⁻¹)
end

# P₁ A Q₁ = B, 
# P₂ B Q₂ = C 
# -> P₂ P₁ A Q₁ Q₂ = C
function compose(T₁::Transform{M}, T₂::Transform{M}) :: Transform{M} where {M <: AbstractMatrix} 
    P   = (T₁.flags[1] && T₂.flags[1]) ? 
        T₂.P * T₁.P : 
        nothing
    P⁻¹ = (T₁.flags[2] && T₂.flags[2]) ? 
        T₁.P⁻¹ * T₂.P⁻¹ : 
        nothing
    Q   = (T₁.flags[3] && T₂.flags[3]) ? 
        T₁.Q   * T₂.Q :
        nothing
    Q⁻¹ = (T₁.flags[4] && T₂.flags[4]) ? 
        T₂.Q⁻¹ * T₁.Q⁻¹ :
        nothing
    Transform(M, P, P⁻¹, Q, Q⁻¹)
end

# p⁻¹ B q = (p⁻¹P) A (Qq)
function permute(T::Transform{M}, p::Permutation, q::Permutation) :: Transform{M} where {M <: SparseMatrix} 
    P   = T.flags[1] ?
        permute_col(T.P, inv(p)) : 
        nothing
    P⁻¹ = T.flags[2] ?
        permute_row(T.P⁻¹, inv(p)) : 
        nothing
    Q   = T.flags[3] ?
        permute_row(T.Q, inv(q)) : 
        nothing
    Q⁻¹ = T.flags[4] ?
        permute_col(T.Q⁻¹, inv(q)) : 
        nothing
    Transform(M, P, P⁻¹, Q, Q⁻¹)
end

Base.convert(::Type{T}, a::T) where {T <: Transform} = a
Base.convert(::Type{Transform{M}}, T::Transform) where {M<:AbstractMatrix} = 
    Transform(
        M,
        isnothing(T.P)   ? nothing : convert(M, T.P),
        isnothing(T.P⁻¹) ? nothing : convert(M, T.P⁻¹),
        isnothing(T.Q)   ? nothing : convert(M, T.Q),
        isnothing(T.Q⁻¹) ? nothing : convert(M, T.Q⁻¹),
    )
