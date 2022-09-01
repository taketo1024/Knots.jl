using SparseArrays: sparse
using Permutations
using ..Matrix: snf_preprocess

abstract type AbstractComplex{R, X} end

function hDegRange(C::AbstractComplex) :: UnitRange{Int}
    throw(MethodError(hDegRange, (C,)))
end

function generators(C::AbstractComplex{R, X}, k::Int) :: Vector{X} where {R, X}
    throw(MethodError(generators, (C, k)))
end

function drop_generators!(C::AbstractComplex, k::Int, r::Int, p::Permutation)
    throw(MethodError(drop_generators!, (C, k, r, p)))
end

# override if necessary
function original_generators(C::AbstractComplex{R, X}, k::Int) :: Vector{X} where {R, X}
    generators(C, k)
end

function differentialDegree(C::AbstractComplex) :: Int
    +1
end

function differentiate(C::AbstractComplex{R, X}, k::Int, x::X) :: Vector{Pair{X, R}} where {R, X}
    throw(MethodError(differentiate, (C, k, x)))
end

# override if necessary
function differential(C::AbstractComplex{R, X}, k::Int) :: AbstractMatrix{R} where {R, X}
    generate_differential(C, k)
end

function set_differential!(C::AbstractComplex{R}, k::Int, A::AbstractMatrix{R}) where {R}
    throw(MethodError(set_differential!, (C, k, A)))
end

function generate_differential(C::AbstractComplex{R, X}, k::Int) :: AbstractMatrix{R} where {R, X}
    Gₖ   = generators(C, k)
    Gₖ₊₁ = generators(C, k + differentialDegree(C))
    generate_matrix(Gₖ, Gₖ₊₁, R) do x 
        differentiate(C, k, x)
    end
end

function generate_matrix(f, Gₖ::Vector{X}, Gₖ₊₁::Vector{X}, ::Type{R}) :: AbstractMatrix{R} where {R, X}
    n = length(Gₖ)
    m = length(Gₖ₊₁)

    gDict = Dict( Gₖ₊₁[i] => i for i in 1 : length(Gₖ₊₁) )
    
    Is = Vector{Int}()
    Js = Vector{Int}()
    Vs = Vector{R}()

    for j in 1 : n 
        x = Gₖ[j]
        ys = f(x)

        for (y, r) in ys
            y ∉ keys(gDict) && continue
            
            i = gDict[y]
            push!(Is, i)
            push!(Js, j)
            push!(Vs, r)
        end
    end
    
    sparse(Is, Js, Vs, m, n)
end

function transform(::AbstractComplex{R}, ::Int) :: Union{AbstractMatrix{R}, Nothing} where {R}
    nothing
end

function set_transform!(C::AbstractComplex{R}, k::Int, T::AbstractMatrix{R}) where {R}
    throw(MethodError(set_differential!, (C, k, T)))
end

function reduce_all!(C::AbstractComplex; with_transform=false)
    for k in hDegRange(C)
        reduce!(C, k; with_transform=with_transform)
    end
end

function reduce!(C::AbstractComplex, k::Int; with_transform=false)
    @debug "reduce C[$k]."

    A = differential(C, k)
    flags = (with_transform, false, false, with_transform) # P and Q⁻¹

    (F, S, p, q) = snf_preprocess(A; flags=flags)
    r = length(F.factors) # diagonal entries of units.

    if r == 0 
        @debug "nothing to reduce."
        return
    end

    k₊₁ = k + differentialDegree(C)
    nₖ   = length(generators(C, k))
    nₖ₊₁ = length(generators(C, k₊₁))

    drop_generators!(C, k, r, q)
    drop_generators!(C, k₊₁, r, p)
    set_differential!(C, k, S)

    if with_transform
        Tₖ = F.T.Q⁻¹[r + 1 : nₖ, :]
        Tₖ₊₁ = F.T.P[r + 1 : nₖ₊₁, :]
        set_transform!(C, k, Tₖ)
        set_transform!(C, k₊₁, Tₖ₊₁)
    end

    @debug """cancelled $r pairs of generators.
      C[$k]: $nₖ -> $(nₖ - r)
      C[$k₊₁]: $nₖ₊₁ -> $(nₖ₊₁ - r)
    """
end

function vectorize(C::AbstractComplex{R, X}, k::Int, z::Dict{X, R}, ::Type{V}) :: V where {R, X, V <: AbstractVector{R}}
    gens = original_generators(C, k)
    n = length(gens)
    gDict = Dict(gens[i] => i for i in 1 : n)

    v = fill(zero(R), n)
    for (x, r) in z
        i = gDict[x]
        v[i] = r
    end

    v = V(v)
    T = transform(C, k)
    if !isnothing(T)
        v = T * v
    end

    v
end