using SparseArrays: sparse
using ..Matrix: snf_preprocess

abstract type AbstractComplex{R, X} end

function hDegRange(C::AbstractComplex) :: UnitRange{Int}
    throw(MethodError(hDegRange, (C,)))
end

function generators(C::AbstractComplex{R, X}, k::Int) :: Vector{X} where {R, X}
    throw(MethodError(generators, (C, k)))
end

function differentialDegree(C::AbstractComplex) :: Int
    +1
end

function differentiate(C::AbstractComplex{R, X}, k::Int, x::X) :: Vector{Pair{X, R}} where {R, X}
    throw(MethodError(differentiate, (C, k, x)))
end

function differential(C::AbstractComplex{R, X}, k::Int) :: AbstractMatrix{R} where {R, X}
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

function set_generators(C::AbstractComplex{R, X}, ::Int, ::Vector{X}) where {R, X}
    @warn "`set_generators` not implemented in $(typeof(C))."
end

function set_differential(C::AbstractComplex{R}, ::Int, ::AbstractMatrix{R}) where {R}
    @warn "`set_differential` not implemented in $(typeof(C))."
end

function reduce_all!(C::AbstractComplex)
    for k in hDegRange(C)
        reduce!(C, k)
    end
end

function reduce!(C::AbstractComplex, k::Int)
    @debug "reduce C[$k]."

    A = differential(C, k)

    (F, S, p, q) = snf_preprocess(A; flags=(false, false, false, false))
    r = length(F.factors) # diagonal entries of units.

    if r == 0 
        @debug "nothing to reduce."
        return
    end

    k₊₁ = k + differentialDegree(C)

    Gₖ   = generators(C, k)
    Gₖ₊₁ = generators(C, k₊₁)

    for j in reverse( sort([ q(j) for j in 1:r ]) )
        deleteat!(Gₖ, j)
    end

    for i in reverse( sort([ p(i) for i in 1:r ]) )
        deleteat!(Gₖ₊₁, i)
    end

    set_generators(C, k, Gₖ)
    set_generators(C, k₊₁, Gₖ₊₁)
    set_differential(C, k, S)

    @debug """cancelled $r pairs of generators.
    C[$k]: $(length(Gₖ) + r) -> $(length(Gₖ))
    C[$k₊₁]: $(length(Gₖ₊₁) + r) -> $(length(Gₖ₊₁))
    """
end