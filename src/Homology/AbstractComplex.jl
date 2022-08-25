using SparseArrays: sparse

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

function differentiate(C::AbstractComplex{R, X}, k::Int, x::X) :: Vector{Tuple{X, R}} where {R, X}
    throw(MethodError(differentiate, (C, k, x)))
end

function differential(C::AbstractComplex{R, X}, k::Int) :: AbstractMatrix{R} where {R, X}
    Gₖ   = generators(C, k)
    Gₖ₊₁ = generators(C, k + differentialDegree(C))

    n = length(Gₖ)
    m = length(Gₖ₊₁)

    gDict = Dict( Gₖ₊₁[i] => i for i in 1 : length(Gₖ₊₁) )
    
    Is = Vector{Int}()
    Js = Vector{Int}()
    Vs = Vector{R}()

    for j in 1 : n 
        x = Gₖ[j]
        ys = differentiate(C, k, x)

        for (y, r) in ys
            i = gDict[y]
            push!(Is, i)
            push!(Js, j)
            push!(Vs, r)
        end
    end
    
    sparse(Is, Js, Vs, m, n)
end