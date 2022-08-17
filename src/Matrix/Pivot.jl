# Implementation based on:
#
# "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
# https://hal.inria.fr/hal-01646133/document
#
# see also: SpaSM (Sparse direct Solver Modulo p)
# https://github.com/cbouilla/spasm


using AbstractAlgebra: RingElement, is_unit
using SparseArrays

mutable struct Pivot{R<:RingElement}
    size::Tuple{Int, Int}
    nonzero::Vector{Vector{Int}}   # row -> [col]
    rowHead::Vector{Int}           # row -> col
    rowWeight::Vector{Int}         # row -> weight
    colWeight::Vector{Int}         # col -> weight
    candidates::Vector{Set{Int}}   # row -> Set(cols)
    hasPivot::Set{Int}             # Set(rows)
    pivots::Vector{Int}            # col -> row
end

Pivot(A::AbstractSparseMatrix{R}) where {R} = begin
    (m, n) = size(A)
    nonzero = map( _ -> Int[], 1:m)
    rowHead   = fill(0, m)
    rowWeight = fill(0, m)
    colWeight = fill(0, n)
    candidates = map( _ -> Set{Int}(), 1:m)
    hasPivot = Set{Int}()
    pivots = fill(0, n)

    for (i, j, a) in zip(findnz(A)...)
        iszero(a) && continue
        push!(nonzero[i], j)
        rowWeight[i] += 1
        colWeight[j] += 1

        if rowHead[i] == 0 || j < rowHead[i]
            rowHead[i] = j
        end
        if isone(a) || isone(-a)
            push!(candidates[i], j)
        end
    end

    Pivot{R}((m, n), nonzero, rowHead, rowWeight, colWeight, candidates, hasPivot, pivots)
end

function findPivots(A::AbstractSparseMatrix{R}) where {R<:RingElement} 
    @debug "find pivot: A = size$(size(A))"
    findFLPivots!(A)
    # findFLColumnPivots!(A)
    # findCycleFreePivots!(A)
end

function result(piv::Pivot) :: Vector{CartesianIndex{2}}
    []
end

# Faugère-Lachartre pivot search
function findFLPivots!(piv::Pivot)
    m = piv.size[1]
    tmp = Dict()
    for i in 1:m
        j = piv.rowHead[i]
        j == 0 && continue
        if !haskey(tmp, j) || piv.rowWeight[i] < piv.rowWeight[tmp[j]]
            tmp[j] = i
        end
    end

    for (j, i) in tmp
        setPivot!(piv, i, j)
    end

    @debug "FL-pivots: $(length(tmp))"
end

function findFLColumnPivots!(piv::Pivot)
end

function findCycleFreePivots!(piv::Pivot)
end

function setPivot!(piv::Pivot, i::Int, j::Int)
    push!(piv.hasPivot, i)
    piv.pivots[j] = i
end