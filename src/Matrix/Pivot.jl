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
    entries::Vector{Vector{Int}}   # row -> [col]
    rowHead::Vector{Int}           # row -> col
    rowWeight::Vector{Int}         # row -> weight
    colWeight::Vector{Int}         # col -> weight
    candidates::Vector{Set{Int}}   # row -> Set(cols)
    pivotRows::Set{Int}            # Set(rows)
    pivotPos::Vector{Int}          # col -> row
end

Pivot(A::AbstractSparseMatrix{R}) where {R} = begin
    (m, n) = size(A)
    entries = map( _ -> Int[], 1:m)
    rowHead   = fill(0, m)
    rowWeight = fill(0, m)
    colWeight = fill(0, n)
    candidates = map( _ -> Set{Int}(), 1:m)
    pivotRows = Set{Int}()
    pivotPos = fill(0, n)

    for (i, j, a) in zip(findnz(A)...)
        iszero(a) && continue
        push!(entries[i], j)
        rowWeight[i] += 1
        colWeight[j] += 1

        if rowHead[i] == 0 || j < rowHead[i]
            rowHead[i] = j
        end
        if isone(a) || isone(-a)
            push!(candidates[i], j)
        end
    end

    Pivot{R}((m, n), entries, rowHead, rowWeight, colWeight, candidates, pivotRows, pivotPos)
end

function findPivots!(A::AbstractSparseMatrix{R}) where {R<:RingElement} 
    @debug "find pivot: A = size$(size(A))"
    findFLPivots!(A)
    # findFLColumnpivotPos!(A)
    # findCycleFreepivotPos!(A)
end

function result(piv::Pivot) :: Vector{CartesianIndex{2}}
    []
end

# Faugère-Lachartre pivot search
function findFLPivots!(piv::Pivot)
    m = piv.size[1]
    pivots = Dict()

    for i in 1:m
        j = piv.rowHead[i]
        j ∉ piv.candidates[i] && continue

        if !haskey(pivots, j) || piv.rowWeight[i] < piv.rowWeight[pivots[j]]
            pivots[j] = i
        end
    end

    for (j, i) in pivots
        setPivot!(piv, i, j)
    end

    @debug "FL-pivots: $(length(piv.pivotRows))"
end

function findFLColumnPivots!(piv::Pivot)
    before = length(piv.pivotRows)

    m = piv.size[1]
    occupied = reduce(piv.pivotRows; init=Set()) do res, i
        union!(res, piv.entries[i])
    end

    rows = collect( setdiff(Set(1:m), piv.pivotRows) )
    sort!(rows, by=(i -> piv.rowWeight[i]))

    for i in rows
        isempty(piv.entries[i]) && continue # skip empty rows
        
        candidates = Int[]

        for j in piv.entries[i]
            if j ∉ occupied && j ∈ piv.candidates[i] 
                push!(candidates, j)
            end
        end

        isempty(candidates) && continue

        sort!(candidates, by=(j -> colWeight[j]))
        j = first(candidates)

        setPivot!(piv, i, j)
        union!(occupied, piv.entries[i])
    end

    @debug "FL-col-pivots: $(length(piv.pivotRows) - before)"
end

function findCycleFreePivots!(piv::Pivot)
end

function setPivot!(piv::Pivot, i::Int, j::Int)
    push!(piv.pivotRows, i)
    piv.pivotPos[j] = i
end