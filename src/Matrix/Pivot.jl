# Implementation based on:
#
# "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
# https://hal.inria.fr/hal-01646133/document
#
# see also: SpaSM (Sparse direct Solver Modulo p)
# https://github.com/cbouilla/spasm

using AbstractAlgebra: RingElement, is_unit
using SparseArrays
using OrderedCollections
using Permutations: Permutation

export findPivots, pivotPermutations

function findPivots(A::AbstractSparseMatrix{R}) :: Tuple{Vector{Int}, Vector{Int}} where {R<:RingElement} 
    _findPivots(A::AbstractSparseMatrix{R})
end

function pivotPermutations(A::AbstractSparseMatrix{R}) :: Tuple{Permutation, Permutation, Int} where {R<:RingElement} 
    (m, n) = size(A)
    (I, J) = _findPivots(A::AbstractSparseMatrix{R})
    (permutation(I, m), permutation(J, n), length(I))
end

# private

mutable struct Pivot{R<:RingElement}
    size::Tuple{Int, Int}
    entries::Vector{Vector{Int}}   # row -> [col]
    rowHead::Vector{Int}           # row -> col
    rowWeight::Vector{Int}         # row -> weight
    colWeight::Vector{Int}         # col -> weight
    candidates::Vector{Set{Int}}   # row -> Set(cols)
    pivots::OrderedDict{Int, Int} # col -> row
end

Pivot(A::AbstractSparseMatrix{R}) where {R} = begin
    (m, n) = size(A)
    entries = map( _ -> Int[], 1:m)
    rowHead   = fill(0, m)
    rowWeight = fill(0, m)
    colWeight = fill(0, n)
    candidates = map( _ -> Set{Int}(), 1:m)
    pivots = OrderedDict{Int, Int}()

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

    Pivot{R}((m, n), entries, rowHead, rowWeight, colWeight, candidates, pivots)
end

function isCandidate(piv::Pivot, i::Int, j::Int) :: Bool 
    j ∈ piv.candidates[i]
end

function setPivot!(piv::Pivot, i::Int, j::Int)
    @debug "\tpivot: ($i, $j)"
    piv.pivots[j] = i
end

function occupiedCols(piv::Pivot) :: Set{Int}
    reduce(piv.pivots; init=Set()) do res, (_, i)
        union!(res, piv.entries[i])
    end
end

function _findPivots(A::AbstractSparseMatrix) :: Tuple{Vector{Int}, Vector{Int}}
    @debug "find pivots: A = size$(size(A))"

    piv = Pivot(A)
    findFLPivots!(piv)
    findFLColumnPivots!(piv)
    findCycleFreePivots!(piv)

    sortPivots(piv)
end

# Faugère-Lachartre pivot search
function findFLPivots!(piv::Pivot)
    m = piv.size[1]
    pivots = OrderedDict()

    for i in 1:m
        j = piv.rowHead[i]
        !isCandidate(piv, i, j) && continue

        if !haskey(pivots, j) || piv.rowWeight[i] < piv.rowWeight[pivots[j]]
            pivots[j] = i
        end
    end

    for (j, i) in pivots
        setPivot!(piv, i, j)
    end

    @debug "FL-pivots: $(length(piv.pivots))"
end

function findFLColumnPivots!(piv::Pivot)
    before = length(piv.pivots)

    m = piv.size[1]
    occupied = occupiedCols(piv) 
    rows = let 
        remain = Set(1:m)
        for (_, i) in piv.pivots
            delete!(remain, i)
        end
        sort!(collect(remain), by=(i -> piv.rowWeight[i]))
    end

    for i in rows
        isempty(piv.entries[i]) && continue # skip empty rows
        
        candidates = Int[]

        for j in piv.entries[i]
            if j ∉ occupied && isCandidate(piv, i, j)
                push!(candidates, j)
            end
        end

        isempty(candidates) && continue

        sort!(candidates, by=(j -> piv.colWeight[j]))
        j = first(candidates)

        setPivot!(piv, i, j)
        union!(occupied, piv.entries[i])
    end

    @debug "FL-col-pivots: $(length(piv.pivots) - before)"
end

function findCycleFreePivots!(piv::Pivot)
end

function sortPivots(piv::Pivot) :: Tuple{Vector{Int}, Vector{Int}}
    npiv = length(piv.pivots)
    cols = sort!(
        collect(keys(piv.pivots)), 
        by=( j -> piv.rowWeight[ piv.pivots[j] ] )
    )
    dict = Dict( cols[j] => j for j in 1:npiv) # reverse

    data = map(1 : npiv) do c
        targets = Set{Int}() 
        j = cols[c]
        i = piv.pivots[j]
        entries = piv.entries[i]
        for k in entries 
            if k ≠ j && k ∈ keys(piv.pivots)
                push!(targets, dict[k])
            end
        end
        sort!(collect(targets))
    end

    sorted = topsort(data)

    Is = Int[]
    Js = Int[]

    map(sorted) do k 
        j = cols[k]
        i = piv.pivots[j]
        push!(Is, i)
        push!(Js, j)
    end

    (Is, Js)
end

# TODO: move to Utils
function topsort(data::Vector{Vector{Int}}) :: Vector{Int}
    n = length(data)
    n == 0 && return Int[]

    weight = fill(0, n)
    for targets in data
        for i in targets
            weight[i] += 1
        end
    end

    result = Int[]
    queue = filter( i -> weight[i] == 0, 1:n )
    isempty(queue) && error("detect loop")

    while !isempty(queue)
        i = popfirst!(queue)
        push!(result, i)

        for j in data[i]
            weight[j] -= 1
            if weight[j] == 0
                push!(queue, j)
            end
        end
    end

    @assert all( map( w -> w == 0, weight) )
    @assert length(result) == n

    result
end

# TODO: move to Utils
function permutation(indices::Vector{Int}, length::Int) :: Permutation
    result = Int[]
    remain = OrderedSet(1:length)
    for i in indices
        push!(result, i)
        delete!(remain, i)
    end
    for i in remain
        push!(result, i)
    end
    Permutation(result)
end