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
using Base.Threads
using ..Utils: ReadWriteLock, read_lock, write_lock

export Pivot, pivot, coordinates, permutations

mutable struct Pivot{R<:RingElement}
    size::Tuple{Int, Int}
    entries::Vector{Vector{Int}}   # row -> [col]
    rowHead::Vector{Int}           # row -> col
    rowWeight::Vector{Int}         # row -> weight
    colWeight::Vector{Int}         # col -> weight
    candidates::Vector{Set{Int}}   # row -> Set(cols)
    pivots::OrderedDict{Int, Int}  # col -> row
end

Pivot(A::SparseMatrix{R}) where {R<:RingElement} = begin
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

function pivot(A::SparseMatrix{R}) :: Pivot{R} where {R}
    piv = Pivot(A)
    findPivots!(piv)
end

function npivots(piv::Pivot) :: Int
    length(piv.pivots)
end

function coordinates(piv::Pivot) :: Vector{CartesianIndex{2}}
    map(collect(piv.pivots)) do (j, i)
        CartesianIndex(i, j)
    end
end

function permutations(piv::Pivot) :: Tuple{Permutation, Permutation, Int}
    (m, n) = piv.size
    I = collect(values(piv.pivots))
    J = collect(keys(piv.pivots))
    (permutation(I, m), permutation(J, n), length(I))
end

# private

function isCandidate(piv::Pivot, i::Int, j::Int) :: Bool 
    j ∈ piv.candidates[i]
end

function setPivot!(piv::Pivot, i::Int, j::Int)
    piv.pivots[j] = i
end

function remainingRows(piv::Pivot) :: Vector{Int}
    m = piv.size[1]
    remain = Set(1:m)
    for (_, i) in piv.pivots
        delete!(remain, i)
    end
    sort!(collect(remain), by=(i -> piv.rowWeight[i]))
end

function occupiedCols(piv::Pivot) :: Set{Int}
    reduce(piv.pivots; init=Set()) do res, (_, i)
        union!(res, piv.entries[i])
    end
end

function findPivots!(piv::Pivot)
    @debug "find pivots: A = size$(piv.size)"

    findFLPivots!(piv)
    findFLColumnPivots!(piv)
    findCycleFreePivots!(piv)

    sortPivots!(piv)
    piv
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

    npiv = length(piv.pivots)
    @debug "found FL-pivots: $npiv."
end

function findFLColumnPivots!(piv::Pivot)
    before = length(piv.pivots)
    (before == min(piv.size...)) && return

    remain = remainingRows(piv)
    occupied = occupiedCols(piv) 

    for i in remain
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

    npiv = length(piv.pivots)
    @debug "found FL-col-pivots: $(npiv - before), total: $npiv."
end

function findCycleFreePivots!(piv::Pivot) :: Union{Int, Nothing}
    before = length(piv.pivots)
    (before == min(piv.size...)) && return

    remain = remainingRows(piv)

    @debug "start findCycleFreePivots, remain: $(length(remain))"

    if Threads.nthreads() > 0
        findCycleFreePivots_p!(piv, remain)
    else 
        findCycleFreePivots_s!(piv, remain)
    end

    npiv = length(piv.pivots)
    @debug "cycle-free-pivots: $(npiv - before), total: $npiv."
end

function findCycleFreePivots_s!(piv::Pivot, remain::Vector{Int}) :: Union{Int, Nothing}
    while !isempty(remain)
        i = popfirst!(remain)
        j = findCycleFreePivot(piv, i)
        isnothing(j) && continue
        setPivot!(piv, i, j)
    end
end

function findCycleFreePivots_p!(piv::Pivot, remain::Vector{Int}) :: Union{Int, Nothing}
    lock = ReadWriteLock()

    @threads for i in remain
        while true 
            (npiv, j) = read_lock(lock) do 
                npiv = length(piv.pivots)
                j = findCycleFreePivot(piv, i)
                (npiv, j)
            end

            isnothing(j) && break

            success = write_lock(lock) do 
                if npiv == length(piv.pivots)
                    setPivot!(piv, i, j)
                    true
                else
                    false
                end
            end

            success && break
        end
    end
end

function findCycleFreePivot(piv::Pivot, i::Int) :: Union{Int, Nothing}
    #          j1          j2
    #     i [  o   #   #   #      # ]    *: pivot,
    #          |           :             #: candidate,
    #          V           : rmv         o: queued,
    #    i2 [  * --> o --> x        ]    x: entry
    #                |            
    #                V            
    #       [        * ------> o    ]
    #                          |  
    candidates = Set{Int}()
    queue = Int[]
    added = Set{Int}()

    for j in piv.entries[i]
        if j ∈ keys(piv.pivots)
            push!(queue, j)
            push!(added, j)
        elseif isCandidate(piv, i, j)
            push!(candidates, j)
        end
    end

    while !isempty(queue) && !isempty(candidates)
        j1 = popfirst!(queue)
        i2 = piv.pivots[j1]

        for j2 in piv.entries[i2]
            if j2 ∈ keys(piv.pivots) && j2 ∉ added
                push!(queue, j2)
                push!(added, j2)
            elseif j2 ∈ candidates
                delete!(candidates, j2)
                isempty(candidates) && break
            end
        end
    end

    if isempty(candidates)
        nothing
    else
        argmin(j -> piv.colWeight[j], candidates)
    end
end

function sortPivots!(piv::Pivot)
    (tree, cols) = makeTree(piv)
    sorted = topsort(tree)
    piv.pivots = OrderedDict( 
        cols[idx] => piv.pivots[ cols[idx] ] 
        for idx in sorted 
    )
end

function makeTree(piv::Pivot) :: Tuple{ Vector{Set{Int}}, Vector{Int} }
    npiv = length(piv.pivots)
    cols = sort!(
        collect( keys(piv.pivots) ), 
        by=( j -> piv.rowWeight[ piv.pivots[j] ] )
    )
    reindex = Dict( cols[idx] => idx for idx in 1:npiv)
    tree = map(1:npiv) do idx 
        j = cols[idx]
        i = piv.pivots[j]

        targets = Set{Int}() 
        for k in piv.entries[i]
            if k ≠ j && k ∈ keys(piv.pivots)
                push!(targets, reindex[k])
            end
        end

        targets
    end
    
    (tree, cols)
end

# TODO: move to Utils
function topsort(data::Vector{Set{Int}}) :: Vector{Int}
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