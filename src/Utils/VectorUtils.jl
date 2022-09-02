function count_occurrence(arr::Vector{X}) :: Dict{X, Int} where {X}
    reduce(arr; init=Dict{X, Int}()) do res, x
        push!(res, x => get!(res, x, 0) + 1)
    end
end

function findfirst_elm(arr::AbstractArray{T}, elm::T) where {T}
    findfirst( x -> x == elm, arr )
end

function delete_where!(f, v::Vector{X}) where X 
    i = 1
    while i <= length(v)
        f(v[i]) ? deleteat!(v, i) : (i += 1)
    end
end
