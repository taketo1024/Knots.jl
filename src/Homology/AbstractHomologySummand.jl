using ..Extensions: symbol
using ..Utils: superscript, count_occurrence

abstract type AbstractHomologySummand{R} end

function h_degree(s::AbstractHomologySummand) :: Int
    throw(MethodError(h_degree, (s,)))
end

function parent(s::AbstractHomologySummand{R}) :: AbstractHomology{R} where {R}
    throw(MethodError(parent, (s,)))
end

function rank(s::AbstractHomologySummand) :: Int
    throw(MethodError(rank, (s,)))
end

function torsions(s::AbstractHomologySummand{R}) :: Vector{R} where {R}
    throw(MethodError(torsions, (s,)))
end

function transform(s::AbstractHomologySummand{R}) :: AbstractMatrix{R} where {R}
    throw(MethodError(transform, (s,)))
end

function vectorize(s::AbstractHomologySummand{R}, z::Dict{X, R}, v_type::Type{V}) :: V where {R, X, V <: AbstractVector{R}}
    C = complex(parent(s))
    k = h_degree(s)
    v = vectorize(C, k, z, v_type)
    T = transform(s)
    T * v
end

function asString(s::AbstractHomologySummand{R}) :: String where {R}
    iszero(s) && return "⋅"
    
    res = String[]
    R_str = symbol(R)

    r = rank(s)
    if r > 0
        r_str = r > 1 ? R_str * superscript(r) : R_str
        push!(res, r_str)
    end

    tor_count = count_occurrence(torsions(s))
    for (t, r) in tor_count
        t_str = "$R_str/$t"
        t_str = r > 1 ? "($t_str)" * superscript(r) : t_str
        push!(res, t_str)
    end

    join(res, " ⊕ ")
end

Base.iszero(s::AbstractHomologySummand) = 
    (rank(s) == 0 && isempty(torsions(s)))

Base.show(io::IO, s::AbstractHomologySummand) = 
    print(io, asString(s))