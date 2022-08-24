using ..Utils

abstract type AbstractHomologySummand{R} end

function rank(s::AbstractHomologySummand) :: Int
    throw(MethodError(rank, (s,)))
end

function torsions(s::AbstractHomologySummand{R}) :: Vector{R} where {R}
    throw(MethodError(torsions, (s,)))
end

function asString(s::AbstractHomologySummand{R}) :: String where {R}
    iszero(s) && return "⋅"
    
    symbol = Utils.symbol(R)

    r = rank(s)
    r_str = r > 1 ? Utils.superscript(r) : ""

    res = (r > 0) ? ["$symbol$r_str"] : []
    for t in torsions(s)
        push!(res, "$symbol/$t")
    end

    join(res, " ⊕ ")
end

Base.iszero(s::AbstractHomologySummand) = 
    (rank(s) == 0 && isempty(torsions(s)))

Base.show(io::IO, s::AbstractHomologySummand) = 
    print(io, asString(s))