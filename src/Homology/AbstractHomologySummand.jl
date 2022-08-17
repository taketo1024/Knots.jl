using AbstractAlgebra: RingElement, Ring

abstract type AbstractHomologySummand{R<:RingElement, RR<:Ring} end

function baseRing(s::AbstractHomologySummand{R, RR}) :: RR where {R, RR <: Ring}
    throw(MethodError(AbstractHomologySummand, (s,)))
end

function rank(s::AbstractHomologySummand) :: Int
    throw(MethodError(rank, (s,)))
end

function torsions(s::AbstractHomologySummand{R}) :: Vector{R} where {R}
    throw(MethodError(torsions, (s,)))
end

function asString(s::AbstractHomologySummand) :: String
    iszero(s) && return "⋅"
    
    R = baseRing(s)
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