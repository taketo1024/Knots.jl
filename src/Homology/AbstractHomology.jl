using AbstractAlgebra: RingElement, Ring

abstract type AbstractHomology{R<:RingElement, RR<:Ring} end

# must override 
function complex(H::AbstractHomology{R, RR}) :: AbstractComplex{R, RR} where {R, RR}
    throw(MethodError(complex, (H,)))
end

# must override 
function makeSummand(H::AbstractHomology{R, RR}, rank::Int, torsions::Vector{R}) :: AbstractHomologySummand{R, RR} where {R, RR}
    throw(MethodError(_makeSummand, (H, rank, torsions)))
end

function baseRing(H::AbstractHomology{R, RR}) :: RR where {R, RR <: Ring}
    baseRing(complex(H))
end

function hDegRange(H::AbstractHomology) :: UnitRange{Int}
    hDegRange(complex(H))
end

# can override (e.g. change order of computation, etc)
function computeAll(H::AbstractHomology{R, RR}) :: Dict{Int, AbstractHomologySummand{R, RR}} where {R, RR}
    Dict( k => compute(H, k) for k in hDegRange(H) )
end

# can override (e.g. use different algorithm / cache result etc)
function compute(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    compute_single(H, k)
end

function Base.getindex(H::AbstractHomology{R, RR}, k::Int) :: AbstractHomologySummand{R, RR} where {R, RR}
    compute(H, k)
end

function asString(H::AbstractHomology) :: String
    lines = String[]
    summands = computeAll(H)
    for i in hDegRange(H)
        Hi = asString(summands[i])
        push!(lines, "H[$i] = $Hi")
    end
    join(lines, "\n")
end

function Base.show(io::IO, H::AbstractHomology) 
    print(io, asString(H))
end