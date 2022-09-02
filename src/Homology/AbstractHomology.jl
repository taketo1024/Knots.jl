abstract type AbstractHomology{R} end

# must override 
function complex(H::AbstractHomology{R}) :: AbstractComplex{R} where {R}
    throw(MethodError(complex, (H,)))
end

# must override 
function makeSummand(H::AbstractHomology{R}, k::Int, rank::Int, torsions::Vector{R}) :: AbstractHomologySummand{R} where {R}
    throw(MethodError(_makeSummand, (H, k, rank, torsions)))
end

function hDegRange(H::AbstractHomology) :: UnitRange{Int}
    hDegRange(complex(H))
end

# can override (e.g. change order of computation, etc)
function computeAll(H::AbstractHomology{R}) :: Dict{Int, AbstractHomologySummand{R}} where {R}
    Dict( k => compute(H, k) for k in hDegRange(H) )
end

# can override (e.g. use different algorithm / cache result etc)
function compute(H::AbstractHomology{R}, k::Int) :: AbstractHomologySummand{R} where {R}
    compute_single(H, k)
end

function Base.getindex(H::AbstractHomology{R}, k::Int) :: AbstractHomologySummand{R} where {R}
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