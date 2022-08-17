module Homology

include("AbstractComplex.jl")
include("AbstractHomologySummand.jl")
include("AbstractHomology.jl")
include("HomologyComputations.jl")

export AbstractHomology, AbstractHomologySummand, AbstractComplex
export baseRing, hDegRange, generators, differential, differentialDegree

end