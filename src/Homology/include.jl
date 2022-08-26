module Homology

include("AbstractComplex.jl")
include("AbstractHomologySummand.jl")
include("AbstractHomology.jl")
include("HomologyComputations.jl")

export AbstractHomology, AbstractHomologySummand, AbstractComplex
export hDegRange, generators, differential, differentialDegree, reduce_all!, reduce!

end