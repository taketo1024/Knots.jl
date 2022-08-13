using Test

@testset verbose = true "Khovanov" begin
    include("LinkTests.jl")
    include("JonesTests.jl")
end