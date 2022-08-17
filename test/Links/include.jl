using Test

@testset verbose=true "Links" begin
    include("LinkTests.jl")
    include("JonesTests.jl")
end