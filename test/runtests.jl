using Test

@testset verbose = true "Khovanov" begin
    include("LinkTests.jl")
    include("JonesTests.jl")
    include("KhAlgebraTest.jl")
    include("KhChainTests.jl")
    include("KhCubeTests.jl")
    include("KhComplexTests.jl")
    include("KhHomologyTests.jl")
end