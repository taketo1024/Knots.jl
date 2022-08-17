using Test

@testset verbose = true "Khovanov" begin
    include("Links/LinkTests.jl")
    include("Links/JonesTests.jl")
    include("Matrix/SNFTests.jl")
    include("KhAlgebraTest.jl")
    include("KhChainTests.jl")
    include("KhCubeTests.jl")
    include("KhComplexTests.jl")
    include("KhHomologyTests.jl")
end