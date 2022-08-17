using Test

include("Links/include.jl")
include("Matrix/include.jl")

@testset verbose = true "Khovanov" begin
    include("KhAlgebraTest.jl")
    include("KhChainTests.jl")
    include("KhCubeTests.jl")
    include("KhComplexTests.jl")
    include("KhHomologyTests.jl")
end