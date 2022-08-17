using Test

@testset verbose = true "Kh" begin
    include("KhAlgebraTest.jl")
    include("KhChainTests.jl")
    include("KhCubeTests.jl")
    include("KhComplexTests.jl")
    include("KhHomologyTests.jl")
end