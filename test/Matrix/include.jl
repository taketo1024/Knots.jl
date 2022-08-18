using Test

@testset verbose=true "Matrix" begin
    include("PivotTests.jl")
    include("SNFTests.jl")
end