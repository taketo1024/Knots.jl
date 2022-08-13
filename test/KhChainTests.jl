using Test

@testset "KhChain.jl" begin
    using ComputationalHomology
    using Khovanov: KhAlgStructure, Kh, X, I
    using Khovanov: KhChain, dim, length

    R = Int
    C = KhChain{R}

    @testset "empty" begin
        z = C(0, Dict())
        @test dim(z) == 0
        @test length(z) == 0
        println(z)
    end
end