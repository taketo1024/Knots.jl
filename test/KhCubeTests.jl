using Test

@testset "KhCube.jl" begin
    using Khovanov: KhCubeVertex
    using Khovanov: Link, unknot, trefoil

    @testset "vertex" begin
        l = Link([0, 0, 1, 1])
        v0 = KhCubeVertex(l, [0])
        v1 = KhCubeVertex(l, [1])

        @test v0.state == [0]
        @test length(v0.circles) == 2
        @test length(v0.generators) == 4
        @test v1.state == [1]
        @test length(v1.circles) == 1
        @test length(v1.generators) == 2
    end
end