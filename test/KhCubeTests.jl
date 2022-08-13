using Test

@testset "KhCube.jl" begin
    using Khovanov: KhCube, KhCubeVertex, KhAlgGenerator, Kh, vertex, edge, mergeEdge, splitEdge
    using Khovanov: Link, unknot, trefoil

    R = Int

    @testset "raw-vertex" begin
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

    @testset "vertex" begin
        l = Link([0, 0, 1, 1])
        A = KhAlgStructure{R}(Kh)
        cube = KhCube(A, l)
        v0 = vertex(cube, [0])
        v1 = vertex(cube, [1])

        @test v0.state == [0]
        @test length(v0.circles) == 2
        @test length(v0.generators) == 4
        @test v1.state == [1]
        @test length(v1.circles) == 1
        @test length(v1.generators) == 2
    end

    @testset "edge-merge" begin
        l = Link([0, 0, 1, 1])
        A = KhAlgStructure{R}(Kh)
        cube = KhCube(A, l)
        e = edge(cube, [0], [1])

        if isa(e, mergeEdge)
            @test e.from == (1, 2)
            @test e.to == 1
        else
            @test false
        end
    end

    @testset "edge-split" begin
        l = Link([0, 1, 1, 0])
        A = KhAlgStructure{R}(Kh)
        cube = KhCube(A, l)
        e = edge(cube, [0], [1])

        if isa(e, splitEdge)
            @test e.from == 1
            @test e.to == (1, 2)
        else
            @test false
        end
    end

end