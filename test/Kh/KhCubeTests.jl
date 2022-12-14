using Test

@testset "KhCube" begin
    using Knots.Khovanov: KhCube, KhCubeVertex, vertex, nextVertices, edge, edgeSign, mergeEdge, splitEdge, differentiate
    using Knots.Khovanov: KhAlgStructure, KhAlgGenerator
    using Knots.Links

    R = Int
    A = KhAlgStructure(0, 0)

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

    @testset "empty" begin
        l = emptyLink 
        cube = KhCube(l, A)
        v = vertex(cube, Int[])

        @test length(v.generators) == 1
    end

    @testset "unknot" begin
        l = unknot 
        cube = KhCube(l, A)
        v = vertex(cube, Int[])

        @test length(v.generators) == 2
    end

    @testset "one-crossing" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(l, A)
        v0 = vertex(cube, [0])
        v1 = vertex(cube, [1])

        @test v0.state == [0]
        @test length(v0.circles) == 2
        @test length(v0.generators) == 4
        @test v1.state == [1]
        @test length(v1.circles) == 1
        @test length(v1.generators) == 2
    end

    @testset "nextVertices" begin
        l = trefoil
        cube = KhCube(l, A)
        @test nextVertices(cube, [0, 0, 0]) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        @test nextVertices(cube, [1, 0, 0]) == [[1, 1, 0], [1, 0, 1]]
        @test nextVertices(cube, [1, 1, 0]) == [[1, 1, 1]]
        @test nextVertices(cube, [1, 1, 1]) == []
    end

    @testset "edge-merge" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(l, A)
        e = edge(cube, [0], [1])
        @test e == mergeEdge([0], [1], +1, ((1, 2), 1))
    end

    @testset "edge-split" begin
        l = Link([0, 1, 1, 0])
        cube = KhCube(l, A)
        e = edge(cube, [0], [1])
        @test e == splitEdge([0], [1], +1, (1, (1, 2)))
    end

    @testset "2-crossing-edge-sign" begin
        l = Link([0, 0, 1, 2], [1, 3, 3, 2])
        cube = KhCube(l, A)
        e1 = edgeSign(cube, [0, 0], [1, 0])
        e2 = edgeSign(cube, [0, 0], [0, 1])
        e3 = edgeSign(cube, [1, 0], [1, 1])
        e4 = edgeSign(cube, [0, 1], [1, 1])

        @test e1 == +1
        @test e2 == +1
        @test e3 == -1
        @test e4 == +1
    end

    @testset "2-crossing-edge" begin
        l = Link([0, 0, 1, 2], [1, 3, 3, 2])
        cube = KhCube(l, A)
        e1 = edge(cube, [0, 0], [1, 0])
        e2 = edge(cube, [0, 0], [0, 1])
        e3 = edge(cube, [1, 0], [1, 1])
        e4 = edge(cube, [0, 1], [1, 1])

        @test e1 == mergeEdge([0, 0], [1, 0], +1, ((1, 2), 1))
        @test e2 == splitEdge([0, 0], [0, 1], +1, (2, (2, 3)))
        @test e3 == splitEdge([1, 0], [1, 1], -1, (1, (1, 2)))
        @test e4 == mergeEdge([0, 1], [1, 1], +1, ((1, 2), 1))
    end

    @testset "differential-merge" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(l, A)

        u = [0]
        v = [1]
        x = vertex(cube, u).generators
        y = vertex(cube, v).generators

        @test differentiate(cube, x[1]) == []
        @test differentiate(cube, x[2]) == [(y[1] => 1)]
        @test differentiate(cube, x[3]) == [(y[1] => 1)]
        @test differentiate(cube, x[4]) == [(y[2] => 1)]
    end

    @testset "differential-split" begin
        l = Link([0, 1, 1, 0])
        cube = KhCube(l, A)

        u = [0]
        v = [1]
        x = vertex(cube, u).generators
        y = vertex(cube, v).generators

        @test differentiate(cube, x[1]) == [(y[1] => 1)]
        @test differentiate(cube, x[2]) == [(y[2] => 1), (y[3] => 1)]
    end

    @testset "cache" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(l, A)

        @test length(cube.vertices) == 0

        edge(cube, [0], [1])

        @test length(cube.vertices) == 2
    end
end