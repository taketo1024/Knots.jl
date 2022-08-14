using Test

@testset "KhCube.jl" begin
    using Khovanov: KhCube, KhCubeVertex, vertex, nextVertices, edge, edgeSign, edgeMap, mergeEdge, splitEdge
    using Khovanov: KhAlgStructure, KhAlgGenerator, Kh
    using Khovanov: Link, unknot, emptyLink, trefoil

    R = Int
    A = KhAlgStructure{R}(Kh)

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
        cube = KhCube(A, l)
        v = vertex(cube, Int[])

        @test length(v.generators) == 1
    end

    @testset "unknot" begin
        l = unknot 
        cube = KhCube(A, l)
        v = vertex(cube, Int[])

        @test length(v.generators) == 2
    end

    @testset "one-crossing" begin
        l = Link([0, 0, 1, 1])
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

    @testset "nextVertices" begin
        l = trefoil
        cube = KhCube(A, l)
        @test nextVertices(cube, [0, 0, 0]) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        @test nextVertices(cube, [1, 0, 0]) == [[1, 1, 0], [1, 0, 1]]
        @test nextVertices(cube, [1, 1, 0]) == [[1, 1, 1]]
        @test nextVertices(cube, [1, 1, 1]) == []
    end

    @testset "edge-merge" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(A, l)
        e = edge(cube, [0], [1])
        @test e == mergeEdge(+1, (1, 2), 1)
    end

    @testset "edge-split" begin
        l = Link([0, 1, 1, 0])
        cube = KhCube(A, l)
        e = edge(cube, [0], [1])
        @test e == splitEdge(+1, 1, (1, 2))
    end

    @testset "2-crossing-edge-sign" begin
        l = Link([0, 0, 1, 2], [1, 3, 3, 2])
        cube = KhCube(A, l)
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
        cube = KhCube(A, l)
        e1 = edge(cube, [0, 0], [1, 0])
        e2 = edge(cube, [0, 0], [0, 1])
        e3 = edge(cube, [1, 0], [1, 1])
        e4 = edge(cube, [0, 1], [1, 1])

        @test e1 == mergeEdge(+1, (1, 2), 1)
        @test e2 == splitEdge(+1, 2, (2, 3))
        @test e3 == splitEdge(-1, 1, (1, 2))
        @test e4 == mergeEdge(+1, (1, 2), 1)
    end

    @testset "edgemap-merge" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(A, l)

        u = [0]
        v = [1]
        x = vertex(cube, u).generators
        y = vertex(cube, v).generators

        @test edgeMap(cube, u, v, x[1]) == []
        @test edgeMap(cube, u, v, x[2]) == [(y[1], 1)]
        @test edgeMap(cube, u, v, x[3]) == [(y[1], 1)]
        @test edgeMap(cube, u, v, x[4]) == [(y[2], 1)]
    end

    @testset "edgemap-split" begin
        l = Link([0, 1, 1, 0])
        cube = KhCube(A, l)

        u = [0]
        v = [1]
        x = vertex(cube, u).generators
        y = vertex(cube, v).generators

        @test edgeMap(cube, u, v, x[1]) == [(y[1], 1)]
        @test edgeMap(cube, u, v, x[2]) == [(y[2], 1), (y[3], 1)]
    end

    @testset "cache" begin
        l = Link([0, 0, 1, 1])
        cube = KhCube(A, l)

        @test length(cube._vertexCache) == 0
        @test length(cube._edgeCache) == 0

        edge(cube, [0], [1])

        @test length(cube._vertexCache) == 2
        @test length(cube._edgeCache) == 1
    end
end