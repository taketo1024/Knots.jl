using Test

@testset verbose=true "KhComplex.jl" begin
    using Khovanov: KhComplex, hDegRange, _collectGenerators, _matrix
    using Khovanov: KhCube, KhAlgStructure, Kh
    using Khovanov: Link, emptyLink, unknot, trefoil

    R = Int
    A = KhAlgStructure{R}(Kh)

    @testset "empty" begin
        l = emptyLink
        C = KhComplex(A, l)

        @test C.degShift == (0, 0)
        @test hDegRange(C) == 0:0
    end

    @testset "unknot" begin
        l = unknot
        C = KhComplex(A, l)

        @test C.degShift == (0, 0)
        @test hDegRange(C) == 0:0
    end

    @testset "pos-cross" begin
        l = l = Link([0, 0, 1, 1])
        C = KhComplex(A, l)

        @test C.degShift == (0, 1)
        @test hDegRange(C) == 0:1
    end

    @testset "neg-cross" begin
        l = l = Link([0, 1, 1, 0])
        C = KhComplex(A, l)

        @test C.degShift == (-1, -2)
        @test hDegRange(C) == -1:0
    end

    @testset "collect-generators" begin
        l = Link([0, 0, 1, 3], [2, 3, 1, 2])
        cube = KhCube(A, l)

        G₀ = Dict(_collectGenerators(cube, 0))
        G₁ = Dict(_collectGenerators(cube, 1))
        G₂ = Dict(_collectGenerators(cube, 2))

        @test length( G₀[[0, 0]] ) == 4
        @test length( G₁[[1, 0]] ) == 2
        @test length( G₁[[0, 1]] ) == 8
        @test length( G₂[[1, 1]] ) == 4
    end

    @testset "matrix" begin
        l = Link([0, 0, 1, 3], [2, 3, 1, 2])
        cube = KhCube(A, l)

        A = _matrix(cube, 0)
        B = _matrix(cube, 1)
        
        @test size(A) == (10, 4)
        @test size(B) == (4, 10)
        @test count(!iszero, B * A) == 0
    end
end