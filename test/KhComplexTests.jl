using Test

@testset "KhComplex.jl" begin
    using Khovanov: KhComplex, hDegRange, chainGenerators, differential, _chainGenerators
    using Khovanov: KhCube, KhAlgStructure
    using Khovanov: Link, emptyLink, unknot, trefoil, hopfLink, mirror

    R = Int
    A = KhAlgStructure(0, 0)

    @testset "private-chain-generators" begin
        l = Link([0, 0, 1, 2], [1, 3, 3, 2])
        cube = KhCube(A, l)

        G₀ = Dict(_chainGenerators(cube, 0))
        G₁ = Dict(_chainGenerators(cube, 1))
        G₂ = Dict(_chainGenerators(cube, 2))

        @test length( G₀[[0, 0]] ) == 4
        @test length( G₁[[1, 0]] ) == 2
        @test length( G₁[[0, 1]] ) == 8
        @test length( G₂[[1, 1]] ) == 4
    end

    @testset "empty" begin
        l = emptyLink
        C = KhComplex(A, l)

        @test C.degShift == (0, 0)
        @test hDegRange(C) == 0:0
        @test length(chainGenerators(C, 0)) == 1
        @test size(differential(C, 0)) == (0, 1)
    end

    @testset "unknot" begin
        l = unknot
        C = KhComplex(A, l)

        @test C.degShift == (0, 0)
        @test hDegRange(C) == 0:0
        @test length(chainGenerators(C, 0)) == 2
        @test size(differential(C, 0)) == (0, 2)
    end

    @testset "pos-twist" begin
        l  = Link([0, 0, 1, 1])
        C = KhComplex(A, l)

        @test C.degShift == (0, 1)
        @test hDegRange(C) == 0:1
        @test length(chainGenerators(C, 0)) == 4
        @test length(chainGenerators(C, 1)) == 2
        @test size(differential(C, 0)) == (2, 4)
    end

    @testset "neg-twist" begin
        l  = Link([0, 1, 1, 0])
        C = KhComplex(A, l)

        @test C.degShift == (-1, -2)
        @test hDegRange(C) == -1:0
        @test length(chainGenerators(C, -1)) == 2
        @test length(chainGenerators(C,  0)) == 4
        @test size(differential(C, -1)) == (4, 2)
    end

    @testset "pos-neg-twist" begin
        l = Link([0, 0, 1, 2], [1, 3, 3, 2])
        C = KhComplex(A, l)

        @test C.degShift == (-1, -1)
        @test hDegRange(C) == -1:1
        @test length(chainGenerators(C, -1)) == 4
        @test length(chainGenerators(C,  0)) == 10
        @test length(chainGenerators(C,  1)) == 4

        D₋₁ = differential(C, -1)
        D₀  = differential(C,  0)
        
        @test size(D₋₁) == (10, 4)
        @test size(D₀)  == (4, 10)
        @test count(!iszero, D₀ * D₋₁) == 0
    end

    @testset "hopfLink" begin
        l = hopfLink
        C = KhComplex(A, l)

        @test C.degShift == (-2, -4)
        @test hDegRange(C) == -2:0
        @test length(chainGenerators(C, -2)) == 4
        @test length(chainGenerators(C, -1)) == 4
        @test length(chainGenerators(C,  0)) == 4

        D₋₂ = differential(C, -2)
        D₋₁ = differential(C, -1)
        
        @test size(D₋₂) == (4, 4)
        @test size(D₋₁) == (4, 4)
        @test count(!iszero, D₋₁ * D₋₂) == 0
    end

    @testset "trefoil" begin
        l = mirror(trefoil)
        C = KhComplex(A, l)

        @test C.degShift == (0, 3)
        @test hDegRange(C) == 0:3
        @test length(chainGenerators(C, 0)) == 4
        @test length(chainGenerators(C, 1)) == 6
        @test length(chainGenerators(C, 2)) == 12
        @test length(chainGenerators(C, 3)) == 8

        D₀ = differential(C, 0)
        D₁ = differential(C, 1)
        D₂ = differential(C, 2)
        
        @test size(D₀) == (6, 4)
        @test size(D₁) == (12, 6)
        @test size(D₂) == (8, 12)
        @test count(!iszero, D₁ * D₀) == 0
        @test count(!iszero, D₂ * D₁) == 0
    end
end