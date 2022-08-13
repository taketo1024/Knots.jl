using Test

@testset "KhComplex.jl" begin
    using Khovanov: _collectGenerators, _matrix
    using Khovanov: KhCube, KhAlgStructure, Kh
    using Khovanov: Link, unknot, trefoil

    R = Int
    A = KhAlgStructure{R}(Kh)

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