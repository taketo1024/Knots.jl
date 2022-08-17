using Test

@testset "KhHomology.jl" begin
    using Khovanov: KhAlgStructure, KhHomology
    using Khovanov.Links: Link, emptyLink, unknot, trefoil, figure8, hopfLink, mirror

    R = Int
    A = KhAlgStructure(0, 0)

    @testset "emptyLink" begin
        l = emptyLink
        H = KhHomology(A, l)

        @test iszero(H[-1])

        @test H[0].rank == 1
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end
    
    @testset "unknot" begin
        l = unknot
        H = KhHomology(A, l)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "unknot-RM1" begin
        l = Link([1,2,2,1])
        H = KhHomology(A, l)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "unknot-RM2" begin
        l = Link([1,4,2,1], [2,4,3,3])
        H = KhHomology(A, l)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "trefoil" begin
        l = trefoil
        H = KhHomology(A, l)

        @test H[-3].rank == 1
        @test H[-3].torsions == []

        @test H[-2].rank == 1
        @test H[-2].torsions == [2]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "trefoil-mirror" begin
        l = mirror(trefoil)
        H = KhHomology(A, l)

        @test H[0].rank == 2
        @test H[0].torsions == []

        @test H[1].rank == 0
        @test H[1].torsions == []
        
        @test H[2].rank == 1
        @test H[2].torsions == []

        @test H[3].rank == 1
        @test H[3].torsions == [2]
    end

    @testset "figure8" begin
        l = figure8
        H = KhHomology(A, l)

        @test H[-2].rank == 1
        @test H[-2].torsions == []

        @test H[-1].rank == 1
        @test H[-1].torsions == [2]

        @test H[0].rank == 2
        @test H[0].torsions == []

        @test H[1].rank == 1
        @test H[1].torsions == []

        @test H[2].rank == 1
        @test H[2].torsions == [2]
    end
end
