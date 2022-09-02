using Test

@testset "KhHomology" begin
    using Knots.Khovanov: KhAlgStructure, KhHomology
    using Knots.Links

    @testset "emptyLink" begin
        l = emptyLink
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

        @test iszero(H[-1])

        @test H[0].rank == 1
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end
    
    @testset "unknot" begin
        l = unknot
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "unknot-RM1" begin
        l = Link([1,2,2,1])
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "unknot-RM2" begin
        l = Link([1,4,2,1], [2,4,3,3])
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

        @test iszero(H[-1])

        @test H[0].rank == 2
        @test H[0].torsions == []
        
        @test iszero(H[1])
    end

    @testset "trefoil" begin
        l = trefoil
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

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
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

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
        A = KhAlgStructure("Z-Kh")
        H = KhHomology(l, A)

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

    @testset "Q-trefoil" begin
        l = trefoil
        A = KhAlgStructure("Q-Kh")
        H = KhHomology(l, A)

        @test H[-3].rank == 1
        @test H[-3].torsions == []

        @test H[-2].rank == 1
        @test H[-2].torsions == []

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F2-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F2-Kh")
        H = KhHomology(l, A)

        @test H[-3].rank == 2
        @test H[-3].torsions == []

        @test H[-2].rank == 2
        @test H[-2].torsions == []

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F3-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F3-Kh")
        H = KhHomology(l, A)

        @test H[-3].rank == 1
        @test H[-3].torsions == []

        @test H[-2].rank == 1
        @test H[-2].torsions == []

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "Q[h]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("Q[h]-bigraded")
        H = KhHomology(l, A)

        h = A.h

        @test H[-3].rank == 0
        @test H[-3].torsions == []

        @test H[-2].rank == 0
        @test H[-2].torsions == [h^2]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F2[h]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F2[h]-bigraded")
        H = KhHomology(l, A)

        h = A.h

        @test H[-3].rank == 0
        @test H[-3].torsions == []

        @test H[-2].rank == 0
        @test H[-2].torsions == [h, h]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F3[h]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F3[h]-bigraded")
        H = KhHomology(l, A)

        h = A.h

        @test H[-3].rank == 0
        @test H[-3].torsions == []

        @test H[-2].rank == 0
        @test H[-2].torsions == [h^2]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "Q[t]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("Q[t]-bigraded")
        H = KhHomology(l, A)

        t = A.t

        @test H[-3].rank == 0
        @test H[-3].torsions == []

        @test H[-2].rank == 0
        @test H[-2].torsions == [t]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F2[t]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F2[t]-bigraded")
        H = KhHomology(l, A)

        t = A.t

        @test H[-3].rank == 2
        @test H[-3].torsions == []

        @test H[-2].rank == 2
        @test H[-2].torsions == []

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end

    @testset "F3[h]-trefoil" begin
        l = trefoil
        A = KhAlgStructure("F3[t]-bigraded")
        H = KhHomology(l, A)

        t = A.t

        @test H[-3].rank == 0
        @test H[-3].torsions == []

        @test H[-2].rank == 0
        @test H[-2].torsions == [t]

        @test H[-1].rank == 0
        @test H[-1].torsions == []

        @test H[0].rank == 2
        @test H[0].torsions == []
    end
end
