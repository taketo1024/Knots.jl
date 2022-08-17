using Test

@testset "Link" begin
    using Khovanov.Links: Link, isEmpty, crossingNum, signedCrossingNums, writhe, components, mirror, emptyLink, unknot, trefoil, figure8, hopfLink

    @testset "emptyLink" begin
        l = emptyLink
        @test isEmpty(l)
        @test crossingNum(l) == 0
        @test writhe(l) == 0
        @test length(components(l)) == 0
    end

    @testset "unknot" begin
        l = unknot
        @test !isEmpty(l)
        @test crossingNum(l) == 0
        @test writhe(l) == 0
        @test length(components(l)) == 1
    end

    @testset "orientation" begin 
        l1 = Link([0, 0, 1, 1])
        @test signedCrossingNums(l1) == (1, 0)

        l2 = Link([0, 1, 1, 0])
        @test signedCrossingNums(l2) == (0, 1)
    end

    @testset "mirror" begin
        l = Link([0, 0, 1, 1])
        m = mirror(l)

        @test signedCrossingNums(l) == (1, 0)
        @test signedCrossingNums(m) == (0, 1)
    end

    @testset "trefoil" begin
        k = trefoil
        @test length(components(k)) == 1
        @test crossingNum(k) == 3
        @test writhe(k) == -3
    end

    @testset "trefoilMirror" begin
        k = mirror(trefoil)
        @test length(components(k)) == 1
        @test crossingNum(k) == 3
        @test writhe(k) == 3
    end

    @testset "figureEight" begin
        k = figure8
        @test length(components(k)) == 1
        @test crossingNum(k) == 4
        @test writhe(k) == 0
    end

    @testset "figureEightMirror" begin
        k = mirror(figure8)
        @test length(components(k)) == 1
        @test crossingNum(k) == 4
        @test writhe(k) == 0
    end

    @testset "hopfLink" begin 
        l = hopfLink
        @test length(components(l)) == 2
        @test crossingNum(l) == 2
        @test writhe(l) == -2
    end

    @testset "2compUnlink" begin 
        l = Link([1,2,3,4], [3,2,1,4])
        @test length(components(l)) == 2
        @test crossingNum(l) == 2
        @test writhe(l) == 0
    end
end