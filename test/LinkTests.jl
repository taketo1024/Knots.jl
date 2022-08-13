using Test
using Khovanov: Link, isEmpty, crossingNum, signedCrossingNums, writhe, components, mirror, emptyLink, unknot, trefoil, figure8, hopfLink

function testEmptyLink()
    l = emptyLink
    @test isEmpty(l)
    @test crossingNum(l) == 0
    @test writhe(l) == 0
    @test length(components(l)) == 0
end

function testUnknot()
    l = unknot
    @test !isEmpty(l)
    @test crossingNum(l) == 0
    @test writhe(l) == 0
    @test length(components(l)) == 1
end

function testOrientation() 
    l1 = Link([0, 0, 1, 1])
    @test signedCrossingNums(l1) == (1, 0)

    l2 = Link([0, 1, 1, 0])
    @test signedCrossingNums(l2) == (0, 1)
end

function testMirror()
    l = Link([0, 0, 1, 1])
    m = mirror(l)

    @test signedCrossingNums(l) == (1, 0)
    @test signedCrossingNums(m) == (0, 1)
end

function testTrefoil()
    k = trefoil
    @test length(components(k)) == 1
    @test crossingNum(k) == 3
    @test writhe(k) == -3
end

function testTrefoilMirror()
    k = mirror(trefoil)
    @test length(components(k)) == 1
    @test crossingNum(k) == 3
    @test writhe(k) == 3
end

function testFigureEight()
    k = figure8
    @test length(components(k)) == 1
    @test crossingNum(k) == 4
    @test writhe(k) == 0
end

function testFigureEightMirror()
    k = mirror(figure8)
    @test length(components(k)) == 1
    @test crossingNum(k) == 4
    @test writhe(k) == 0
end

function testHopfLink() 
    l = hopfLink
    @test length(components(l)) == 2
    @test crossingNum(l) == 2
    @test writhe(l) == -2
end

function test2CompUnlink()
    l = Link([1,2,3,4], [3,2,1,4])
    @test length(components(l)) == 2
    @test crossingNum(l) == 2
    @test writhe(l) == 0
end

@testset "Links.jl" begin
    testEmptyLink()
    testUnknot()
    testOrientation()
    testMirror()
    testTrefoil()
    testTrefoilMirror()
    testFigureEight()
    testFigureEightMirror()
    testHopfLink()
    test2CompUnlink()
end