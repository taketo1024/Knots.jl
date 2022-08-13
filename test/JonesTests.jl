using Test
using Khovanov: Link, emptyLink, unknot, trefoil, figure8, hopfLink
using Khovanov: jonesPolynomial, P

q = P([0, 1], 0, "q")

function testEmpty()
    l = emptyLink
    p = jonesPolynomial(l)
    @test p == 1
end

function testUnknot()
    l = unknot
    p = jonesPolynomial(l)
    @test p == q + q^(-1)
end

function testTrefoil()
    l = trefoil
    p = jonesPolynomial(l)
    @test p == -q^(-9) + q^(-5) + q^(-3) + q^(-1)
end

function testFigure8()
    l = figure8
    p = jonesPolynomial(l)
    @test p == q^5 + q^(-5)
end

function testHopfLink()
    l = hopfLink
    p = jonesPolynomial(l)
    @test p == 1 + q^(-2) + q^(-4) + q^(-6)
end

@testset "Jones.jl" begin
    testEmpty()
    testUnknot()
    testTrefoil()
    testFigure8()
    testHopfLink()
end