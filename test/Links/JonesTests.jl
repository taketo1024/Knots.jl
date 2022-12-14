using Test

@testset "Jones" begin
    using Polynomials
    using Knots.Links: Link, emptyLink, unknot, trefoil, figure8, hopfLink
    using Knots.Links: JonesPolynomial, jonesPolynomial

    q = variable(JonesPolynomial)

    @testset "empty" begin
        l = emptyLink
        p = jonesPolynomial(l)
        @test p == 1
        x = 1
    end

    @testset "unknot" begin
        l = unknot
        p = jonesPolynomial(l)
        @test p == q + q^(-1)
    end

    @testset "trefoil" begin
        l = trefoil
        p = jonesPolynomial(l)
        @test p == -q^(-9) + q^(-5) + q^(-3) + q^(-1)
    end

    @testset "figure8" begin
        l = figure8
        p = jonesPolynomial(l)
        @test p == q^5 + q^(-5)
    end

    @testset "hopfLink" begin
        l = hopfLink
        p = jonesPolynomial(l)
        @test p == 1 + q^(-2) + q^(-4) + q^(-6)
    end
end