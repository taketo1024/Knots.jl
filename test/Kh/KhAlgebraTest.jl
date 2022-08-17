using Test

@testset "KhAlgebra" begin
    using Khovanov.Kh: KhAlgStructure, product, coproduct, X, I
    
    @testset "type:Kh" begin
        A = KhAlgStructure(0, 0)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == []
        @test Δ(I) == [(I, X, 1), (X, I, 1)]
        @test Δ(X) == [(X, X, 1)]
    end

    @testset "type:BN" begin
        A = KhAlgStructure(1, 0)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == [(X, 1)]
        @test Δ(I) == [(I, X, 1), (X, I, 1), (I, I, -1)]
        @test Δ(X) == [(X, X, 1)]
    end

    @testset "type:Lee" begin
        A = KhAlgStructure(0, 1)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == [(I, 1)]
        @test Δ(I) == [(I, X, 1), (X, I, 1)]
        @test Δ(X) == [(X, X, 1), (I, I, 1)]
    end
end