using Test

@testset "KhAlgebra.jl" begin
    using Khovanov: KhAlgStructure, product, coproduct, X, I, Kh, BN, Lee
    
    @testset "type:Kh" begin
        A = KhAlgStructure{Int}(Kh)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == []
        @test Δ(I) == [(I, X, 1), (X, I, 1)]
        @test Δ(X) == [(X, X, 1)]
    end

    @testset "type:BN" begin
        A = KhAlgStructure{Int}(BN)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == [(X, 1)]
        @test Δ(I) == [(I, X, 1), (X, I, 1), (I, I, -1)]
        @test Δ(X) == [(X, X, 1)]
    end

    @testset "type:Lee" begin
        A = KhAlgStructure{Int}(Lee)
        m = product(A)
        Δ = coproduct(A)

        @test m(I, I) == [(I, 1)]
        @test m(I, X) == m(I, X) == [(X, 1)]
        @test m(X, X) == [(I, 1)]
        @test Δ(I) == [(I, X, 1), (X, I, 1)]
        @test Δ(X) == [(X, X, 1), (I, I, 1)]
    end
end