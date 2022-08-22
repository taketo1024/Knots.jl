using Test

@testset "Env" begin
    using Knots.Env: set_base_ring, current_base_ring, check_base_ring, get_base_ring
    using AbstractAlgebra

    set_base_ring(nothing)

    @testset "basering-Z" begin
        @assert isnothing(current_base_ring())
        @test get_base_ring(Int) == parent(1)
    end

    @testset "basering-Q" begin
        @assert isnothing(current_base_ring())
        @test get_base_ring(Rational) == parent(1//1)
    end

    @testset "basering-ZZ" begin
        R = elem_type(ZZ)
        @assert isnothing(current_base_ring())
        @test get_base_ring(R) == ZZ
    end

    @testset "basering-QQ" begin
        R = elem_type(QQ)
        @assert isnothing(current_base_ring())
        @test get_base_ring(R) == QQ
    end

    @testset "basering-F2" begin
        F = GF(2)
        R = elem_type(F)
        @assert isnothing(current_base_ring())

        set_base_ring(F)
        @test get_base_ring(R) == F
        set_base_ring(nothing)
    end

    @testset "basering-Z[x]" begin
        RR, _ = PolynomialRing(ZZ, :x)
        R = elem_type(RR)
        @assert isnothing(current_base_ring())

        set_base_ring(RR)
        @test get_base_ring(R) == RR
        set_base_ring(nothing)
    end

    @testset "basering-not-set" begin
        F = GF(2)
        R = elem_type(F)
        @assert isnothing(current_base_ring())

        @test check_base_ring(R, warn=false) == false
    end

    @testset "basering-incompatible" begin
        F = GF(2)
        R = elem_type(F)
        @assert isnothing(current_base_ring())

        set_base_ring(parent(1))
        @test check_base_ring(R, warn=false) == false
        set_base_ring(nothing)
    end

    set_base_ring(nothing)
end