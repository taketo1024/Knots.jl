using Test
using GaloisFields

@testset "SNF" begin
    using SparseArrays: sparse, dropzeros, blockdiag, spdiagm
    using Polynomials
    using GaloisFields
    using Knots.Matrix: SNF, snf, print_matrix
    using Knots.Khovanov: KhAlgStructure

    @testset "snf-no-preprocess-simple" begin 
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 3 1 0 1 0 2 0;
            0 0 1 1 0 0 0 5 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 2 0 1 1 0 1 1
        ])
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        @test d == [1, 1, 1, 1, 1, 1]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-no-preprocess-random" begin 
        density = 0.1
        (m, n) = (30, 40)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-no-preprocess-Z" begin 
        # d[3] of CKh(R-trefoil)
        A = sparse(
            [1, 2, 3, 5, 4, 6, 1, 2, 5, 3, 4, 7, 1, 2, 3, 5, 6, 7], 
            [1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12], 
            [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1], 
            8, 12
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        @test d == [1, 1, 1, 1, 1, 1, 2]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-no-preprocess-Q[h]" begin
        R = Polynomial{Rational{Int}, :h}
        h = variable(R)

        A = sparse(
            [1, 5, 9, 1, 5, 10, 1, 7, 9, 2, 7, 10, 3, 5, 9, 3, 6, 10, 3, 7, 11, 4, 8, 12], 
            [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8], 
            [h, h, h, R(1), R(1), h, R(1), h, R(1), R(1), R(1), R(1), h, R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1)], 
            12, 8
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        @test d == [R(1), R(1), R(1), R(1), R(1), R(1), R(1), h^2]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-no-preprocess-F2[h]" begin
        F2 = GaloisField(2)
        R = Polynomial{F2, :h}
        h = variable(R)

        A = sparse(
            [1, 2, 3, 5, 7, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8], 
            [1, 2, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 12, 12], 
            [R(1), R(1), R(1), R(1), h, R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h], 
            8, 12
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        @test d == [R(1), R(1), R(1), R(1), R(1), R(1), h, h]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-preprocess-simple" begin 
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 3 1 0 1 0 2 0;
            0 0 1 1 0 0 0 5 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 2 0 1 1 0 1 1
        ])
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        @test d == [1, 1, 1, 1, 1, 1]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-preprocess-random" begin 
        density = 0.1
        (m, n) = (30, 40)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-preprocess-Z" begin 
        # d[3] of CKh(R-trefoil)
        A = sparse(
            [1, 2, 3, 5, 4, 6, 1, 2, 5, 3, 4, 7, 1, 2, 3, 5, 6, 7], 
            [1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12], 
            [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1], 
            8, 12
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        @test d == [1, 1, 1, 1, 1, 1, 2]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-preprocess-Q[h]" begin
        R = Polynomial{Rational{Int}, :h}
        h = variable(R)

        A = sparse(
            [1, 5, 9, 1, 5, 10, 1, 7, 9, 2, 7, 10, 3, 5, 9, 3, 6, 10, 3, 7, 11, 4, 8, 12], 
            [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8], 
            [h, h, h, R(1), R(1), h, R(1), h, R(1), R(1), R(1), R(1), h, R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1), R(1)], 
            12, 8
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        @test d == [R(1), R(1), R(1), R(1), R(1), R(1), R(1), h^2]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-preprocess-F2[h]" begin
        F2 = GaloisField(2)
        R = Polynomial{F2, :h}
        h = variable(R)

        A = sparse(
            [1, 2, 3, 5, 7, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8], 
            [1, 2, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 12, 12], 
            [R(1), R(1), R(1), R(1), h, R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h, R(1), R(1), R(1), h], 
            8, 12
        )
        (m, n) = size(A)
        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        @test d == [R(1), R(1), R(1), R(1), R(1), R(1), h, h]
        @test P * A * Q == spdiagm(m, n, d)
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end
end

nothing