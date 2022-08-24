using Test

@testset "SNF" begin
    using SparseArrays: sparse
    using Knots.Matrix: SNF, snf, print_matrix
    
    # HNF

    @testset "snf-no-preprocess" begin 
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 3 1 0 1 0 2 0;
            0 0 1 1 0 0 0 5 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 2 0 1 1 0 1 1
        ])

        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=false, flags=(true, true, true, true))

        B = P * A * Q

        @test isone(P * Pinv)
        @test isone(Q * Qinv)
        @test all( map(i -> B[i, i] == d[i], 1:length(d) ))
    end

    @testset "snf-preprocess" begin 
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 3 1 0 1 0 2 0;
            0 0 1 1 0 0 0 5 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 2 0 1 1 0 1 1
        ])

        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        B = P * A * Q

        @test isone(P * Pinv)
        @test isone(Q * Qinv)
        @test all( map(i -> B[i, i] == d[i], 1:length(d) ))
    end

    @testset "snf-preprocess-random" begin 
        density = 0.1
        (m, n) = (60, 50)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        (d, P, Pinv, Q, Qinv) = snf(A; preprocess=true, flags=(true, true, true, true))

        B = P * A * Q

        @test isone(P * Pinv)
        @test isone(Q * Qinv)
        @test all( map(i -> B[i, i] == d[i], 1:length(d) ))
    end

end

nothing