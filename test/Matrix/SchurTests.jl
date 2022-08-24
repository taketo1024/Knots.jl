using Test

@testset "Schur" begin
    using SparseArrays
    using Knots.Matrix: pivot, schur_complement, print_matrix
    
    @testset "test" begin
        density = 0.1
        (m, n) = (50, 50)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        piv = pivot(A)
        (r, S, T) = schur_complement(A, piv; flags=(true,true,true,true))
        (P, Pinv, Q, Qinv) = T

        B = dropzeros(P * A * Q)

        @test isone(B[1:r, 1:r])
        @test iszero(B[1:r, r+1:n])
        @test iszero(B[r+1:m, 1:r])
        @test B[r+1:m, r+1:n] == S
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end
end