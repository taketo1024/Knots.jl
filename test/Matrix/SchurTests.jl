using Test

@testset "Schur" begin
    using SparseArrays
    using Knots.Matrix: pivot, schur_complement, print_matrix, is_identity, is_zero
    
    @testset "test" begin
        density = 0.1
        (m, n) = (50, 50)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        piv = pivot(A)
        (r, S, T) = schur_complement(A, piv; flags=(true,true,true,true))
        (P, Pinv, Q, Qinv) = T

        B = dropzeros(P * A * Q)

        @test is_identity(B[1:r, 1:r])
        @test is_zero(B[1:r, r+1:n])
        @test is_zero(B[r+1:m, 1:r])
        @test B[r+1:m, r+1:n] == S
        @test is_identity(P * Pinv)
        @test is_identity(Q * Qinv)
    end
end