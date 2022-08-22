using Test

@testset "SNF" begin
    using AbstractAlgebra
    using SparseArrays: sparse
    using Knots.Matrix: SNF, snf, print_matrix
    using Knots.Matrix: kb_canonical_row_x!, kb_reduce_column_x!, hnf_x, snf_kb_x_step1!, snf_kb_x_step2!, snf_x
    
    Z = zz
    Q = qq

    # HNF

    @testset "kb_canonical_row_x!" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([
            2 -1 -2 -2 -3;
            0 0 (-1) 2 5;
            2 -2 -4 -3 -6;
            1 7 1 5 3;
            1 -12 -6 -10 -11
        ])
        B = deepcopy(A)
        P = identity_matrix(B, 5)
        Pinv = identity_matrix(B, 5)
        kb_canonical_row_x!(B, P, Pinv, 2, 3)

        @test P * A == B
        @test Pinv * B == A
        @test isone(P * Pinv)
    end

    @testset "kb_reduce_column_x!" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([
            0  (3) 1  5  0;
            (1) 2  2 -7  2;
            0   0  0  0 (2);
            0   0  0 (3) 1;
            1   2  3  4  5
        ])

        B = deepcopy(A)
        P = identity_matrix(B, 5)
        Pinv = identity_matrix(B, 5)
        pivot = [2, 1, 0, 4, 5]
        kb_reduce_column_x!(B, P, Pinv, pivot, 4)

        @test P * A == B
        @test Pinv * B == A
        @test isone(P * Pinv)
    end

    @testset "hnf_x" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([2 -1 -2 -2 -3;
            1 2 -1 1 -1;
            2 -2 -4 -3 -6;
            1 7 1 5 3;
            1 -12 -6 -10 -11])
        (H, P, Pinv) = hnf_x(A)

        @test is_hnf(H)
        @test P * A == H
        @test Pinv * H == A
        @test isone(P * Pinv)
    end

    @testset "hnf_x_ZZ" begin
        M = MatrixSpace(ZZ, 3, 3)
        A = M([4 6 2; 0 0 10; 0 5 3])
        H, P, Pinv = hnf_x(A)
     
        @test H == M([4 1 9; 0 5 3; 0 0 10])
        @test is_hnf(H)
        @test is_unit(det(P))
        @test is_unit(det(Pinv))
        @test P * A == H
        @test Pinv * H == A
        @test isone(P * Pinv)
     end

     @testset "hnf_x_Qx" begin
        R, x = PolynomialRing(QQ, "x")
        M = MatrixSpace(R, 4, 3)
        A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))
        H, U, Uinv = hnf_x(A)
        
        @test is_hnf(H)
        @test is_unit(det(U))
        @test is_unit(det(Uinv))
        @test U * A == H
        @test Uinv * H == A
        @test isone(U * Uinv)
     end

     @testset "hnf_F7[x]" begin
        # Fake up finite field of char 7, degree 2
        R, x = PolynomialRing(GF(7), "x")
        F = ResidueField(R, x^2 + 6x + 3)
        a = F(x)
     
        S, y = PolynomialRing(F, "y")
        N = MatrixSpace(S, 3, 4)
        A = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))
     
        H, U, Uinv = hnf_x(A)

        @test is_hnf(H)
        @test is_unit(det(U))
        @test is_unit(det(Uinv))
        @test U * A == H
        @test Uinv * H == A
        @test isone(U * Uinv)
    end

    # SNF

    @testset "snf_kb_x_step1!" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([4, 6, -18, -15, -46, -1, 0, 6, 4, 13, -13, -12, 36, 30, 97, -7, -6, 18, 15, 49, -6, -6, 18, 15, 48])

        B = deepcopy(A)
        P = identity_matrix(B, 5)
        Pinv = identity_matrix(B, 5)
        Q = identity_matrix(B, 5)
        Qinv = identity_matrix(B, 5)
        snf_kb_x_step1!(B, P, Pinv, Q, Qinv)

        @test P * A * Q == B
        @test Pinv * B * Qinv == A
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf_kb_x_step2!" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([10, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 42, 0, 0, 0, 0, 0, 80, 0, 0, 0, 0, 0, 30])

        B = deepcopy(A)
        P = identity_matrix(B, 5)
        Pinv = identity_matrix(B, 5)
        Q = identity_matrix(B, 5)
        Qinv = identity_matrix(B, 5)
        snf_kb_x_step2!(B, P, Pinv, Q, Qinv)

        @test(is_snf(B))
        @test P * A * Q == B
        @test Pinv * B * Qinv == A
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf_x" begin
        M = MatrixSpace(Z, 5, 5)
        A = M([2 -1 -2 -2 -3;
            1 2 -1 1 -1;
            2 -2 -4 -3 -6;
            1 7 1 5 3;
            1 -12 -6 -10 -11])
        (B, P, Pinv, Q, Qinv) = snf_x(A)

        @test is_snf(B)
        @test P * A * Q == B
        @test Pinv * B * Qinv == A
        @test isone(P * Pinv)
        @test isone(Q * Qinv)
    end

    @testset "snf-Q[x]" begin
        R, x = PolynomialRing(QQ, "x")
        M = MatrixSpace(R, 4, 3)
        A = M(map(R, Any[0 0 0; x^3+1 x^2 0; 0 x^2 x^5; x^4+1 x^2 x^5+x^3]))
     
        S, P, Pinv, Q, Qinv = snf_x(A)
        @test is_snf(S)
        @test is_unit(det(P))
        @test is_unit(det(Pinv))
        @test is_unit(det(Q))
        @test is_unit(det(Qinv))
        @test P * A * Q == S
        @test Pinv * S * Qinv == A
     end

     @testset "snf-F7[x]" begin
        # Fake up finite field of char 7, degree 2
        R, x = PolynomialRing(GF(7), "x")
        F = ResidueField(R, x^2 + 6x + 3)
        a = F(x)
     
        S, y = PolynomialRing(F, "y")
        N = MatrixSpace(S, 3, 4)
        A = N(map(S, Any[1 0 a 0; a*y^3 0 3*a^2 0; y^4+a 0 y^2+y 5]))
     
        S, P, Pinv, Q, Qinv = snf_x(A)
        @test is_snf(S)
        @test is_unit(det(P))
        @test is_unit(det(Pinv))
        @test is_unit(det(Q))
        @test is_unit(det(Qinv))
        @test P * A * Q == S
        @test Pinv * S * Qinv == A
     end

     @testset "hnf_x-flags" begin
        M = MatrixSpace(Z, 2, 2)
        A = M([1, 0, 0, 1])

        (H, P, Pinv) = hnf_x(A; flags=(true, true))
        @test !isnothing(P)
        @test !isnothing(Pinv)

        (H, P, Pinv) = hnf_x(A; flags=(true, false))
        @test !isnothing(P)
        @test isnothing(Pinv)

        (H, P, Pinv) = hnf_x(A; flags=(false, true))
        @test isnothing(P)
        @test !isnothing(Pinv)

        (H, P, Pinv) = hnf_x(A; flags=(false, false))
        @test isnothing(P)
        @test isnothing(Pinv)
    end

    @testset "snf_x-flags" begin
        M = MatrixSpace(Z, 2, 2)
        A = M([1, 0, 0, 1])

        (S, P, Pinv, Q, Qinv) = snf_x(A; flags=(true, true, true, true))
        @test !isnothing(P)
        @test !isnothing(Pinv)
        @test !isnothing(Q)
        @test !isnothing(Qinv)

        (S, P, Pinv, Q, Qinv) = snf_x(A; flags=(true, true, false, false))
        @test !isnothing(P)
        @test !isnothing(Pinv)
        @test isnothing(Q)
        @test isnothing(Qinv)

        (S, P, Pinv, Q, Qinv) = snf_x(A; flags=(false, false, true, true))
        @test isnothing(P)
        @test isnothing(Pinv)
        @test !isnothing(Q)
        @test !isnothing(Qinv)

        (S, P, Pinv, Q, Qinv) = snf_x(A; flags=(false, false, true, false))
        @test isnothing(P)
        @test isnothing(Pinv)
        @test !isnothing(Q)
        @test isnothing(Qinv)
    end

    @testset "snf-no-preprocess" begin 
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 3 1 0 1 0 2 0;
            0 0 1 1 0 0 0 5 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 2 0 1 1 0 1 1
        ])

        F = snf(A; preprocess=false, flags=(true, true, true, true))
        (d, P, Pinv, Q, Qinv) = (F.S, F.P, F.P⁻¹, F.Q, F.Q⁻¹)

        B = P * A * Q

        @test is_one(P * Pinv)
        @test is_one(Q * Qinv)
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

        F = snf(A; preprocess=true, flags=(true, true, true, true))
        (d, P, Pinv, Q, Qinv) = (F.S, F.P, F.P⁻¹, F.Q, F.Q⁻¹)

        B = P * A * Q

        @test is_one(P * Pinv)
        @test is_one(Q * Qinv)
        @test all( map(i -> B[i, i] == d[i], 1:length(d) ))
    end

    @testset "snf-preprocess-random" begin 
        density = 0.1
        (m, n) = (60, 50)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        F = snf(A; preprocess=true, flags=(true, true, true, true))

        (d, P, Pinv, Q, Qinv) = (F.S, F.P, F.P⁻¹, F.Q, F.Q⁻¹)

        B = P * A * Q

        @test is_one(P * Pinv)
        @test is_one(Q * Qinv)
        @test all( map(i -> B[i, i] == d[i], 1:length(d) ))
    end

end

nothing