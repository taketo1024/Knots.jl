using Test

@testset "Pivot" begin
    using SparseArrays
    using OrderedCollections
    using Permutations: Permutation
    using Knots.Matrix: Pivot, pivot, npivots, permutations, findFLPivots!, findFLColumnPivots!, findCycleFreePivots!, sortPivots!, occupiedCols, permutation

    function ok(B, r) 
        for i in 1 : r
            for j in 1 : i - 1
                !iszero(B[i, j]) && return false
            end
        end
        true
    end

    @testset "initialize" begin
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 1 1 0 1 0 2 0;
            0 0 1 1 0 0 0 1 1;
            0 1 1 0 3 0 0 0 0;
            0 1 0 1 0 0 1 0 1;
            1 0 1 0 1 1 0 1 1
        ])

        piv = Pivot(A)

        @test piv.entries == [[1,3,6,7,9], [2,3,4,6,8], [3,4,8,9], [2,3,5], [2,4,7,9],[1,3,5,6,8,9]]
        @test piv.rowHead == [1,2,3,2,2,1]
        @test piv.rowWeight == [5,5,4,3,4,6]
        @test piv.colWeight == [2,3,5,3,2,3,2,3,4]
        @test piv.candidates == [Set([1,3,6,7,9]), Set([2,3,4,6]), Set([3,4,8,9]), Set([2,3]), Set([2,4,7,9]), Set([1,3,5,6,8,9])]
        @test isempty(piv.pivots)
    end

    @testset "findFLPivots" begin
        A = sparse([
            1 0 1 0 0 1 1 0 1; # 
            0 1 1 1 0 1 0 1 0; 
            0 0 1 1 0 0 0 1 1;
            0 1 1 0 1 0 0 0 0; #
            0 0 1 1 0 0 0 0 0; #
            0 0 0 0 0 1 0 1 1  #
        ])
        piv = Pivot(A)

        findFLPivots!(piv)

        @test piv.pivots == OrderedDict(1 => 1, 2 => 4, 3 => 5, 6 => 6)
        @test occupiedCols(piv) == Set([1, 2, 3, 4, 5, 6, 7, 8, 9])
    end

    @testset "findFLColumnPivots" begin
        A = sparse([
            1 0 0 0 0 1 0 0 1; #
            0 1 1 1 0 1 0 1 0; 
            0 0 1 1 0 0 0 1 1; #
            0 1 0 0 1 0 0 0 0; #
            0 0 1 0 0 0 0 0 0; #
            0 1 0 0 0 1 0 1 0  #
        ])
        piv = Pivot(A)

        findFLPivots!(piv)

        @test piv.pivots == OrderedDict(1 => 1, 2 => 4, 3 => 5)
        @test occupiedCols(piv) == Set([1, 2, 3, 5, 6, 9])

        findFLColumnPivots!(piv)

        @test piv.pivots == OrderedDict(1 => 1, 2 => 4, 3 => 5, 8 => 6, 4 => 3)
        @test occupiedCols(piv) == Set([1, 2, 3, 4, 5, 6, 8, 9])
    end


    function pa(A)
        Base.print_matrix(stdout, Array(A), "[", " ", "]")
        println()
    end

    @testset "findCycleFreePivots" begin
        A = sparse([
            1 0 0 0 0 1 0 0 1; #
            0 1 1 1 0 1 0 1 0; 
            0 0 1 1 0 0 0 1 1; #
            0 1 0 0 1 0 0 0 0; #
            0 0 1 0 0 0 0 0 0; #
            0 1 0 0 0 1 0 1 0  #
        ])
        piv = Pivot(A)

        findCycleFreePivots!(piv)
        sortPivots!(piv)

        r = npivots(piv)
        (p, q) = permutations(piv)

        B = permute(A, p.data, q.data)

        @test r >= 5
        @test ok(B, r)
    end

    @testset "test-random" begin
        density = 0.1
        (m, n) = (100, 100)
        A = sparse([ rand() < density ? 1 : 0 for i in 1:m, j in 1:n])

        piv = pivot(A)
        r = npivots(piv)
        (p, q) = permutations(piv)
        
        B = permute(A, p.data, q.data)

        @test r > 30
        @test ok(B, r)
    end

    @testset "permutation" begin
        p = permutation([7, 3, 4, 1, 5], 10)
        @test p == Permutation([7, 3, 4, 1, 5, 2, 6, 8, 9, 10])
    end
end

nothing