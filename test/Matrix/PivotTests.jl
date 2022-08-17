using Test

@testset "Pivot" begin
    using SparseArrays
    # using Knots.Matrix: Pivot, findPivots, findFLPivots

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

        @test piv.nonzero == [[1,3,6,7,9], [2,3,4,6,8], [3,4,8,9], [2,3,5], [2,4,7,9],[1,3,5,6,8,9]]
        @test piv.rowHead == [1,2,3,2,2,1]
        @test piv.rowWeight == [5,5,4,3,4,6]
        @test piv.colWeight == [2,3,5,3,2,3,2,3,4]
        @test piv.candidates == [Set([1,3,6,7,9]), Set([2,3,4,6]), Set([3,4,8,9]), Set([2,3]), Set([2,4,7,9]), Set([1,3,5,6,8,9])]
        @test isempty(piv.hasPivot)
        @test piv.pivots == fill(0, 9)
    end

    @testset "findFLPivots" begin
        A = sparse([
            1 0 1 0 0 1 1 0 1;
            0 1 1 1 0 1 0 1 0;
            0 0 1 1 0 0 0 1 1;
            0 1 1 0 1 0 0 0 0;
            0 0 1 1 0 0 0 0 0;
            0 0 0 0 0 1 0 1 1
        ])
        piv = Pivot(A)

        findFLPivots!(piv)

        @test piv.hasPivot == Set([1,4,5,6])
        @test piv.pivots == [1,4,5,0,0,6,0,0,0]
    end
end

nothing