module LinkTests
    using Test
    using Khovanov: Link, emptyLink, unknot, isEmpty, crossingNum, signedCrossingNums, writhe, components

    function testAll()
        @testset "Links.jl" begin
            testEmptyLink()
            testUnknot()
            testOrientation()
        end
    end

    function testEmptyLink()
        l = emptyLink
        @test isEmpty(l)
        @test crossingNum(l) == 0
        @test writhe(l) == 0
        @test length(components(l)) == 0
    end

    function testUnknot()
        l = unknot
        println(l)
        @test !isEmpty(l)
        @test crossingNum(l) == 0
        @test writhe(l) == 0
        @test length(components(l)) == 1
    end

    function testOrientation() 
        @test signedCrossingNums(Link([0, 0, 1, 1])) == (1, 0)
        @test signedCrossingNums(Link([0, 1, 1, 0])) == (0, 1)
    end
end