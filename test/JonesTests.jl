module JonesTests
    using Test
    using Khovanov.Links: Link, unknot
    using Khovanov.Jones: kauffmanBracket

    function testAll()
        @testset "Jones.jl" begin
            test()
        end
    end

    function test()
        u = unknot
        b = kauffmanBracket(u)
        println(b)
    end

end # module JonesTests