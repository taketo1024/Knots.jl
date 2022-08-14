using Test

@testset "KhChain.jl" begin
    using Khovanov: KhChainGenerator, KhChain, Kh, X, I, simplify!

    R = Int
    C = KhChain{R}

    @testset "generator-equality" begin
        x = KhChainGenerator([0], [X])
        y = KhChainGenerator([0], [X])
        z = KhChainGenerator([0], [I])
        @test x == y
        @test hash(x) == hash(y)
    end

    @testset "generator-inequality" begin
        x = KhChainGenerator([0], [X])
        y = KhChainGenerator([0], [I])
        z = KhChainGenerator([1], [X])
        @test x ≠ y
        @test x ≠ z
    end

    @testset "zero" begin
        c = zero(C)
        @test length(c) == 0
    end

    @testset "singleton" begin
        x = KhChainGenerator([0], [X])
        c = C(Dict(x => 1))
        @test length(c) == 1
    end

    @testset "index-access" begin
        x = KhChainGenerator([0], [I])
        c = zero(C)
        @test c[x] == 0
        c[x] = 1
        @test c[x] == 1
        c[x] = 0
        @test c[x] == 0
    end

    @testset "add" begin
        x = KhChainGenerator([0], [I])
        y = KhChainGenerator([0], [X])
        c = C(Dict(x => 1))
        d = C(Dict(x => 3, y => 2))
        e = c + d
        @test length(e) == 2
        @test e[x] == 4
        @test e[y] == 2
    end

    @testset "sub" begin
        x = KhChainGenerator([0], [I])
        y = KhChainGenerator([0], [X])
        c = C(Dict(x => 1))
        d = C(Dict(x => 3, y => 3))
        e = c - d
        @test length(e) == 2
        @test e[x] == -2
        @test e[y] == -3
    end

    @testset "scalar-mul" begin
        x = KhChainGenerator([0], [I])
        y = KhChainGenerator([0], [X])
        c = C(Dict(x => 4, y => -3))
        d = 2c
        @test length(d) == 2
        @test d[x] == 8
        @test d[y] == -6
    end
    
    @testset "equal" begin
        x = KhChainGenerator([0], [I])
        y = KhChainGenerator([0], [X])
        c = C(Dict(x => 3, y => 2))
        d = C(Dict(x => 3, y => 2))
        e = C(Dict(x => 3, y => 1))
        @test c == d
        @test c ≠ e
    end
    
    @testset "simplify" begin
        x = KhChainGenerator([0], [I])
        c = C(Dict(x => 0))

        @test length(c) == 1
        @test c ≠ zero(C)
        
        simplify!(c)
        
        @test length(c) == 0
        @test c == zero(C)
    end
end