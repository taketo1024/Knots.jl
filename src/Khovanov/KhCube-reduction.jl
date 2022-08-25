using SparseArrays: findnz
using ..Matrix: snf_iterate_pivots
using ..Homology: generate_matrix

# Cube reduction 

function reduce!(cube::KhCube)
    n = dim(cube)
    for k in 0:n-1
        reduce!(cube, k)
    end
end

function reduce!(cube::KhCube{R}, k::Int) where {R}
    @debug "cube-reduce (k: $k)"

    Gₖ = chain_generators(cube, k)
    Gₖ₊₁ = chain_generators(cube, k + 1)

    A = generate_matrix(Gₖ, Gₖ₊₁, R) do x 
        differentiate(cube, x)
    end

    (F, S, p, q) = snf_iterate_pivots(A; flags=(false, false, false, false))
    r = length(F.factors)
    r == 0 && return 

    targets = map( i -> (p(i), q(i)), 1:r)

    for (i, j) in targets
        x = Gₖ[j]
        y = Gₖ₊₁[i]
        remove_generator!(cube, x)
        remove_generator!(cube, y)
    end

    @debug "cancelled $(length(targets)) pairs."

    n = length(Gₖ)
    for j in 1 : n - r
        x = Gₖ[q(j + r)]
        ys = Pair{KhChainGenerator, R}[]

        Sⱼ = S[:, j]
        for (i, a) in zip(findnz(Sⱼ)...)
            iszero(a) && continue
            y = Gₖ₊₁[p(i + r)]

            push!(ys, (y => a))
        end

        cube.targets[x] = ys
    end
end

function remove_generator!(cube::KhCube, x::KhChainGenerator)
    u = x.state
    V = vertex(cube, u)
    i = findfirst(V.generators, x)
    deleteat!(V.generators, i)
end