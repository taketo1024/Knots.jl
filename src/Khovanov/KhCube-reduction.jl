using SparseArrays: findnz
using ..Matrix: pivot, npivots, coordinates, permutations, schur_complement, print_matrix
using ..Homology: generate_matrix

# Cube reduction 

function reduce!(cube::KhCube)
    n = dim(cube)
    for k in 0:n-1
        reduce!(cube, k)
    end
end

function reduce!(cube::KhCube{R}, k::Int) where {R}
    @debug "cube-reduce, k = $k"

    Gₖ = chain_generators(cube, k)
    Gₖ₊₁ = chain_generators(cube, k + 1)

    A = generate_matrix(Gₖ, Gₖ₊₁, R) do x 
        differentiate(cube, x)
    end

    piv = pivot(A)
    r = npivots(piv)

    r == 0 && return 

    (p, q) = permutations(piv)
    targets = coordinates(piv)

    for t in targets
        (i, j) = Tuple(t)
        x = Gₖ[j]
        y = Gₖ₊₁[i]
        remove_generator!(cube, x)
        remove_generator!(cube, y)
    end

    @debug "cancelled $(length(targets)) pairs."

    (_, S, _) = schur_complement(A, piv; flags=(false, false, false, false))

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

    reduce!(cube, k)
end

function remove_generator!(cube::KhCube, x::KhChainGenerator)
    u = x.state
    V = vertex(cube, u)
    i = findfirst(V.generators, x)
    deleteat!(V.generators, i)
end