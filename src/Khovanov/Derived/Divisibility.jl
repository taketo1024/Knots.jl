export s_c

using ..Links: Link, components, writhe, n_seifert_circles
using ..Homology: vectorize, reduce!, compute_single, rank
using ..Extensions: symbol, isunit

function divisibility(a::R, c::R) :: Int where {R}
    @assert !iszero(c) && !isunit(c)
    iszero(a) && return typemax(Int)

    v = copy(a)
    k = 0
    while v % c == 0 
        v = div(v, c)
        k += 1
    end
    k
end

function s_c(l::Link, c::R; reduced=false) where {R}
    @debug "compute s_$c" l (R, c) reduced

    @assert !iszero(c)
    @assert !isunit(c)
    @assert length(components(l)) == 1 "only knots are supported."

    A = KhAlgStructure(c, zero(R))
    α = canonical_cycle(l, A)

    @debug "canonical cycle" α

    C = KhComplex(l, A; reduced=reduced, perform_reduction=false, with_transform=false)

    reduce!(C, -2)
    reduce!(C, -1; flags=(true, false, false, false)) # P
    reduce!(C,  0; flags=(false, false, false, true)) # Q⁻¹
    reduce!(C,  1)

    H = KhHomology(C)
    H0 = compute_single(H, 0; preprocess=false, with_transform=true)

    @debug "homology computed" H0
    @assert rank(H0) == (reduced ? 1 : 2) "invalid homology" 

    v = vectorize(H0, α)
    
    @debug "vectorized" v

    k = reduced ?
        divisibility(v[1], c) : 
        min(divisibility(v[1], c), divisibility(v[2], c))

    w = writhe(l)
    r = n_seifert_circles(l)
    s = 2k + w - r + 1

    @debug "result" k w r s

    s
end