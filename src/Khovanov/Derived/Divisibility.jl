export s_c

using Knots.Links: Link, writhe, n_seifert_circles
using Knots.Homology: vectorize, reduce!, compute_single

function divisibility(a::R, c::R) :: Int where {R}
    v = copy(a)
    k = 0
    while v % c == 0 
        v = div(v, c)
        k += 1
    end
    k
end

function s_c(l::Link, c::R; reduced=false) where {R}
    A = KhAlgStructure(c, zero(R))
    α = canonical_cycle(A, l)
    C = KhComplex(A, l; reduced=reduced, perform_reduction=false, with_transform=false)

    reduce!(C, -2)
    reduce!(C, -1; flags=(true, false, false, false)) # P
    reduce!(C,  0; flags=(false, false, false, true)) # Q⁻¹
    reduce!(C,  1)
    
    H = KhHomology(C)
    H0 = compute_single(H, 0; preprocess=false, with_transform=true)

    v = vectorize(H0, α)[1]
    k = divisibility(v, c)

    w = writhe(l)
    r = n_seifert_circles(l)

    2*k + w - r + 1
end