using ..Links: State, Component, crossings, edges as l_edges, resolve, canonical_state

@enum SeifertCircleColor SeifertCircleColor_a SeifertCircleColor_b
const SeifertCircleColors = (SeifertCircleColor_a, SeifertCircleColor_b)

function canonical_cycle(A::KhAlgStructure{R}, l::Link; positive=true) :: KhChain{R} where {R}
    @assert iszero(A.t)

    s = canonical_state(l)
    circles = components(resolve(l, s))
    colors = _color_circles(l, circles)

    _canonical_cycle(A, s, colors; positive=positive)
end

function canonical_cycles(A::KhAlgStructure{R}, l::Link) :: Tuple{KhChain{R}, KhChain{R}} where {R}
    (canonical_cycle(A, l), canonical_cycle(A, l; positive=false))
end

function _color_circles(l::Link, circles::Vector{Component}) :: Vector{SeifertCircleColor}
    isempty(circles) && return SeifertCircleColor[]

    (a, b) = SeifertCircleColors
    e0 = minimum(l_edges(l))
    i0 = findfirst( c -> e0 in c.edges, circles )
    c0 = circles[i0]
    
    colors = Dict(c0 => a)
    queue = [c0]

    remain = copy(circles)
    deleteat!(remain, i0)

    while !isempty(queue)
        c1 = popfirst!(queue)

        # crossings that touch `circle`.
        xs = filter(crossings(l)) do x 
            any(map(x.edges) do e
                e in c1.edges
            end)
        end

        # find circles that touch `circle`
        for x in xs 
            e_i = findfirst( e -> !(e in c1.edges), x.edges )
            e = x.edges[e_i]

            c2_i = findfirst( c -> e in c.edges, remain )
            isnothing(c2_i) && continue

            c2 = remain[c2_i]
            deleteat!(remain, c2_i)

            colors[c2] = (colors[c1] == a) ? b : a
            push!(queue, c2)
        end
    end

    @assert isempty(remain)

    map(c -> colors[c], circles)
end

function _canonical_cycle(A::KhAlgStructure{R}, s::State, colors::Vector{SeifertCircleColor}; positive=true) :: KhChain{R} where {R}
    a = SeifertCircleColor_a
    h = A.h
    one = Base.one(R)

    # collect factors
    factors = map(colors) do c
        (c == a) == positive ? 
        [(X, one)] : 
        [(X, one), (I, -h)]
    end

    # construct tensor product
    labels = [KhAlgGenerator[]]
    coeffs = [one]

    for factor in factors
        l = length(labels)
        t = length(factor) # 1 or 2
        if t == 2 
            append!(labels, deepcopy(labels))
            append!(coeffs, coeffs)
        end

        for i in 1 : t
            (x, r) = factor[i]
            for j in 1 : l
                k = (i - 1) * l + j # ∈ 1 : 2l
                push!(labels[k], x)
                coeffs[k] *= r
            end
        end
    end
    
    KhChain(Dict(
        KhChainGenerator(s, label) => coeff 
        for (label, coeff) in zip(labels, coeffs)
    ))
end