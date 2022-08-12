module Links
    export Link

    Edge = Int
    Resolution = Int

    # Crossing

    @enum CrossingType X⁺ X⁻ V H 
    
    struct Crossing 
        type :: CrossingType
        edges :: Vector{Edge}
    end

    function isResolved(x::Crossing) :: Bool
        x.type ∈ [V, H]
    end

    function resolve(x::Crossing, r::Resolution) :: Crossing
        @assert r ∈ [0, 1]
        
        if (x.type, r) ∈ [(X⁺, 0), (X⁻, 1)]
            Crossing(V, x.edges)
        elseif  (x.type, r) ∈ [(X⁺, 1), (X⁻, 0)]
            Crossing(H, x.edges)
        else
            throw(Exception)
        end
    end

    function mirror(x::Crossing) :: Crossing
        t = if x.type == X⁺
            X⁻
        elseif x.type == X⁻
            X⁺
        else
            x.type
        end
        Crossing(t, x.edges)
    end

    function pass(x::Crossing, i::Int) :: Int
        p(j) = if x.type ∈ [X⁺, X⁻]
            (j + 2) % 4
        elseif x.type == H
            3 - j
        else
            (5 - j) % 4
        end
        p(i - 1) + 1
    end

    function strands(x::Crossing) :: Vector{Vector{Edge}}
        if x.type ∈ [X⁺, X⁻]
            [[x.edges[1], x.edges[3]], [x.edges[2], x.edges[4]]]
        elseif x.type == H
            [[x.edges[1], x.edges[4]], [x.edges[2], x.edges[3]]]
        else
            [[x.edges[1], x.edges[2]], [x.edges[3], x.edges[4]]]
        end
    end

    Base.show(io::IO, x::Crossing) = print(io, "$(x.type)$(x.edges)")

    # Component

    struct Component
        edges::Vector{Edge}
        isClosed::Bool
    end

    function isConnectable(c::Component, arc::Vector{Edge}) :: Bool
        if c.isClosed
            return false
        end
        ends = [first(c.edges), last(c.edges)]
        first(arc) ∈ ends || last(arc) ∈ ends
    end

    function isConnectable(c1::Component, c2::Component) :: Bool
        isConnectable(c1, c2.edges)
    end

    function connect(c::Component, arc::Vector{Edge}) :: Component
        if c.isClosed
            throw(Exception)
        end

        edges = if last(c.edges) == first(arc) 
            vcat(c.edges,  arc[2:end])
        elseif last(c.edges) == last(arc)
            vcat(c.edges, reverse(arc[1:end-1]))
        elseif first(c.edges) == last(arc)
            vcat(arc[1:end-1], c.edges)
        elseif first(c.edges) == first(c.edges)
            vcat(reverse(arc[2:end]), c.edges)
        else 
            throw(Exception)
        end

        if first(edges) != last(edges)
            Component(edges, false)
        else
            Component(edges[1:end-1], true)
        end
    end

    function connect(c1::Component, c2::Component) :: Component
        connect(c1, c2.edges)
    end

    Base.show(io::IO, c::Component) = begin
        str = join(c.edges, "-")
        print(io, c.isClosed ? "[$str]" : "($str)")
    end 

    # Link

    struct Link
        data::Vector{Crossing}
    end

    emptyLink = Link([])

    # Planer Diagram code, represented by crossings:
    # 
    #     3   2
    #      \ /
    #       \      = (0, 1, 2, 3)
    #      / \
    #     0   1
    # 
    # The lower edge has direction 0 -> 2.
    # The crossing is +1 if the upper goes 3 -> 1.
    # see: http://katlas.math.toronto.edu/wiki/Planar_Diagrams

    function Link(pdCode::Vector{Vector{Int}})
        crossings = map(pdCode) do edges
            @assert length(edges) == 4
            Crossing(X⁻, edges)
        end
        Link(crossings)
    end

    function Link(pdCode::Vector{Int}...)
        Link(collect(pdCode))
    end

    function isEmpty(l::Link) :: Bool
        isempty(l.data)
    end

    function crossingNum(l::Link) :: Int 
        count(l.data) do x
            !isResolved(x)
        end
    end

    function crossings(l::Link) :: Vector{Crossing}
        filter(l.data) do x 
            !isResolved(x)
        end
    end

    function writhe(l::Link) :: Int 
        sum(+, _crossingSigns(l))
    end

    function signedCrossingNums(l::Link) :: Tuple{Int, Int}
        e = _crossingSigns(l)
        (count(x -> x > 0, e), -count(x -> x < 0, e))
    end

    # search for the next crossing
    Optional{T} = Union{T, Nothing}
    function _next(l::Link, x::Crossing, i::Int) :: Optional{Tuple{Crossing, Int}} 
        e = x.edges[i]
        for y in l.data 
            for (j, f) in enumerate(y.edges)
                if e == f && (x ≠ y || (x == y && i ≠ j))
                    return (y, j)
                end
            end
        end
        nothing
    end

    function _crossingSigns(l::Link) :: Vector{Int}
        # traverse edges and determine orientation
        result = Dict{Crossing, Int}()

        function _traverse!(queue::Vector{Crossing}, sindex::Int)
            while !isempty(queue)
                start = popfirst!(queue)
                x = start
                i = sindex

                while true
                    if i ∈ [2, 4]
                        result[x] = ((x.type == X⁺) == (i == 2)) ? 1 : -1
                        deleteat!(queue, findall(queue) do y y == x end)
                    end
                    next = _next(l, x, pass(x, i))
                    if next == (start, sindex) || next ≡ nothing
                        break
                    else
                        (x, i) = next
                    end
                end
            end
        end

        queue = l.data[:]
        _traverse!(queue, 1)

        if 0 ∈ values(result)
            queue = findall(keys(result)) do x
                result[x] == 0
            end
            _traverse!(queue, 2)
        end

        @assert 0 ∉ values(result)

        map(crossings(l)) do x
            result[x]
        end
    end

    function components(l::Link) :: Vector{Component}
        comps = Vector{Component}()
        for x in l.data
            for s in strands(x)
                i = findfirst(comps) do c 
                    isConnectable(c, s) 
                end
                if i ≠ nothing 
                    comps[i] = connect(comps[i], s)

                    j = findfirst(1 : length(comps)) do j
                        i ≠ j && isConnectable(comps[i], comps[j])
                    end
                    if j ≠ nothing 
                        comps[i] = connect(comps[i], comps[j])
                        deleteat!(comps, j)
                    end
                else 
                    c = Component(s, false)
                    push!(comps, c)
                end
            end
        end
        comps
    end

    function mirror(l::Link) :: Link
        Link(map(l.data) do x
            mirror(x)
        end)
    end

    function _actualIndex(l::Link, i::Int) :: Int
        c = 0
        for (j, x) in enumerate(l.data)
            if !isResolved(x)
                c += 1
                if c >= i 
                    return j 
                end
            end
        end
        throw(Exception)
    end

    #     \ /           \ /
    #      /     <-->    \
    #     / \           / \
   
   function crossingChange(l::Link, i::Int) :: Link
        idx = _actualIndex(l, i)
        data = map(enumerate(l.data)) do (j, x)
            j == idx ? mirror(x) : x
        end
        Link(data)
    end

    #     \ /     0     \ /
    #      /     --->   | |
    #     / \           / \
    #
    #     \ /     1     \_/
    #      /     --->
    #     / \           /‾\

    function resolve(l::Link, i::Int, r::Resolution) :: Link
        idx = _actualIndex(l, i)
        data = map(enumerate(l.data)) do (j, x)
            j == idx ? resolve(x, r) : x
        end
        Link(data)
     end

    function resolve(l::Link, s::Vector{Resolution}) :: Link
        i = 0
        data = map(enumerate(l.data)) do (j, x) 
            if isResolved(x)
                x
            else
                i += 1
                resolve(x, s[i])
            end
        end
        Link(data)
    end

    Base.show(io::IO, l::Link) = print(io, "Link($(join(l.data, ", ")))")
end