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

    # function pass(x::Crossing, i::Int) :: Int
    #     if x.type ∈ [X⁺, X⁻]
    #         (i + 2) % 4
    #     elseif x.type == H
    #         3 - i
    #     else
    #         (5 - i) % 4
    #     end
    # end

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

    function crossings(l::Link) :: Vector{Crossing}
        filter(l.data) do x 
            !isResolved(x)
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

using .Links
l = Link([1,4,2,5],[3,6,4,1],[5,2,6,3])
l2 = Links.resolve(l, [1,1,1])
println(Links.components(l2))