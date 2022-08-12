module Links
    export Link

    Edge = Int
    @enum Resolution R0 R1

    # Crossing

    @enum CrossingType X⁺ X⁻ V H 
    
    struct Crossing 
        type :: CrossingType
        edges :: Vector{Edge}
    end

    function isResolved(x::Crossing)
        x.type ∈ [V, H]
    end

    function resolve(x::Crossing, r::Resolution)
        if (x.type, r) ∈ [(X⁺, R0), (X⁻, R1)]
            Crossing(V, x.edges)
        elseif  (x.type, r) ∈ [(X⁺, R1), (X⁻, R0)]
            Crossing(H, x.edges)
        else
            throw(Exception)
        end
    end

    function mirror(x::Crossing)
        t = if x.type == X⁺
            X⁻
        elseif x.type == X⁻
            X⁺
        else
            x.type
        end
        Crossing(t, x.edges)
    end

    function pass(x::Crossing, i::Int) 
        if x.type ∈ [X⁺, X⁻]
            (i + 2) % 4
        elseif x.type == H
            3 - i
        else
            (5 - i) % 4
        end
    end

    Base.show(io::IO, x::Crossing) = print(io, "$(x.type)$(x.edges)")

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

    Base.show(io::IO, l::Link) = print(io, "Link($(join(l.data, ", ")))")
end