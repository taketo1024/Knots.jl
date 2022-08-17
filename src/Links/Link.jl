#typealiases

const Edge = Int
const Resolution = Int
const State = Vector{Resolution}

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
    elseif x.type == V
        3 - j
    else
        (5 - j) % 4
    end
    p(i - 1) + 1
end

function strands(x::Crossing) :: Vector{Vector{Edge}}
    if x.type ∈ [X⁺, X⁻]
        [[x.edges[1], x.edges[3]], [x.edges[2], x.edges[4]]]
    elseif x.type == V
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

function isConnectable(c1::Component, c2::Component) :: Bool
    if c1.isClosed || c2.isClosed
        return false
    end

    (c1_head, c1_tail) = (first(c1.edges), last(c1.edges))
    (c2_head, c2_tail) = (first(c2.edges), last(c2.edges))
    ends = Set([c1_head, c1_tail, c2_head, c2_tail])

    length(ends) < 4
end

function connect(c1::Component, c2::Component) :: Component
    if c1.isClosed || c2.isClosed
        throw(Exception)
    end

    (c1_head, c1_tail) = (first(c1.edges), last(c1.edges))
    (c2_head, c2_tail) = (first(c2.edges), last(c2.edges))

    edges = if c1_tail == c2_head
        vcat(c1.edges,  c2.edges[2:end])
    elseif c1_tail == c2_tail
        vcat(c1.edges, reverse(c2.edges[1:end-1]))
    elseif c1_head == c2_tail
        vcat(c2.edges[1:end-1], c1.edges)
    elseif c1_head == c2_head
        vcat(reverse(c2.edges[2:end]), c1.edges)
    else 
        throw(Exception)
    end

    if first(edges) != last(edges)
        Component(edges, false)
    else
        Component(edges[1:end-1], true)
    end
end

Base.show(io::IO, c::Component) = begin
    str = join(c.edges, "-")
    print(io, c.isClosed ? "[$str]" : "($str)")
end 

Base.:(==)(c1::Component, c2::Component) = begin 
    (c1.edges, c1.isClosed) == (c2.edges, c2.isClosed)
end


# Link

struct Link
    data::Vector{Crossing}
end

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
    sum(+, _crossingSigns(l), init=0)
end

function signedCrossingNums(l::Link) :: Tuple{Int, Int}
    e = _crossingSigns(l)
    (count(x -> x > 0, e), count(x -> x < 0, e))
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
    result = Dict(map(crossings(l)) do x 
        (x, 0)
    end)

    function _traverse!(queue::Vector{Crossing}, sindex::Int)
        done = Set{Crossing}()
        while !isempty(queue)
            start = popfirst!(queue)
            x = start
            i = sindex

            while true
                if i ∈ [2, 4]
                    result[x] = ((x.type == X⁺) == (i == 2)) ? 1 : -1
                end

                next = _next(l, x, pass(x, i))
                if next == (start, sindex) || next ≡ nothing
                    break
                else
                    (x, i) = next
                    if x ∈ done 
                        break
                    end
                end
            end
        end
    end

    queue = l.data[:]
    _traverse!(queue, 1)

    if 0 ∈ values(result)
        remain = filter(result) do (_, val)
            val == 0
        end
        queue = collect(keys(remain))
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
        for str in strands(x)
            c = (str[1] == str[2]) ? 
                Component([str[1]], true) : 
                Component(str, false)

            i = findfirst(comps) do c2
                isConnectable(c, c2) 
            end

            if i ≠ nothing 
                comps[i] = connect(comps[i], c)

                j = findfirst(1 : length(comps)) do j
                    i ≠ j && isConnectable(comps[i], comps[j])
                end
                if j ≠ nothing 
                    comps[i] = connect(comps[i], comps[j])
                    deleteat!(comps, j)
                end
            else 
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

function resolve(l::Link, s::State) :: Link
    @assert length(s) <= crossingNum(l)
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

# constants
const emptyLink = Link([])
const unknot = resolve(Link([1, 2, 2, 1]), [0])
const trefoil = Link([1,4,2,5],[3,6,4,1],[5,2,6,3])
const figure8 = Link([4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8])
const hopfLink = Link([4,1,3,2],[2,3,1,4])