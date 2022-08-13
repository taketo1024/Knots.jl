# KhCube

struct KhCubeVertex
    state::State
    circles::Vector{Component}
    generators::Vector{KhChainGenerator}
end

KhCubeVertex(l::Link, s::State) = begin
    circles = components(resolve(l, s))
    r = length(circles)
    generators = map(0 : 2^r - 1) do i in
        bits = digits(i, base=2, pad=r)
        label = map(b -> (b == 0) ? X : I, bits)
        KhChainGenerator(s, label)
    end
    KhCubeVertex(s, circles, generators)
end