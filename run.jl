using Knots
using Knots.Links
using Knots.Khovanov
using Knots.Homology: compute_single, compute_incremental, compute_reverse_incremental

function run(name::String, str::String="Kh"; mode::Int=1, only_zero::Bool=false)
    println("run: $name, mode: $mode")
    l = Link(name)
    A = KhAlgStructure(str)
    H = KhHomology(A, l)
    prev = nothing
    range = only_zero ? 
        (0:0) : 
        (mode in [1, 2]) ? hDegRange(H) : reverse(hDegRange(H))

    for i in range
        if mode == 1
            Hi = compute_single(H, i)
        elseif mode == 2
            Hi, prev = compute_incremental(H, i; previous=prev)
        elseif mode == 3
            Hi, prev = compute_reverse_incremental(H, i; previous=prev)
        end
        println("H[$i] = ", Hi)
    end
end

@show Threads.nthreads()
@time run("3_1")

nothing

