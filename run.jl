using Knots
using Knots.Links
using Knots.Khovanov
using Knots.Homology: compute_single, compute_incremental, compute_reverse_incremental

function run1(H, msg="")
    isempty(msg) || println(msg)
    # for i in hDegRange(H)
    i = 0
        Hi = compute_single(H, i)
        println("H[$i] = ", Hi)
    # end
end

function run2(H, msg="")
    isempty(msg) || println(msg)
    prev = nothing
    # for i in hDegRange(H)
    i = 0
        Hi, prev = compute_incremental(H, i; previous=prev)
        println("H[$i] = ", Hi)
    # end
end

function run3(H, msg="")
    isempty(msg) || println(msg)
    H = KhHomology(A, l)
    prev = nothing
    # for i in reverse(hDegRange(H))
    i = 0
        Hi, prev = compute_reverse_incremental(H, i; previous=prev)
        println("H[$i] = ", Hi)
    # end
end

@show Threads.nthreads()

A = KhAlgStructure("Z-Kh")

@time run1(KhHomology(A, unknot), "test-run")

l = Link([4,2,5,1],[8,4,9,3],[5,12,6,13],[2,8,3,7],[9,19,10,18],[11,6,12,7],[13,22,14,1],[15,20,16,21],[17,11,18,10],[19,16,20,17],[21,14,22,15])
H = KhHomology(A, l)

@profview @time run1(H, "compute-single")
@profview @time run2(H, "compute-incremental")
@profview @time run3(H, "compute-reverse-incremental")

nothing

