module Utils
    function bitseq(length::Int, degree::Int) :: Vector{Int} 
        if length <= 0 || degree < 0 || length < degree 
            []
        elseif length > 1
            s₀ = map( b -> (b << 1) | 1, bitseq(length - 1, degree - 1) )
            s₁ = map( b -> b << 1, bitseq(length - 1, degree) )
            append!(s₀, s₁)
        else # 0 ≤ degree ≤ length == 1
            [degree] 
        end
    end
end