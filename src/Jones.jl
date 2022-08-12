module Jones
    using Polynomials
    using ..Links: Link, unknot
    
    function kauffmanBracket(l::Link) :: LaurentPolynomial
        A = LaurentPolynomial([0, 1], 0, "A")
    end

end # module Jones