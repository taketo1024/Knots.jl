function subscript(i::Int) :: String
    c = i < 0 ? [Char(0x208B)] : []
    for d in reverse(digits(abs(i)))
        push!(c, Char(0x2080+d))
    end
    join(c)
end

function superscript(i::Int) :: String
    c = i < 0 ? c = [Char(0x207B)] : []
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    join(c)
end

using Polynomials

function symbol(R) :: String
    if R <: Integer
        "Z"
    elseif R <: Rational
        "Q"
    elseif R <: Polynomial
        S = typeof(coeffs(R(0))[1]) # isn't there a better way?
        x = variable(R)
        "$(symbol(S))[$x]"
    # elseif isa(R, AbstractAlgebra.PolyRing) || isa(R, AbstractAlgebra.MPolyRing)
    #     S = base_ring(R)
    #     vars = collect(AbstractAlgebra.symbols(R))
    #     S_str = symbol(S)
    #     vars_str = join( map( x -> string(x), vars ), ",")
    #     "$(S_str)[$(vars_str)]"
    else
        "R"
    end
end