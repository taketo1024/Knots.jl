using SparseArrays: AbstractSparseMatrix
using AbstractAlgebra

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

function symbol(R::RR) :: String where {RR <: AbstractAlgebra.Ring}
    if isa(R, AbstractAlgebra.Integers)
        "Z"
    elseif isa(R, AbstractAlgebra.Rationals)
        "Q"
    elseif isa(R, AbstractAlgebra.GFField)
        p = characteristic(R)
        "F$(subscript(p))"
    elseif isa(R, AbstractAlgebra.PolyRing) || isa(R, AbstractAlgebra.MPolyRing)
        S = base_ring(R)
        vars = collect(AbstractAlgebra.symbols(R))
        S_str = symbol(S)
        vars_str = join( map( x -> string(x), vars ), ",")
        "$(S_str)[$(vars_str)]"
    else
        "R"
    end
end

function print_matrix(A::AbstractMatrix)
    Base.print_matrix(stdout, A, "[", " ", "]")
end

function print_matrix(A::AbstractSparseMatrix)
    print_matrix(Array(A))
end