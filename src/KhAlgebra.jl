using AbstractAlgebra

# KhAlgGenerator
# Generators of A = R[X]/(X^2 - hX - t) ≅ R<1, X>

@enum KhAlgGenerator I X

function degree(x::KhAlgGenerator) 
    (x == I) ? 0 : -2
end

function asString(x::KhAlgGenerator) 
    (x == I) ? "1" : "X"
end

Base.show(io::IO, x::KhAlgGenerator) = 
    print(io, asString(x))

# KhAlgStructure
# Structure of A = R[X]/(X^2 - hX - t)

struct KhAlgStructure{R <: RingElement}
    h::R
    t::R
    zero::R
    one::R
end

KhAlgStructure(h::R, t::R) where {R <: RingElement} = begin
    _R = parent(h)
    KhAlgStructure(h, t, zero(_R), one(_R))
end

function product(A::KhAlgStructure{R}) where {R <: RingElement}
    (x::KhAlgGenerator, y::KhAlgGenerator) -> begin
        # 1)  1^2 = 1
        # 2)  X1 = 1X = X
        # 3)  X^2 = hX + t
        res = if (x, y) == (I, I)         
            [(I, A.one)]           
        elseif  (x, y) ∈ [(I, X), (X, I)] 
            [(X, A.one)]           
        else
            [(X, A.h), (I, A.t)]
        end
        filter(x -> !iszero(x[2]), res)
    end
end

function coproduct(A::KhAlgStructure{R}) where {R <: RingElement} 
    (x::KhAlgGenerator) -> begin
        # 1)  Δ1 = 1X + X1 - h(11)
        # 2)  ΔX = XX + t(11)
        res = if x == I
            [(I, X, A.one), (X, I, A.one), (I, I, -A.h)]
        else
            [(X, X, A.one), (I, I, A.t)]
        end
        filter(x -> !iszero(x[3]), res)
    end
end