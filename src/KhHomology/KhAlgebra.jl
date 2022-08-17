using AbstractAlgebra
using AbstractAlgebra: Ring
using .Utils

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

struct KhAlgStructure{R <: RingElement, RR <: Ring}
    h::R
    t::R
    baseRing::RR
    R_symbol::String
end

KhAlgStructure(h::R, t::R; R_symbol="R") where {R <: RingElement} = begin
    KhAlgStructure(h, t, parent(h), R_symbol)
end

KhAlgStructure(name::String) = begin
    (h, t, baseRing, R_symbol) = _selectAlgStructure(name)
    KhAlgStructure(h, t, baseRing, R_symbol)
end

function product(A::KhAlgStructure{R}) where {R <: RingElement}
    (x::KhAlgGenerator, y::KhAlgGenerator) -> begin
        (h, t) = (A.h, A.t)
        one = Base.one(A.baseRing)

        # 1)  1^2 = 1
        # 2)  X1 = 1X = X
        # 3)  X^2 = hX + t

        res = if (x, y) == (I, I)         
            [(I, one)]
        elseif  (x, y) ∈ [(I, X), (X, I)] 
            [(X, one)]           
        else
            [(X, h), (I, t)]
        end

        filter(x -> !iszero(x[2]), res)
    end
end

function coproduct(A::KhAlgStructure{R}) where {R <: RingElement} 
    (x::KhAlgGenerator) -> begin
        (h, t) = (A.h, A.t)
        one = Base.one(A.baseRing)

        # 1)  Δ1 = 1X + X1 - h(11)
        # 2)  ΔX = XX + t(11)
        res = if x == I
            [(I, X, one), (X, I, one), (I, I, -h)]
        else
            [(X, X, one), (I, I, t)]
        end
        filter(x -> !iszero(x[3]), res)
    end
end

function asString(A::KhAlgStructure{R}) :: String where {R <: RingElement} 
    "R = $(A.R_symbol), (h, t) = ($(A.h), $(A.t))"
end

function _selectAlgStructure(name::String) 
    Z = parent(1)
    Q = parent(1//1)
    F2 = AbstractAlgebra.GF(2)
    F3 = AbstractAlgebra.GF(3)

    (R, R_symbol) = if startswith(name, "Z-")
        (Z, "Z")
    elseif startswith(name, "Q-")
        (Q, "Q")
    elseif startswith(name, "F2-")
        (F2, "F₂")
    elseif startswith(name, "F3-")
        (F3, "F₃")
    elseif startswith(name, "Z[h]-")
        (PolynomialRing(ZZ, :h)[1], "Z[h]")
    elseif startswith(name, "Q[h]-")
        (PolynomialRing(QQ, :h)[1], "Q[h]")
    elseif startswith(name, "F2[h]-")
        (PolynomialRing(F2, :h)[1], "F₂[h]")
    elseif startswith(name, "F3[h]-")
        (PolynomialRing(F3, :h)[1], "F₃[h]")
    elseif startswith(name, "Z[t]-")
        (PolynomialRing(ZZ, :t)[1], "Z[t]")
    elseif startswith(name, "Q[t]-")
        (PolynomialRing(QQ, :t)[1], "Q[t]")
    elseif startswith(name, "F2[t]-")
        (PolynomialRing(F2, :t)[1], "F₂[t]")
    elseif startswith(name, "F3[t]-")
        (PolynomialRing(F3, :t)[1], "F₃[t]")
    elseif startswith(name, "Z[h,t]-")
        (PolynomialRing(ZZ, [:h, :t])[1], "Z[h,t]")
    elseif startswith(name, "Q[h,t]-")
        (PolynomialRing(QQ, [:h, :t])[1], "Q[h,t]")
    elseif startswith(name, "F2[h,t]-")
        (PolynomialRing(F2, [:h, :t])[1], "F₂[h,t]")
    elseif startswith(name, "F3[h,t]-")
        (PolynomialRing(F3, [:h, :t])[1], "F₃[h,t]")
    else
        (Z, "Z")
    end

    (h, t) = if endswith(name, "Kh")
        (R(0), R(0))
    elseif endswith(name, "BN")
        (R(1), R(0))
    elseif endswith(name, "Lee")
        (R(0), R(1))
    elseif endswith(name, "[h]-bigraded")
        (R([0, 1]), R(0))
    elseif endswith(name, "[t]-bigraded")
        (R(0), R([0, 1]))
    elseif endswith(name, "[h,t]-bigraded")
        (R([1], [[1, 0]]), R([1], [[0, 1]]))
    else
        throw(Exception)
    end

    (h, t, parent(h), R_symbol)
end

Base.show(io::IO, A::KhAlgStructure{R}) where {R} = 
    print(io, asString(A))
