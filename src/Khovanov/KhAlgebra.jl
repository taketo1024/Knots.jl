using ..Extensions: symbol

# KhAlgGenerator
# Generators of A = R[X]/(X^2 - hX - t) ≅ R<1, X>

@enum KhAlgGenerator I X

function q_degree(x::KhAlgGenerator) 
    (x == I) ? 0 : -2
end

function asString(x::KhAlgGenerator) 
    (x == I) ? "1" : "X"
end

Base.show(io::IO, x::KhAlgGenerator) = 
    print(io, asString(x))

# KhAlgStructure
# Structure of A = R[X]/(X^2 - hX - t)

struct KhAlgStructure{R}
    h::R
    t::R
end

KhAlgStructure(name::String) = begin
    (h, t) = _selectAlgStructure(name)
    KhAlgStructure(h, t)
end

function product(A::KhAlgStructure)
    (x::KhAlgGenerator, y::KhAlgGenerator) -> begin
        (h, t) = (A.h, A.t)
        one = Base.one(typeof(h))

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

function coproduct(A::KhAlgStructure)
    (x::KhAlgGenerator) -> begin
        (R, h, t) = (typeof(A.h), A.h, A.t)
        one = Base.one(R)

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

function asString(A::KhAlgStructure) :: String
    (R, h, t) = (typeof(A.h), A.h, A.t)
    "R = $(symbol(R)), (h, t) = ($h, $t)"
end

using GaloisFields
using Polynomials

function _selectAlgStructure(name::String) 
    Z = Int
    Q = Rational{Int}
    F2 = @GaloisField(2)
    F3 = @GaloisField(3)

    R = if startswith(name, "Z-")
        Z
    elseif startswith(name, "Q-")
        Q
    elseif startswith(name, "F2-")
        F2
    elseif startswith(name, "F3-")
        F3
    elseif startswith(name, "Z[h]-")
        Polynomial{Z, :h}
    elseif startswith(name, "Q[h]-")
        Polynomial{Q, :h}
    elseif startswith(name, "F2[h]-")
        Polynomial{F2, :h}
    elseif startswith(name, "F3[h]-")
        Polynomial{F3, :h}
    elseif startswith(name, "Z[t]-")
        Polynomial{Z, :t}
    elseif startswith(name, "Q[t]-")
        Polynomial{Q, :t}
    elseif startswith(name, "F2[t]-")
        Polynomial{F2, :t}
    elseif startswith(name, "F3[t]-")
        Polynomial{F3, :t}
    else
        Z
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
    else
        throw(Exception)
    end

    (h, t)
end

Base.show(io::IO, A::KhAlgStructure) = 
    print(io, asString(A))
