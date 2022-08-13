# KhAlgGenerator
# Generators of A = R[X]/(X^2 - hX - t) ≅ R<1, X>

@enum KhAlgGenerator I X

degree(x::KhAlgGenerator) = (x == I) ? 0 : -2
asString(x::KhAlgGenerator) = (x == I) ? "1" : "X"
# Base.<(x::KhAlgGenerator, y::KhAlgGenerator) = degree(x) < degree(y)
Base.show(io::IO, x::KhAlgGenerator) = print(io, asString(x))

# KhAlgStructure
# Structure of A = R[X]/(X^2 - hX - t)

struct KhAlgStructure{R}
    h::R
    t::R
end

@enum KhAlgStructureType Kh Lee BN

KhAlgStructure{R}(type::KhAlgStructureType) where {R} = begin
    (type == Kh ) ? KhAlgStructure(zero(R), zero(R)) : 
    (type == Lee) ? KhAlgStructure(zero(R), one(R)) : 
    (type == BN ) ? KhAlgStructure(one(R), zero(R)) : 
    KhAlgStructure(zero(R), zero(R))
end

product(A::KhAlgStructure{R}) where {R} = begin
    (x::KhAlgGenerator, y::KhAlgGenerator) -> begin
        # 1)  1^2 = 1
        # 2)  X1 = 1X = X
        # 3)  X^2 = hX + t
        res = if (x, y) == (I, I)         
            [(I, one(R))]           
        elseif  (x, y) ∈ [(I, X), (X, I)] 
            [(X, one(R))]           
        else
            [(X, A.h), (I, A.t)]
        end
        filter(x -> x[2] != zero(R), res)
    end
end

coproduct(A::KhAlgStructure{R}) where {R} = begin
    (x::KhAlgGenerator) -> begin
        # 1)  Δ1 = 1X + X1 - h(11)
        # 2)  ΔX = XX + t(11)
        res = if x == I
            [(I, X, one(R)), (X, I, one(R)), (I, I, -A.h)]
        else
            [(X, X, one(R)), (I, I, A.t)]
        end
        filter(x -> x[3] != zero(R), res)
    end
end