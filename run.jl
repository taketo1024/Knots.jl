using Knots.Links
using Knots.Khovanov

types = [
    "Z-Kh", 
    "Q-Kh", 
    "F2-Kh", 
    "F3-Kh", 
    "Q[h]-bigraded", 
    "F2[h]-bigraded", 
    "F3[h]-bigraded", 
    "Q[t]-bigraded", 
    "F2[t]-bigraded", 
    "F3[t]-bigraded"
]

for t in types
    l = trefoil
    A = KhAlgStructure(t)
    H = KhHomology(A, l)
    println(H, "\n")
end

nothing