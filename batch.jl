using Knots.Links, Knots.Khovanov
using Knots.Extensions: GaussInt
using PyCall
using CSV
using DataFrames

ENV["JULIA_DEBUG"] = "Knots"

sg = pyimport("spherogram")

function run(name::String, file="result.csv") :: Bool
    pdcode = try
        sg.Link(name).PD_code()
    catch
        return false
    end

    println("start: $name")

    L = Link(pdcode)

    @time s_2  = s_c(L, 2; reduced=true)
    @time s_3  = s_c(L, 3; reduced=true)
    @time s_1i = s_c(L, GaussInt(1, 1); reduced=true)
    @time s_1i_unred = s_c(L, GaussInt(1, 1); reduced=false)

    nontrivial = !(s_2 == s_3 == s_1i == s_1i_unred)
    
    @show name s_2 s_3 s_1i s_1i_unred nontrivial

    result = DataFrame(
        name = name, 
        s_2 = s_2, 
        s_3 = s_3, 
        s_1i = s_1i,
        s_1i_unred = s_1i_unred,
        nontrivial = nontrivial
    )

    CSV.write(file, result; append = isfile(file))
    true
end

for i in 11:14
    for j in 1:10000
        name = "K$(i)n$(j)"
        next = run(name)
        !next && break 
        break
    end
    break
end


nothing