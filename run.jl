using Knots.Links, Knots.Khovanov

ENV["JULIA_DEBUG"] = "Knots"

L = trefoil
c = 2
r = true

@time s = s_c(L, c; reduced=r)
@show s

nothing