using Khovanov: Link
using Khovanov: jonesPolynomial

l = Link([0, 0, 1, 1])
f = jonesPolynomial(l)
println(f)