using Khovanov.Links: Link
using Khovanov.Jones: kauffmanBracket

l = Link([0, 0, 1, 1])
f = kauffmanBracket(l)
println(f)