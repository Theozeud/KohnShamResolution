using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = logmesh(1,10,5)
normalize = true
ordermin = 2
ordermax = 3

intleg = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)

left = false
right = true

p1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)


intlegp1 = CombineShortPolynomialBasis(intleg, p1)

mass_matrix(intlegp1)


