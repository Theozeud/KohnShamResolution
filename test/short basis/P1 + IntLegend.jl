using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
Rmin = 0.00001
Rmax = 100
Nmesh = 40
m = logmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 3

intleg = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)

left = false
right = true

p1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)

intlegp1 = CombineShortPolynomialBasis(p1, intleg)

mass_matrix(intlegp1)