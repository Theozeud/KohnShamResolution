using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = logmesh(1,10,10)
normalize = true
ordermin = 2
ordermax = 4

basis = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)
@time "Short Integrated Legendre Mass matrix" short_mm = mass_matrix(basis)
display(short_mm)
