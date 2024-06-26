using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = mesh(1:5)
normalize = true

shortp1 = ShortP1Basis(m, T; normalize = normalize, left = true, right = true)
short_mm = mass_matrix(shortp1)

p1 = P1Basis(m, T; left = false, right = false)
mm = mass_matrix(p1)