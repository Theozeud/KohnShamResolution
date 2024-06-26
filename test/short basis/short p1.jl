using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = logmesh(1,10,10)
normalize = false
left = false
right = true

p1 = P1Basis(m, T; left = left, right = right)
@time "Complete P1 Mass matrix" mm = mass_matrix(p1)
display(mm)

shortp1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)
@time "Short P1 Mass matrix" short_mm = mass_matrix(shortp1)
display(short_mm)





