using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
Rmin = 0.00001
Rmax = 100
Nmesh = 5
m = logmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 3

diffleg = ShortDiffLegendreBasis(m, T; normalize = normalize, ordermax = ordermax)

left = false
right = true

p1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)

difflegp1 = CombineShortPolynomialBasis(p1, diffleg)

mass_matrix(difflegp1)


# Plot of basis
X = LinRange(Rmin,Rmax ,1000)
plt = plot()
for i ∈ 1:length(difflegp1)
    p = build_basis(difflegp1, i)
    plot!(X, p.(X))
end
plt
