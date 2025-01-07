using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
Rmin = 0.00001
Rmax = 100
Nmesh = 100
m = logmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 3

basis = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)
@time "Short Integrated Legendre Mass matrix" short_mm = mass_matrix(basis)

# Plot of basis
X = LinRange(1,5,1000)
plt = plot()
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(X, p.(X))
end
plt
