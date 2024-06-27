using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = logmesh(1,5,10)
normalize = true
ordermin = 2
ordermax = 3

basis = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)
@time "Short Integrated Legendre Mass matrix" short_mm = mass_matrix(basis)
display(short_mm)

X = LinRange(1,5,1000)
plt = plot()
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(X, p.(X))
end
plt