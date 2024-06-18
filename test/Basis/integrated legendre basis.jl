using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Rational{Int}
m = mesh(1:5)
order = 3

pb = IntLegendreBasis(m, T; order = order, left = false, right = false)

X = LinRange(1,5,1000)
plt = plot()
for p ∈ pb
    plot!(X,p.(X))
end
plt

A = mass_matrix(pb,1,5)
