using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
m = mesh(1:10)

pb = IntLegendreBasis(m, T; order = 8, left = false, right = false)

X = LinRange(1,10,1000)
plt = plot()
for p ∈ pb
    plot!(X,p.(X))
end
plt

mass_matrix(pb,1,10)