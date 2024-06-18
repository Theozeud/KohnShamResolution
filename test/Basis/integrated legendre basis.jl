using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
m = mesh(1:5)
order = 2

pb = IntLegendreBasis(m, T; order = order, left = false, right = false)

X = LinRange(1,5,1000)
plt = plot()
for p ∈ pb
    plot!(X,p.(X))
end
plt


pb_hat = P1Basis(m, T; left = false, right = false)

A = mass_matrix(pb,1,5)
Ahat = mass_matrix(pb_hat,1,5)